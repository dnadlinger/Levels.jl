# Numerical methods on the assembled master equations, independent of any
# solver package: the Liouvillian of an internal ⊗ Fock model restricted to a
# band of motional coherences, its time evolution, thermal states, the
# integrated-transient solver for steady states and integral relaxation times,
# and the slowest relaxation modes (asymptotic cooling time).
# Everything here works on plain (unit-stripped) matrices over the internal ⊗
# Fock layout of `MotionalModel`: internal index slowest, one Fock ladder per
# mode, last mode fastest — the `tensor` order of the solvers.

"""
    BandLiouvillian(H, c_ops, num_states, num_fock; bandwidth = nothing, time_unit = u"µs")
    BandLiouvillian(mm::MotionalModel; bandwidth = nothing)       (extension)
    BandLiouvillian(model::LindbladModel; time_unit = u"µs")      (extension)

The Lindblad generator on internal ⊗ Fock space, restricted to density-matrix
elements ``ρ_{(i\\vec n),(j\\vec m)}`` with
``|n_ℓ-m_ℓ|≤\\mathtt{bandwidth}`` for every mode ``ℓ``. `bandwidth = nothing`
retains all coherences; `num_fock = Int[]` gives an internal-only model.

Plain `H` and `c_ops` matrices use an internal factor of `num_states`, followed
by Fock factors of sizes `num_fock` (last mode fastest). Their entries use
inverse `time_unit` and inverse square-root `time_unit`, respectively.
The model constructors supply these matrices and their time scale.

For bandwidth ``K`` small compared with the Fock sizes, the generator's
dimension falls from ``(N_\\mathrm{int}∏_ℓN_ℓ)^2`` to about
``N_\\mathrm{int}^2∏_ℓ[(2K+1)N_ℓ]``. This is a projection of the Liouvillian;
increase bandwidth and Fock sizes to check convergence.

`L` stores the sparse generator. `entries` maps band indices to density-matrix
indices `(r, c)`; `index` maps back, with zero outside the band.
[`band_vector`](@ref) and [`full_matrix`](@ref) convert representations.
[`cooling_curve`](@ref) propagates the equation,
[`IntegratedTransientSolver`](@ref) finds steady states and integral
relaxation times, and [`relaxation_modes`](@ref) the slowest eigenmodes.
"""
struct BandLiouvillian{T<:Unitful.Units}
    L::SparseMatrixCSC{ComplexF64,Int}
    index::Matrix{Int}
    entries::Vector{Tuple{Int,Int}}
    num_states::Int
    num_fock::Vector{Int}
    bandwidth::Int
    time_unit::T
end

fock_sizes(num_fock) = num_fock isa Integer ? [Int(num_fock)] : collect(Int, num_fock)

"""Returns the dimension of the internal ⊗ Fock space of the layout."""
full_dimension(bl::BandLiouvillian) = bl.num_states * prod(bl.num_fock)

function BandLiouvillian(
    H::AbstractMatrix,
    c_ops,
    num_states::Integer,
    num_fock;
    bandwidth::Union{Nothing,Integer}=nothing,
    time_unit::Unitful.Units=u"µs",
)
    dimension(time_unit) == Unitful.𝐓 ||
        throw(ArgumentError("time_unit must be a time unit"))
    nf = fock_sizes(num_fock)
    all(>=(1), nf) || throw(ArgumentError("Each Fock ladder needs at least one state"))
    d = num_states * prod(nf)
    size(H) == (d, d) || throw(
        ArgumentError(
            "The Hamiltonian is $(size(H)) but the layout has $d states " *
            "($num_states internal × $(nf) Fock)",
        ),
    )
    K = isnothing(bandwidth) ? (isempty(nf) ? 0 : maximum(nf) - 1) : Int(bandwidth)
    K >= 0 || throw(ArgumentError("The bandwidth must be non-negative, got $K"))
    fock = [fock_numbers(r, nf) for r in 1:d]
    index = zeros(Int, d, d)
    entries = Tuple{Int,Int}[]
    for c in 1:d, r in 1:d
        if all(abs(fock[r][m] - fock[c][m]) <= K for m in eachindex(nf))
            push!(entries, (r, c))
            index[r, c] = length(entries)
        end
    end
    BandLiouvillian(
        band_generator(H, c_ops, index, entries),
        index,
        entries,
        num_states,
        nf,
        K,
        time_unit,
    )
end

# Assembles the generator over the band index set directly from the sparse
# operators: every term is of the form ρ ↦ A ρ B, whose matrix element from
# source (r′, c′) to target (r, c) is A[r, r′] B[c′, c].
function band_generator(H, c_ops, index, entries)
    d = size(index, 1)
    num_entries = length(entries)
    S = SparseMatrixCSC{ComplexF64,Int}
    Hs = S(H)
    Ls = [S(L) for L in c_ops]
    Id = sparse(ComplexF64(1) * I, d, d)
    G = isempty(Ls) ? spzeros(ComplexF64, d, d) : sum(L' * L for L in Ls)
    terms = Tuple{S,S}[(-im * Hs, Id), (Id, im * Hs), (-G / 2, Id), (Id, -G / 2)]
    for L in Ls
        push!(terms, (L, sparse(L')))
    end
    total = spzeros(ComplexF64, num_entries, num_entries)
    for (A, B) in terms
        Bt = sparse(transpose(B))   # column c′ of Bᵀ holds row c′ of B
        Is, Js, Vs = Int[], Int[], ComplexF64[]
        for (src, (r′, c′)) in enumerate(entries)
            for p in nzrange(A, r′)
                r = rowvals(A)[p]
                a = nonzeros(A)[p]
                for q in nzrange(Bt, c′)
                    tgt = index[r, rowvals(Bt)[q]]
                    tgt == 0 && continue
                    push!(Is, tgt)
                    push!(Js, src)
                    push!(Vs, a * nonzeros(Bt)[q])
                end
            end
        end
        total += sparse(Is, Js, Vs, num_entries, num_entries)
    end
    total
end

"""
    band_vector(bl::BandLiouvillian, ρ::AbstractMatrix) -> Vector{ComplexF64}

Returns the band vector of a full density matrix over the layout of `bl`
(elements outside the band are dropped).
"""
band_vector(bl::BandLiouvillian, ρ::AbstractMatrix) =
    ComplexF64[ρ[r, c] for (r, c) in bl.entries]

"""
    full_matrix(bl::BandLiouvillian, v) -> Matrix{ComplexF64}

Returns the full density matrix of a band vector (zero outside the band).
"""
function full_matrix(bl::BandLiouvillian, v)
    d = full_dimension(bl)
    ρ = zeros(ComplexF64, d, d)
    for (k, (r, c)) in enumerate(bl.entries)
        ρ[r, c] = v[k]
    end
    ρ
end

"""
    expectation_row(bl::BandLiouvillian, A::AbstractMatrix) -> Vector{ComplexF64}

Returns the trace functional `a` with
`transpose(a) * v == tr(A * full_matrix(bl, v))`.
For an arbitrary density matrix, contributions from coherences outside the
band are omitted.
"""
expectation_row(bl::BandLiouvillian, A::AbstractMatrix) =
    ComplexF64[A[c, r] for (r, c) in bl.entries]

"""
    fock_populations(bl::BandLiouvillian, v) -> Vector{Vector{Float64}}
    fock_populations(bl::BandLiouvillian, v, mode) -> Vector{Float64}
    fock_populations(ρ, mm::MotionalModel) -> Vector{Vector{Float64}}   (extension)

Returns the motional populations ``p_n``, ``n = 0 … \\mathtt{num\\_fock} - 1``,
of each mode (or of the given one), traced over the internal states and the
other modes — of a band vector, or in the extension of a QuantumToolbox
density matrix on the space of a [`MotionalModel`](@ref).
"""
function fock_populations(bl::BandLiouvillian, v, mode::Integer)
    p = zeros(bl.num_fock[mode])
    for r in 1:full_dimension(bl)
        p[fock_numbers(r, bl.num_fock)[mode]+1] += real(v[bl.index[r, r]])
    end
    p
end

fock_populations(bl::BandLiouvillian, v) =
    [fock_populations(bl, v, m) for m in eachindex(bl.num_fock)]

"""
    phonon_number_operator(num_states, num_fock[, mode]) -> Diagonal

Returns the phonon-number operator ``1 ⊗ a_m^† a_m`` of the given mode — or,
without one, the total phonon number of all modes — on the internal ⊗ Fock
layout, as a plain diagonal matrix.
"""
function phonon_number_operator(num_states::Integer, num_fock, mode=nothing)
    nf = fock_sizes(num_fock)
    isnothing(mode) || 1 <= mode <= length(nf) || throw(ArgumentError("No mode $mode"))
    d = num_states * prod(nf)
    Diagonal([
        begin
            ns = fock_numbers(r, nf)
            Float64(isnothing(mode) ? sum(ns; init=0) : ns[mode])
        end for r in 1:d
    ])
end

"""
    thermal_populations(nbar, num_fock) -> Vector{Float64}

Returns the Boltzmann populations ``p_n ∝ (\\bar n / (\\bar n + 1))^n`` of a
thermal state of mean occupation `nbar`, truncated to `num_fock` Fock states
and renormalised (cf. [`fock_truncation`](@ref) for the truncation a given
tail population needs).
"""
function thermal_populations(nbar::Real, num_fock::Integer)
    nbar >= 0 || throw(ArgumentError("nbar must be non-negative, got $nbar"))
    num_fock >= 1 || throw(ArgumentError("num_fock must be positive, got $num_fock"))
    p = [(nbar / (nbar + 1))^n for n in 0:(num_fock-1)]
    p ./ sum(p)
end

"""
    thermal_state(internal::AbstractMatrix, num_fock, nbars) -> Matrix{ComplexF64}
    thermal_state(mm::MotionalModel, internal, nbars)                (extension)

Returns the product state ``ρ_\\mathrm{int} ⊗ ρ_\\mathrm{th}(\\bar n_1) ⊗ …`` of an
internal density matrix and truncated thermal states
([`thermal_populations`](@ref)) of the given mean occupations (one per mode,
or one number for all) on the internal ⊗ Fock layout — the usual initial
state of a cooling calculation. The extension form takes the internal state as
a QuantumToolbox ket or density matrix and returns a `QuantumObject` on the
space of the [`MotionalModel`](@ref).
"""
function thermal_state(internal::AbstractMatrix, num_fock, nbars)
    nf = fock_sizes(num_fock)
    nb = nbars isa Real ? fill(Float64(nbars), length(nf)) : collect(Float64, nbars)
    length(nb) == length(nf) ||
        throw(ArgumentError("One mean occupation per mode is required"))
    ρ = Matrix{ComplexF64}(internal)
    for (n, nbar) in zip(nf, nb)
        ρ = kron(ρ, Diagonal(thermal_populations(nbar, n)))
    end
    ρ
end

"""
    cooling_curve(bl::BandLiouvillian, ρ0, ts; fock = false, reltol = 1e-8, abstol = 1e-10)
    cooling_curve(mm::MotionalModel, ρ0, ts; fock = false, kwargs...)   (extension)

Integrates the master equation from the density matrix `ρ0` and returns
`(; t, nbar, fock)`: the mean phonon number of every mode at the times `ts`
(a matrix, one column per mode) and, with `fock = true`, the motional
populations of each mode along the trajectory (one matrix per mode, Fock
number × time; `nothing` otherwise). The band form integrates the band vector
with a Verner 7(6) method (OrdinaryDiffEqVerner); the times are in the unit of
the generator (µs for models from [`motional_model`](@ref)). The extension
form propagates the full master equation of a [`MotionalModel`](@ref) with
QuantumToolbox's `mesolve` (further keywords pass through), with `ts` either
plain numbers in the model's time unit or unitful times.
"""
function cooling_curve(
    bl::BandLiouvillian,
    ρ0::AbstractMatrix,
    ts;
    fock::Bool=false,
    reltol=1e-8,
    abstol=1e-10,
    alg=Vern7(),
    kwargs...,
)
    v0 = band_vector(bl, ρ0)
    L = bl.L
    rhs!(dv, v, _, _) = mul!(dv, L, v)
    times = eltype(ts) <: Quantity ? ustrip.(bl.time_unit, ts) : collect(float.(ts))
    span = (first(times), last(times))
    sol =
        solve(ODEProblem(rhs!, v0, span), alg; saveat=times, reltol, abstol, kwargs...)
    successful_retcode(sol) ||
        error("Master-equation propagation failed: $(sol.retcode)")
    num_modes = length(bl.num_fock)
    rows = [
        expectation_row(bl, phonon_number_operator(bl.num_states, bl.num_fock, m))
        for m in 1:num_modes
    ]
    nbar = [real(transpose(rows[m]) * u) for u in sol.u, m in 1:num_modes]
    populations =
        fock ?
        [
            reduce(hcat, [fock_populations(bl, u, m) for u in sol.u]) for
            m in 1:num_modes
        ] : nothing
    (; t=collect(sol.t) .* bl.time_unit, nbar, fock=populations)
end

# The trace-row system of a band Liouvillian: `M` is `bl.L` with the population
# row `row` (the first diagonal entry) replaced by the trace functional, `F` its
# sparse LU with strict partial pivoting (the default threshold pivoting fails
# silently once the steady state spans ~100 decades), `diagonal` the band
# indices of the populations. The unit vector `e_row` on the right-hand side
# gives the steady state. A traceless right-hand side `b` with `b[row] = 0`
# gives the traceless solution of `L x = b`: as `L` preserves the trace, its
# row `row` is minus the sum of the other population rows, so that equation
# holds there too — `F` applies the group inverse of `L`.
function trace_row_system(bl::BandLiouvillian)
    num_entries = size(bl.L, 1)
    d = full_dimension(bl)
    diagonal = [bl.index[r, r] for r in 1:d]
    row = diagonal[1]
    M = copy(bl.L)
    M[row, :] .= 0
    dropzeros!(M)
    M += sparse(fill(row, d), diagonal, ones(ComplexF64, d), num_entries, num_entries)
    control = SparseArrays.UMFPACK.get_umfpack_control(ComplexF64, Int)
    control[SparseArrays.UMFPACK.JL_UMFPACK_PIVOT_TOLERANCE] = 1.0
    (; M, F=lu(M; control), row, diagonal)
end

# Removes the trace of the band vector `v` (populations at `diagonal`) in place.
function remove_trace!(v::AbstractVector, diagonal)
    t = sum(k -> v[k], diagonal) / length(diagonal)
    for k in diagonal
        v[k] -= t
    end
    v
end

"""
    IntegratedTransientSolver(bl::BandLiouvillian, observables; atol = 1e-9)
    IntegratedTransientSolver(mm::MotionalModel; bandwidth = nothing, observables)   (extension)
    IntegratedTransientSolver(model::LindbladModel, observables; time_unit)       (extension)

Solves for a unique steady state and integrated transients of a relaxing
Liouvillian. The transient integral satisfies

```math
\\hat τ = ∫_0^∞[ρ(t)-ρ_∞]\\,dt, \\qquad
\\mathcal L\\hat τ = ρ_∞-ρ_0, \\qquad \\mathrm{Tr}\\hat τ=0.
```

For an observable ``A``, the integral relaxation time is
``τ_A = \\mathrm{Tr}(A\\hat τ)/[\\mathrm{Tr}(Aρ_0)-\\mathrm{Tr}(Aρ_∞)]``.
This is the τ-matrix method of Il'enkov et al., JETP **123**, 1 (2016), and
Kulosa et al., New J. Phys. **25**, 053008 (2023), Sec. 3, Eqs. (10)–(14);
see also [`integral_relaxation_time`](@ref).

One sparse LU factorisation replaces a population row by the trace constraint.
A forward solve gives the steady state; a transposed solve per observable
stores the linear functional needed for its transient integral. Later
initial states require only dot products. `observables` is one matrix or a
vector of matrices; the motional-model form defaults to phonon numbers.

Strict partial pivoting and iterative refinement handle poorly scaled Fock
ladders. Solves are checked by normwise backward error; steady populations
and trace are checked to `atol`. The factorisation is released after use.
"""
struct IntegratedTransientSolver
    liouvillian::BandLiouvillian
    steady::Vector{ComplexF64}
    rows::Vector{Vector{ComplexF64}}
    steady_values::Vector{Float64}
    adjoints::Vector{Vector{ComplexF64}}
    trace_row::Int
end

function IntegratedTransientSolver(bl::BandLiouvillian, observables; atol::Real=1e-9)
    obs = observables isa AbstractVector ? observables : [observables]
    isempty(obs) && throw(ArgumentError("At least one observable is required"))
    (; M, F, row, diagonal) = trace_row_system(bl)
    e = zeros(ComplexF64, size(M, 1))
    e[row] = 1
    steady = refined_solve(F, M, e)
    rows = [expectation_row(bl, A) for A in obs]
    Mt = sparse(transpose(M))
    adjoints = [refined_solve(transpose(F), Mt, a) for a in rows]
    finalize(F)
    check_solution(M, steady, e)
    for (y, a) in zip(adjoints, rows)
        check_solution(Mt, y, a)
    end
    populations = [steady[k] for k in diagonal]
    if minimum(real, populations) < -atol || abs(sum(populations) - 1) > atol
        error(
            "The steady state is not a density matrix (smallest population " *
            "$(minimum(real, populations)), trace $(sum(populations))); the " *
            "Liouvillian solve is inaccurate",
        )
    end
    values = [real(transpose(a) * steady) for a in rows]
    IntegratedTransientSolver(bl, steady, rows, values, adjoints, row)
end

# Solves M x = b with the factorisation F of M, followed by up to two steps of
# iterative refinement (the Fock ladder makes the system badly scaled).
function refined_solve(F, M, b)
    x = F \ b
    for _ in 1:2
        r = b - M * x
        norm(r) <= 1e-12 * norm(b) && break
        x += F \ r
    end
    x
end

# Checks normwise backward error, allowing for roundoff on the scale ‖M‖∞ ‖x‖∞
# when slow relaxation makes the adjoint solution large. Small backward error
# does not guarantee small forward error for an ill-conditioned system.
function check_solution(M, x, b)
    η = norm(M * x - b, Inf) / (opnorm(M, Inf) * norm(x, Inf) + norm(b, Inf))
    η <= 1e-12 || error("Inaccurate Liouvillian solve (backward error $η)")
end

"""
    steady_density_matrix(s::IntegratedTransientSolver) -> Matrix{ComplexF64}

Returns the steady-state density matrix of the solver's Liouvillian (full
matrix, zero outside the band).
"""
steady_density_matrix(s::IntegratedTransientSolver) =
    full_matrix(s.liouvillian, s.steady)

"""
    steady_value(s::IntegratedTransientSolver, k = 1) -> Float64

Returns the steady-state expectation value ``\\mathrm{Tr}(A_k ρ_∞)`` of the
solver's `k`-th observable.
"""
steady_value(s::IntegratedTransientSolver, k::Integer=1) = s.steady_values[k]

"""
    integral_relaxation_time(s::IntegratedTransientSolver, ρ0, k = 1) -> Unitful.Time

Returns the integral relaxation time
``τ_A = ∫_0^∞ (⟨A⟩(t) - ⟨A⟩_∞) \\, dt \\, / \\, (⟨A⟩(0) - ⟨A⟩_∞)`` of the solver's
`k`-th observable from the initial density matrix `ρ0` (full matrix over the
layout), as a unitful time — the time constant for an
exponential relaxation, and a well-defined "area" time otherwise (for the
phonon number: the cooling time of Kulosa et al. 2023, Eq. (13)). It is
undefined when the observable starts at its steady value; division by zero
then produces `Inf` or `NaN`.
"""
function integral_relaxation_time(
    s::IntegratedTransientSolver,
    ρ0::AbstractMatrix,
    k::Integer=1,
)
    v0 = band_vector(s.liouvillian, ρ0)
    b = s.steady - v0
    b[s.trace_row] = 0   # Tr τ̂ = 0
    (
        real(transpose(s.adjoints[k]) * b) /
        (real(transpose(s.rows[k]) * v0) - s.steady_values[k])
    ) * s.liouvillian.time_unit
end

# The group inverse of a band Liouvillian on its traceless subspace, applied
# through the trace-row factorisation: `y = L⁺ x` solves `L y = P x`, `Tr y =
# 0`, with `P` the projection removing the trace (which the Krylov vectors
# only re-acquire through roundoff). Its eigenvalues are the reciprocals of
# the nonzero eigenvalues of `L`, so the slowest relaxation modes dominate the
# iteration, while the stationary state — whose shift-inverted image at zero
# would swamp them — is excluded exactly.
struct GroupInverse{F}
    factorisation::F
    row::Int
    diagonal::Vector{Int}
    buffer::Vector{ComplexF64}
end

Base.size(G::GroupInverse) = (length(G.buffer), length(G.buffer))
Base.size(G::GroupInverse, i::Integer) = i <= 2 ? length(G.buffer) : 1
Base.eltype(::GroupInverse) = ComplexF64

function LinearAlgebra.mul!(y::AbstractVector, G::GroupInverse, x::AbstractVector)
    b = remove_trace!(copyto!(G.buffer, x), G.diagonal)
    b[G.row] = 0
    ldiv!(y, G.factorisation, b)
end

"""
    relaxation_modes(bl::BandLiouvillian; num_eigenpairs = 12, dense_limit = 2500) -> (; values, vectors, steady)

Returns the `num_eigenpairs` eigenmodes nearest the Liouvillian's origin: its
nonzero eigenvalues ``λ`` of smallest magnitude, in order of increasing
magnitude, as complex rates in inverse `time_unit` (``-\\mathrm{Re}\\,λ`` is
the decay rate, ``\\mathrm{Im}\\,λ`` the oscillation frequency; a positive real
part marks a spurious growing mode of a truncated model), with the
corresponding eigenvectors as unit-norm band-vector columns and the
stationary state as a unit-trace band vector. The modes are traceless: they
describe how an observable's deviation from its steady value decays.
Ordering by ``|λ|`` need not order the decay rates ``-\\mathrm{Re}\\,λ`` when
the modes oscillate.

Generators of dimension up to `dense_limit` are diagonalised densely. Larger
ones are treated by Krylov–Schur iteration (ArnoldiMethod.jl) on the group
inverse of ``\\mathcal L`` — its inverse on the traceless subspace, applied
through the trace-row factorisation of [`IntegratedTransientSolver`](@ref) —
whose dominant eigenvalues are the reciprocals of the slowest nonzero rates.
Unlike shift-and-invert at zero, this excludes the stationary state exactly
instead of letting its numerically singular image swamp the Krylov space,
and it needs no shift on the scale of the unknown relaxation rate. The
iteration starts from a fixed vector and is deterministic, converging the
reciprocal rates to a relative residual of ``10^{-12}``. Eigenvalue errors
can be amplified by their condition numbers. A unique steady state is required.
"""
function relaxation_modes(
    bl::BandLiouvillian;
    num_eigenpairs::Integer=12,
    dense_limit::Integer=2500,
)
    num_entries = size(bl.L, 1)
    1 <= num_eigenpairs <= num_entries - 1 || throw(
        ArgumentError(
            "num_eigenpairs must be between 1 and the Liouvillian dimension minus one",
        ),
    )
    # The stationary state comes from the trace-row solve in either case (an
    # eigenvector of the numerically zero eigenvalue is less accurate).
    (; M, F, row, diagonal) = trace_row_system(bl)
    e = zeros(ComplexF64, num_entries)
    e[row] = 1
    steady = refined_solve(F, M, e)
    check_solution(M, steady, e)
    if num_entries <= dense_limit
        finalize(F)
        eig = eigen(Matrix(bl.L))
        order = sortperm(abs.(eig.values))
        selected = order[2:(num_eigenpairs+1)]
        values = eig.values[selected]
        vectors = eig.vectors[:, selected]
    else
        G = GroupInverse(F, row, diagonal, zeros(ComplexF64, num_entries))
        # A fixed, structureless starting vector (golden-angle phases) makes the
        # iteration deterministic.
        v1 = ComplexF64[
            cis(2.399963229728653 * k) * (1 + 0.5 * cos(0.1 * k)) for k in 1:num_entries
        ]
        remove_trace!(v1, diagonal)
        # Residual tolerance relative to the reciprocal rate; the eigenvalue
        # error is that times the mode's condition number, which grows along
        # the Fock ladder (the slowest modes are the best conditioned).
        decomposition, history =
            partialschur(G; v1, nev=num_eigenpairs, which=:LM, tol=1e-12)
        finalize(F)
        history.converged || error(
            "The Krylov iteration for the slowest Liouvillian modes did not " *
            "converge ($(history.nconverged) of $num_eigenpairs eigenpairs after " *
            "$(history.mvproducts) products); reduce num_eigenpairs",
        )
        reciprocals, vectors = partialeigen(decomposition)
        values = 1 ./ reciprocals
        order = sortperm(abs.(values))
        values = values[order]
        vectors = vectors[:, order]
    end
    for v in eachcol(vectors)
        v ./= norm(v)
    end
    (; values=values .* bl.time_unit^-1, vectors, steady)
end

export BandLiouvillian, band_vector, full_matrix, expectation_row, fock_populations
export phonon_number_operator, thermal_populations, thermal_state, cooling_curve
export IntegratedTransientSolver,
    steady_density_matrix, steady_value, integral_relaxation_time
export relaxation_modes
