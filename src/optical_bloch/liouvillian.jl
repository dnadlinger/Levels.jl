# Numerical methods on the assembled master equations, independent of any
# solver package: the Liouvillian of an internal ⊗ Fock model restricted to a
# band of motional coherences, its time evolution, thermal states, and the
# integrated-transient solver for steady states and integral relaxation times.
# Everything here works on plain (unit-stripped) matrices over the internal ⊗
# Fock layout of `MotionalModel`: internal index slowest, one Fock ladder per
# mode, last mode fastest — the `tensor` order of the solvers.

"""
    BandLiouvillian(H, c_ops, num_states, num_fock; bandwidth = nothing)
    BandLiouvillian(mm::MotionalModel; bandwidth = nothing)       (extension)
    BandLiouvillian(model::LindbladModel; time_unit = u"µs")      (extension)

The master-equation generator
``\\mathcal{L}ρ = -i[H, ρ] + \\sum_k (L_k ρ L_k^† - \\tfrac{1}{2}\\{L_k^† L_k, ρ\\})``
of an internal ⊗ Fock model — sparse `H` and jump operators `c_ops` as plain
matrices over `num_states` internal states and one Fock ladder of
`num_fock[m]` states per mode (internal index slowest, last mode fastest) —
restricted to the density-matrix elements ``ρ_{(i\\,\\vec{n}),(j\\,\\vec{m})}``
with ``|n_m - m_m| ≤ \\mathtt{bandwidth}`` for every mode: a Galerkin
truncation of the motional coherences, the full Liouvillian for
`bandwidth = nothing` (or no modes at all, `num_fock = Int[]`, for a purely
internal model).

The far-off-diagonal coherences are both weakly driven — the displacement
matrix elements of [`displacement_elements`](@ref) fall off as
``(η\\sqrt{n})^{|n - m|}/|n - m|!`` — and detuned by ``ω_m |n - m|``, so the
restriction converges quickly (a bandwidth of 4 reproduces steady states and
relaxation times of the quenched-sideband-cooling example to ``10^{-5}``),
while the dimension drops from ``(N_\\mathrm{int} \\prod_m N_m)^2`` to about
``N_\\mathrm{int}^2 \\prod_m (2K + 1) N_m``, which is what makes Fock ladders
of hundreds of states — a Doppler-cooled ion — affordable.

Fields: `L` (sparse, the restricted generator over the band index set),
`index` (full ``(r, c)`` → band index, `0` outside the band), `entries` (band
index → ``(r, c)``), `num_states`, `num_fock` and `bandwidth`.
[`band_vector`](@ref) and [`full_matrix`](@ref) convert between density
matrices and band vectors, [`expectation_row`](@ref) gives the trace functional
of an observable, [`fock_populations`](@ref) the motional populations;
[`cooling_curve`](@ref) integrates the band equation and
[`IntegratedTransientSolver`](@ref) solves it for steady states and
relaxation times. The extension constructors take the operators from a
[`MotionalModel`](@ref) or, for the internal problem alone, a
[`LindbladModel`](@ref).
"""
struct BandLiouvillian
    L::SparseMatrixCSC{ComplexF64,Int}
    index::Matrix{Int}
    entries::Vector{Tuple{Int,Int}}
    num_states::Int
    num_fock::Vector{Int}
    bandwidth::Int
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
)
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

Returns the row vector `a` with `transpose(a) * v == tr(A * ρ)` for the band
vector `v` of any density matrix `ρ` — the trace functional of the observable
`A` on the band.
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
)
    v0 = band_vector(bl, ρ0)
    L = bl.L
    rhs!(dv, v, _, _) = mul!(dv, L, v)
    span = (float(first(ts)), float(last(ts)))
    sol = solve(ODEProblem(rhs!, v0, span), Vern7(); saveat=ts, reltol, abstol)
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
    (; t=collect(sol.t), nbar, fock=populations)
end

"""
    IntegratedTransientSolver(bl::BandLiouvillian, observables; atol = 1e-9)
    IntegratedTransientSolver(mm::MotionalModel; bandwidth = nothing, observables)   (extension)
    IntegratedTransientSolver(model::LindbladModel, observables; time_unit)       (extension)

The integrated-transient Liouvillian method: the time integral of the
transient towards the steady state,

```math
\\hat{τ} = ∫_0^∞ (ρ(t) - ρ_∞) \\, dt, \\qquad
\\mathcal{L} \\hat{τ} = ρ_∞ - ρ_0, \\quad \\mathrm{Tr}\\, \\hat{τ} = 0,
```

solved as a linear system instead of propagating in time, from which the
*integral relaxation time* of any observable ``A`` follows as
``τ_A = \\mathrm{Tr}(A \\hat{τ}) / (\\mathrm{Tr}(A ρ_0) - \\mathrm{Tr}(A ρ_∞))``
(Il'enkov et al., JETP **123**, 1 (2016); Kulosa et al., New J. Phys. **25**,
053008 (2023), Sec. 3, Eqs. (10)–(14), there the "τ-matrix method"; for the
integral relaxation time cf. Garanin, Phys. Rev. E **54**, 3250 (1996)) —
see [`integral_relaxation_time`](@ref).

One sparse LU factorisation of the Liouvillian with one population row
replaced by the trace functional (the rows of ``\\mathcal{L}`` are linearly
dependent through ``\\mathrm{Tr}\\,\\mathcal{L}ρ = 0``) yields the steady state
(forward solve; [`steady_state`](@ref), [`steady_value`](@ref)) and, as
``\\mathrm{Tr}(A \\hat{τ})`` is linear in ``ρ_0``, the adjoint solution
``y = M^{-T} a`` for the trace functional ``a`` of each observable, after
which the relaxation time from *any* initial state costs one dot product.
`observables` is one matrix on the internal ⊗ Fock space, or a vector of them
(the extension defaults to the phonon-number operators of the modes).

Strict partial pivoting is enforced in UMFPACK: its default threshold
pivoting silently returns wrong solutions once the steady state spans ~100
orders of magnitude along a Fock ladder (``N ≳ 170`` at ``η = 0.1``). The
factorisation is released immediately (it lives outside the Julia heap, so
scans with hundreds of factorisations would otherwise exhaust the memory).
Both solve residuals are checked, and the steady state must be a density
matrix to within `atol`; an inaccurate factorisation raises.
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
    num_entries = size(bl.L, 1)
    d = full_dimension(bl)
    diagonal = [bl.index[r, r] for r in 1:d]
    row = diagonal[1]
    # Replace row `row` by the trace functional.
    M = copy(bl.L)
    M[row, :] .= 0
    dropzeros!(M)
    M += sparse(fill(row, d), diagonal, ones(ComplexF64, d), num_entries, num_entries)
    control = SparseArrays.UMFPACK.get_umfpack_control(ComplexF64, Int)
    control[SparseArrays.UMFPACK.JL_UMFPACK_PIVOT_TOLERANCE] = 1.0
    F = lu(M; control)
    e = zeros(ComplexF64, num_entries)
    e[row] = 1
    steady = refined_solve(F, M, e)
    rows = [expectation_row(bl, A) for A in obs]
    Mt = sparse(transpose(M))
    adjoints = [refined_solve(transpose(F), Mt, a) for a in rows]
    finalize(F)
    check_residual(M * steady - e, e)
    for (y, a) in zip(adjoints, rows)
        check_residual(Mt * y - a, a)
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

function check_residual(residual, rhs)
    norm(residual) <= 1e-8 * norm(rhs) || error(
        "Inaccurate Liouvillian solve (relative residual " *
        "$(norm(residual) / norm(rhs)))",
    )
end

"""
    steady_state(s::IntegratedTransientSolver) -> Matrix{ComplexF64}

Returns the steady-state density matrix of the solver's Liouvillian (full
matrix, zero outside the band).
"""
steady_state(s::IntegratedTransientSolver) = full_matrix(s.liouvillian, s.steady)

"""
    steady_value(s::IntegratedTransientSolver, k = 1) -> Float64

Returns the steady-state expectation value ``\\mathrm{Tr}(A_k ρ_∞)`` of the
solver's `k`-th observable.
"""
steady_value(s::IntegratedTransientSolver, k::Integer=1) = s.steady_values[k]

"""
    integral_relaxation_time(s::IntegratedTransientSolver, ρ0, k = 1) -> Float64

Returns the integral relaxation time
``τ_A = ∫_0^∞ (⟨A⟩(t) - ⟨A⟩_∞) \\, dt \\, / \\, (⟨A⟩(0) - ⟨A⟩_∞)`` of the solver's
`k`-th observable from the initial density matrix `ρ0` (full matrix over the
layout), in the time unit of the Liouvillian — the time constant for an
exponential relaxation, and a well-defined "area" time otherwise (for the
phonon number: the cooling time of Kulosa et al. 2023, Eq. (13)). It is
undefined (and returns ``±∞``) when the observable starts at its steady value.
"""
function integral_relaxation_time(
    s::IntegratedTransientSolver,
    ρ0::AbstractMatrix,
    k::Integer=1,
)
    v0 = band_vector(s.liouvillian, ρ0)
    b = s.steady - v0
    b[s.trace_row] = 0   # Tr τ̂ = 0
    real(transpose(s.adjoints[k]) * b) /
    (real(transpose(s.rows[k]) * v0) - s.steady_values[k])
end

export BandLiouvillian, band_vector, full_matrix, expectation_row, fock_populations
export phonon_number_operator, thermal_populations, thermal_state, cooling_curve
export IntegratedTransientSolver, steady_state, steady_value, integral_relaxation_time
