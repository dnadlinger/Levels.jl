# Conversions from Levels.jl objects to QuantumToolbox.jl `QuantumObject`s, so
# that the QuantumToolbox solvers (`mesolve`, `steadystate_fourier`, …) can be
# applied to models built from Levels data, and the QuantumToolbox-based
# solvers of the Levels.OpticalBloch layer.
#
# QuantumToolbox works on dimensionless matrices, so the conversions fix a
# reference time unit (the `time_unit` keyword; µs by default, matching the
# µs⁻¹ normalisation of the dynamics layer): a unitful quantity of pure time
# dimension 𝐓ᵖ is stripped in units of `time_unit`ᵖ (µs⁻¹ for Hamiltonians in
# angular units, µs^(-1/2) for Lindblad jump operators, …), and times passed
# to the solvers are then in `time_unit`.

module LevelsQuantumToolboxExt

using LinearAlgebra: Diagonal, I, diag, eigen, norm, tr
using Unitful
using Levels: Levels, StateBasis, stateindex
using Levels.PeriodicDynamics: PeriodicDynamics, DrivenTransition
using Levels.OpticalBloch:
    OpticalBloch,
    AdiabaticElimination,
    BandLiouvillian,
    CoolingMetrics,
    DecayLabel,
    DephasingLabel,
    HeatingLabel,
    IntegratedTransient,
    IntegratedTransientSolver,
    LindbladModel,
    MotionalCoupling,
    MotionalModel,
    RecoilLabel,
    LiouvillianSpectrum,
    emission_rule,
    fock_truncation,
    integral_relaxation_time,
    joint_displacement,
    level_projector,
    phonon_number_operator,
    steady_value,
    thermal_populations,
    validate_lamb_dicke_order
using QuantumToolbox:
    QuantumToolbox,
    Ket,
    Operator,
    QuantumObject,
    SteadyStateLinearSolver,
    destroy,
    eigsolve,
    isket,
    ket2dm,
    liouvillian,
    mat2vec,
    mesolve,
    num,
    ptrace,
    qeye,
    steadystate,
    steadystate_fourier,
    tensor,
    to_sparse,
    vec2mat

# Returns the time power p of a quantity of pure time dimension 𝐓ᵖ (0 for a
# dimensionless one), used to fix the stripping unit `time_unit`ᵖ.
function time_power(x)
    dims = typeof(dimension(x)).parameters[1]
    if isempty(dims)
        return 0
    end
    if !(length(dims) == 1 && dims[1] isa Unitful.Dimension{:Time})
        throw(
            ArgumentError(
                "Only quantities of pure time dimension (rates, angular " *
                "frequencies, √rate jump operators, …) map onto the " *
                "dimensionless time-base convention; got $(dimension(x))",
            ),
        )
    end
    dims[1].power
end

# Validates that the reference time unit is indeed one (e.g. `u"µs"`).
function check_time_unit(time_unit)
    if dimension(time_unit) != Unitful.𝐓
        throw(
            ArgumentError(
                "The reference time unit must be a plain unit " *
                "of time; got $time_unit",
            ),
        )
    end
end

strip_time_units(A::AbstractArray{<:Number}, _) = A
strip_time_units(A::AbstractArray{<:Quantity}, time_unit) =
    ustrip.(time_unit^time_power(first(A)), A)

"""
    QuantumObject(A::AbstractMatrix, basis::StateBasis; time_unit = u"µs")
    QuantumObject(v::AbstractVector, basis::StateBasis; time_unit = u"µs")

Builds the `QuantumToolbox.QuantumObject` operator (or ket) representing the
given matrix (or coefficient vector) over the states of `basis`, in
[`Levels.stateindex`](@ref) order.

Unitful elements of pure time dimension ``𝐓^p`` are stripped in units of
``\\mathtt{time\\_unit}^p`` — with the µs default, µs⁻¹ for Hamiltonians in
the angular convention and µs^(-1/2) for Lindblad jump operators — so times
passed to the QuantumToolbox solvers are in `time_unit`.
"""
function QuantumToolbox.QuantumObject(
    A::AbstractMatrix,
    basis::StateBasis;
    time_unit::Unitful.Units=Unitful.µs,
)
    check_time_unit(time_unit)
    if size(A) != (length(basis), length(basis))
        throw(
            ArgumentError(
                "Matrix size $(size(A)) does not match the basis " *
                "($(length(basis)) states)",
            ),
        )
    end
    QuantumObject(strip_time_units(A, time_unit), Operator(), length(basis))
end

function QuantumToolbox.QuantumObject(
    v::AbstractVector,
    basis::StateBasis;
    time_unit::Unitful.Units=Unitful.µs,
)
    check_time_unit(time_unit)
    if length(v) != length(basis)
        throw(
            ArgumentError(
                "Vector length $(length(v)) does not match the basis " *
                "($(length(basis)) states)",
            ),
        )
    end
    QuantumObject(strip_time_units(v, time_unit), Ket(), length(basis))
end

"""
    basis(b::StateBasis, state)
    basis(b::StateBasis, level, m)

Returns the ket `QuantumObject` for the given state of the [`StateBasis`](@ref)
(cf. [`Levels.stateindex`](@ref) for the accepted state forms).
"""
QuantumToolbox.fock(b::StateBasis, state) =
    QuantumToolbox.fock(length(b), stateindex(b, state) - 1)
QuantumToolbox.fock(b::StateBasis, level, m) =
    QuantumToolbox.fock(length(b), stateindex(b, level, m) - 1)

"""
    projection(b::StateBasis, ket[, bra])

Returns the projection operator ``|ket⟩⟨bra|`` between states of the
[`StateBasis`](@ref) (cf. [`Levels.stateindex`](@ref) for the accepted state
forms); `bra` defaults to `ket`, giving the population projector.
"""
QuantumToolbox.projection(b::StateBasis, ket, bra) =
    QuantumToolbox.projection(length(b), stateindex(b, ket) - 1, stateindex(b, bra) - 1)
QuantumToolbox.projection(b::StateBasis, ket) = QuantumToolbox.projection(b, ket, ket)

function PeriodicDynamics.fourier_hamiltonians(
    dt::DrivenTransition;
    δ=zero(dt.drive_frequency),
    time_unit::Unitful.Units=Unitful.µs,
)
    H_0 = (dt.coupling .+ dt.coupling') ./ 2 .+ Diagonal(dt.frame)
    for i in dt.upper_range
        H_0[i, i] -= δ
    end
    H_p = zero(dt.coupling)
    for drive in dt.drives
        H_p .+= (0.5 * cis(drive.phase)) .* drive.amplitude
    end
    (
        H_0=QuantumObject(H_0, dt.basis; time_unit),
        H_p=QuantumObject(H_p, dt.basis; time_unit),
        H_m=QuantumObject(Matrix(H_p'), dt.basis; time_unit),
        ωd=ustrip(time_unit^-1, dt.drive_frequency),
    )
end


# --- Levels.OpticalBloch: internal-state solvers --------------------------------

hamiltonian_qobj(model::LindbladModel; time_unit) =
    QuantumObject(model.hamiltonian, model.basis; time_unit)
jump_qobjs(model::LindbladModel; time_unit) =
    [QuantumObject(L, model.basis; time_unit) for L in model.jump_operators]

function require_static(model::LindbladModel, what)
    isempty(model.harmonics) || throw(
        ArgumentError(
            "$what requires a static model, but this one carries " *
            "$(length(model.harmonics)) beat-note harmonic(s); only " *
            "steadystate(model) handles those (via steadystate_fourier)",
        ),
    )
end

"""
    liouvillian(model::LindbladModel; time_unit = u"µs")
    liouvillian(mm::MotionalModel)

Returns the Liouvillian superoperator of the internal-state master equation of
the [`Levels.OpticalBloch.LindbladModel`](@ref) (static models only), in
inverse `time_unit`, or of the full internal ⊗ motional master equation of a
[`Levels.OpticalBloch.MotionalModel`](@ref).
"""
function QuantumToolbox.liouvillian(
    model::LindbladModel;
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "liouvillian")
    liouvillian(hamiltonian_qobj(model; time_unit), jump_qobjs(model; time_unit))
end

QuantumToolbox.liouvillian(mm::MotionalModel) = liouvillian(mm.H, mm.c_ops)

"""
    steadystate(model::LindbladModel; time_unit = u"µs", num_harmonics = 4, solver, kwargs...)
    steadystate(mm::MotionalModel; solver, kwargs...)

Returns the steady-state internal density matrix of the
[`Levels.OpticalBloch.LindbladModel`](@ref), or the steady state of the full
internal ⊗ motional master equation of a
[`Levels.OpticalBloch.MotionalModel`](@ref). A static model is solved directly
(`QuantumToolbox.steadystate`); an internal model with beat-note harmonics at
one frequency is solved for its periodic steady state with
`steadystate_fourier` (`num_harmonics` harmonics of ρ, cf. its `n_max`), of
which the period-averaged k = 0 component is returned. `solver` defaults to
the direct sparse factorisation `SteadyStateLinearSolver(; alg = nothing)`;
further keywords pass through.
"""
function QuantumToolbox.steadystate(
    model::LindbladModel;
    time_unit::Unitful.Units=Unitful.µs,
    num_harmonics::Int=4,
    solver=SteadyStateLinearSolver(; alg=nothing),
    kwargs...,
)
    H_0 = hamiltonian_qobj(model; time_unit)
    c_ops = jump_qobjs(model; time_unit)
    if isempty(model.harmonics)
        return steadystate(H_0, c_ops; solver, kwargs...)
    end
    length(model.harmonics) == 1 || throw(
        ArgumentError(
            "The periodic steady state is only available for beat notes at a " *
            "single frequency; this model has $(length(model.harmonics))",
        ),
    )
    w, M = model.harmonics[1]
    # Model convention H(t) = H_0 + M e^{-iwt} + M† e^{+iwt} against QuantumToolbox's
    # H_0 + H_p e^{+iωt} + H_m e^{-iωt}.
    H_m = QuantumObject(M, model.basis; time_unit)
    H_p = QuantumObject(Matrix(M'), model.basis; time_unit)
    ρs = steadystate_fourier(
        H_0,
        H_p,
        H_m,
        ustrip(time_unit^-1, w),
        c_ops;
        n_max=num_harmonics,
        solver,
        kwargs...,
    )
    ρs[1]
end

QuantumToolbox.steadystate(
    mm::MotionalModel;
    solver=SteadyStateLinearSolver(; alg=nothing),
    kwargs...,
) = steadystate(mm.H, mm.c_ops; solver, kwargs...)

function OpticalBloch.populations(ρ::QuantumObject, model::LindbladModel)
    num_states = length(model.basis)
    size(ρ.data) == (num_states, num_states) || throw(
        ArgumentError(
            "Density matrix of size $(size(ρ.data)) does not match the model's " *
            "$num_states internal states (reduce a full motional state with " *
            "ptrace first)",
        ),
    )
    real.(diag(ρ.data))
end

function OpticalBloch.populations(ρ::QuantumObject, model::LindbladModel, level)
    P = level_projector(model, level)
    p = OpticalBloch.populations(ρ, model)
    sum(p[i] for i in 1:length(p) if !iszero(P[i, i]); init=0.0)
end

# --- Levels.OpticalBloch: adiabatic-elimination rates ----------------------------

# Adiabatic-elimination rates of one mode from the stripped internal Liouvillian
# `L` and steady state `ρ` (plain matrices, inverse time_unit).
function mode_rates(L, ρ, mc::MotionalCoupling, time_unit)
    H_sb = ustrip.(time_unit^-1, mc.sideband_hamiltonian)
    ω_m = ustrip(time_unit^-1, mc.mode.frequency)
    y = mat2vec(H_sb * ρ)
    spectrum(ω) = tr(H_sb * vec2mat((-L + (im * ω) * I) \ y))
    diffusion = 0.0
    for R in mc.recoil_operators
        Rs = ustrip.(time_unit^(-1 // 2), R)
        diffusion += real(tr(Rs' * Rs * ρ))
    end
    A_plus = 2 * real(spectrum(ω_m)) + diffusion
    A_minus = 2 * real(spectrum(-ω_m)) + diffusion
    A_plus, A_minus
end

function check_heating_rate(heating_rate)
    (heating_rate isa Unitful.Frequency && heating_rate >= zero(heating_rate)) || throw(
        ArgumentError(
            "heating_rate must be a non-negative rate in quanta per time " *
            "(e.g. 1000.0u\"s^-1\"), got $heating_rate",
        ),
    )
end

# One heating rate per mode, from a single rate or one per mode.
function mode_heating_rates(heating_rate, num_modes)
    rates =
        heating_rate isa AbstractVector ? collect(heating_rate) :
        fill(heating_rate, num_modes)
    length(rates) == num_modes ||
        throw(ArgumentError("heating_rate must be one rate or one per mode"))
    foreach(check_heating_rate, rates)
    rates
end

function OpticalBloch.cooling_rates(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    ρ=steadystate(model),
    heating_rate=0.0u"s^-1",
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "cooling_rates")
    check_heating_rate(heating_rate)
    L = liouvillian(model; time_unit).data
    ρm = Matrix(ρ.data)
    h = ustrip(time_unit^-1, heating_rate)
    map(mcs) do mc
        A_plus, A_minus = mode_rates(L, ρm, mc, time_unit)
        net = A_minus - A_plus
        cooling = net > 0
        (;
            A_plus=A_plus * time_unit^-1,
            A_minus=A_minus * time_unit^-1,
            nbar=cooling ? (A_plus + h) / net : NaN,
            τ_c=(cooling ? 1 / net : NaN) * time_unit,
        )
    end
end

OpticalBloch.cooling_rates(model::LindbladModel, mc::MotionalCoupling; kwargs...) =
    only(OpticalBloch.cooling_rates(model, [mc]; kwargs...))

# --- Levels.OpticalBloch: the full internal ⊗ motional model ----------------------

function OpticalBloch.motional_model(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    num_fock,
    lamb_dicke_order::Real=1,
    heating_rate=0.0u"s^-1",
    nodes=nothing,
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "motional_model")
    validate_lamb_dicke_order(lamb_dicke_order)
    order = Float64(lamb_dicke_order)
    num_modes = length(mcs)
    num_modes >= 1 || throw(ArgumentError("At least one mode is required"))
    heating_rates = mode_heating_rates(heating_rate, num_modes)
    nf = num_fock isa Integer ? fill(Int(num_fock), num_modes) : collect(Int, num_fock)
    length(nf) == num_modes ||
        throw(ArgumentError("num_fock must be one truncation or one per mode"))
    all(>(1), nf) ||
        throw(ArgumentError("Each Fock truncation must hold at least 2 states"))
    emission = mcs[1].emission
    if order > 1
        all(mc.emission == emission for mc in mcs) || throw(
            ArgumentError(
                "All modes must share one emission setting (:exact or :isotropic) " *
                "beyond first order in the Lamb–Dicke parameters",
            ),
        )
        emission == :constant && throw(
            ArgumentError(
                "Beyond first order in the Lamb–Dicke parameters the recoil kicks " *
                "need the emission-direction distribution, not only its second " *
                "moment; build the MotionalCoupling with recoil_moment = :exact " *
                "or :isotropic",
            ),
        )
    end
    num_nodes = isnothing(nodes) ? (isinf(order) ? 12 : Int(order) + 1) : Int(nodes)
    num_states = length(model.basis)

    # Sparse throughout: the tensor products (and hence the Liouvillian, of
    # dimension (N Π num_fock)²) would otherwise be dense.
    internal(M) = to_sparse(QuantumObject(M, model.basis; time_unit))
    motional(D) = QuantumObject(D, Operator(), Tuple(nf))
    eye_int = qeye(num_states)
    eyes = [qeye(n) for n in nf]
    # Operator on the full space with `op` in slot `slot` (0 = internal states).
    function embed(op, slot)
        factors = Any[eye_int; eyes...]
        factors[slot+1] = op
        tensor(factors...)
    end
    rate(x) = ustrip(time_unit^-1, x)

    # Hamiltonian: the frame diagonal, the mode energies, and every beam's
    # coupling dressed with the (truncated) joint displacement operator of its
    # projected Lamb–Dicke factors.
    H = embed(internal(Matrix(Diagonal(diag(model.hamiltonian)))), 0)
    for (m, mc) in enumerate(mcs)
        H += rate(mc.mode.frequency) * embed(num(nf[m]), m)
    end
    for (b, C) in enumerate(model.couplings)
        ηs = [mc.projected_lamb_dicke[b] for mc in mcs]
        T = tensor(internal(C), motional(joint_displacement(ηs, nf; order))) / 2
        H += T + T'
    end

    # Jump operators: dephasing on the internal states only; spontaneous
    # emission with its recoil — the separate first-order kick operators, or the
    # emission-direction quadrature set of displaced decay operators.
    c_ops = QuantumObject[]
    labels = Any[]
    directions = [mc.mode.direction for mc in mcs]
    for (k, (L, label)) in enumerate(zip(model.jump_operators, model.jump_labels))
        if label isa DephasingLabel
            push!(c_ops, embed(internal(L), 0))
            push!(labels, label)
            continue
        end
        if order == 1
            push!(c_ops, embed(internal(L), 0))
            push!(labels, label)
            for (m, mc) in enumerate(mcs)
                R = mc.recoil_operators[findfirst(==(k), mc.recoil_indices)]
                a = destroy(nf[m])
                factors = Any[internal(R); eyes...]
                factors[m+1] = a + a'
                push!(c_ops, tensor(factors...))
                push!(labels, RecoilLabel(label, m, 0))
            end
        else
            rank =
                emission == :isotropic ? 0 :
                Levels.multipole_rank(label.lower, label.upper)
            q = emission == :isotropic ? 0 : label.q
            hs, ws = emission_rule(
                rank,
                q,
                num_modes == 1 ? directions[1] : directions;
                nodes=num_nodes,
            )
            η0s = [
                mc.recoil_lamb_dicke[findfirst(==(k), mc.recoil_indices)] for mc in mcs
            ]
            Ls = internal(L)
            for j in eachindex(ws)
                h = num_modes == 1 ? [hs[j]] : hs[j, :]
                D = joint_displacement(-η0s .* h, nf; order)
                push!(c_ops, tensor(sqrt(ws[j]) * Ls, motional(D)))
                push!(labels, RecoilLabel(label, 0, j))
            end
        end
    end

    # Heating bath.
    for m in 1:num_modes
        h_bath = sqrt(rate(heating_rates[m]))
        h_bath > 0 || continue
        a = destroy(nf[m])
        push!(c_ops, h_bath * embed(a, m))
        push!(labels, HeatingLabel(m, false))
        push!(c_ops, h_bath * embed(a', m))
        push!(labels, HeatingLabel(m, true))
    end
    MotionalModel(
        H,
        c_ops,
        labels,
        num_states,
        nf,
        [mc.mode for mc in mcs],
        order,
        time_unit,
    )
end

OpticalBloch.motional_model(model::LindbladModel, mc::MotionalCoupling; kwargs...) =
    OpticalBloch.motional_model(model, [mc]; kwargs...)

# The phonon-number operator of one mode on the full space.
function number_operator(mm::MotionalModel, mode::Integer)
    factors = Any[qeye(mm.num_states); [qeye(n) for n in mm.num_fock]...]
    factors[mode+1] = num(mm.num_fock[mode])
    tensor(factors...)
end

function reduced_populations(ρ::QuantumObject, mm::MotionalModel, mode::Integer)
    p = real.(diag(ptrace(ρ, mode + 1).data))
    p ./ sum(p)
end

function OpticalBloch.fock_populations(ρ::QuantumObject, mm::MotionalModel)
    [reduced_populations(ρ, mm, m) for m in eachindex(mm.num_fock)]
end

function OpticalBloch.mean_phonon_number(
    ρ::QuantumObject,
    mm::MotionalModel;
    method::Symbol=:direct,
)
    method in (:direct, :thermal_ratio) ||
        throw(ArgumentError("method must be :direct or :thermal_ratio, got $method"))
    nbars = map(eachindex(mm.num_fock)) do m
        p = reduced_populations(ρ, mm, m)
        if method == :direct
            sum((n - 1) * p[n] for n in eachindex(p))
        else
            # Ratio estimator of r = n̄/(n̄ + 1), exact for a (truncated) thermal
            # distribution, as p_{n+1}/p_n = r for every n. The rate-equation
            # steady state stays thermal under truncation, as that of any
            # truncated reversible Markov chain is the renormalised original one
            # (Kelly, Reversibility and Stochastic Networks (1979), Sec. 1.6).
            r = sum(p[2:end]) / sum(p[1:(end-1)])
            r < 1 ? r / (1 - r) : Inf
        end
    end
    length(nbars) == 1 ? only(nbars) : nbars
end

function OpticalBloch.thermal_state(mm::MotionalModel, internal::QuantumObject, nbars)
    ρ_int = isket(internal) ? ket2dm(internal) : internal
    nb =
        nbars isa Real ? fill(Float64(nbars), length(mm.num_fock)) :
        collect(Float64, nbars)
    length(nb) == length(mm.num_fock) ||
        throw(ArgumentError("One mean occupation per mode is required"))
    thermals = [
        QuantumObject(Diagonal(thermal_populations(nbar, n)), Operator(), n) for
        (n, nbar) in zip(mm.num_fock, nb)
    ]
    tensor(ρ_int, thermals...)
end

function OpticalBloch.cooling_curve(
    mm::MotionalModel,
    ρ0::QuantumObject,
    ts;
    fock::Bool=false,
    kwargs...,
)
    times = eltype(ts) <: Quantity ? ustrip.(mm.time_unit, ts) : collect(float.(ts))
    num_modes = length(mm.num_fock)
    N_ops = [number_operator(mm, m) for m in 1:num_modes]
    if fock
        sol = mesolve(mm.H, ρ0, times, mm.c_ops; progress_bar=Val(false), kwargs...)
        nbar = [real(tr(N_ops[m].data * ρ.data)) for ρ in sol.states, m in 1:num_modes]
        populations = [
            reduce(hcat, [reduced_populations(ρ, mm, m) for ρ in sol.states]) for
            m in 1:num_modes
        ]
        return (; t=times, nbar, fock=populations)
    end
    sol = mesolve(
        mm.H,
        ρ0,
        times,
        mm.c_ops;
        e_ops=N_ops,
        progress_bar=Val(false),
        kwargs...,
    )
    (; t=times, nbar=Matrix(transpose(real.(sol.expect))), fock=nothing)
end

OpticalBloch.BandLiouvillian(mm::MotionalModel; bandwidth=nothing) = BandLiouvillian(
    mm.H.data,
    [c.data for c in mm.c_ops],
    mm.num_states,
    mm.num_fock;
    bandwidth,
)

function OpticalBloch.BandLiouvillian(
    model::LindbladModel;
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "BandLiouvillian")
    BandLiouvillian(
        hamiltonian_qobj(model; time_unit).data,
        [L.data for L in jump_qobjs(model; time_unit)],
        length(model.basis),
        Int[],
    )
end

function OpticalBloch.IntegratedTransientSolver(
    mm::MotionalModel;
    bandwidth=nothing,
    observables=[
        phonon_number_operator(mm.num_states, mm.num_fock, m) for
        m in eachindex(mm.num_fock)
    ],
    atol::Real=1e-9,
)
    IntegratedTransientSolver(BandLiouvillian(mm; bandwidth), observables; atol)
end

OpticalBloch.IntegratedTransientSolver(
    model::LindbladModel,
    observables;
    time_unit::Unitful.Units=Unitful.µs,
    atol::Real=1e-9,
) = IntegratedTransientSolver(BandLiouvillian(model; time_unit), observables; atol)

function OpticalBloch.cooling_time(
    H::QuantumObject,
    c_ops;
    time_unit::Unitful.Units=Unitful.µs,
    mode_sizes=nothing,
    mode=nothing,
    eigvals::Int=12,
    dense_limit::Int=2500,
    weight_threshold::Real=0.1,
)
    L = liouvillian(H, c_ops)
    dim = size(L.data, 1)
    values, vectors = if dim <= dense_limit
        e = eigen(Matrix(L.data))
        e.values, e.vectors
    else
        r = eigsolve(L; sigma=0.0 + 0.0im, eigvals)
        r.values, r.vectors
    end
    # Phonon-number operator of the selected mode (or all modes) over the
    # motional factors (slot 1 is internal).
    sizes = if isnothing(mode_sizes)
        d = H.dims
        collect(Int, d isa Tuple ? first(d) : d)[2:end]
    else
        collect(Int, mode_sizes)
    end
    num_states = isqrt(dim) ÷ prod(sizes)
    N_op = phonon_number_operator(num_states, sizes, mode)
    weights = [
        abs(tr(N_op * vec2mat(vectors[:, i]))) / norm(vectors[:, i]) for
        i in axes(vectors, 2)
    ]
    order = sortperm(abs.(real.(values)))
    threshold = weight_threshold * maximum(weights)
    # Skip the stationary state (zero rate); the first remaining mode that moves
    # the phonon number appreciably is the slowest relaxation of ⟨n⟩ (coherence
    # modes acquire a weight of order η² through the sideband coupling and are
    # excluded by the relative threshold).
    for i in order[2:end]
        weights[i] >= threshold && return (1 / abs(real(values[i]))) * time_unit
    end
    throw(
        ErrorException(
            "No relaxation mode with phonon-number weight above the threshold " *
            "found; increase eigvals or lower weight_threshold",
        ),
    )
end

OpticalBloch.cooling_time(mm::MotionalModel; kwargs...) = OpticalBloch.cooling_time(
    mm.H,
    mm.c_ops;
    time_unit=mm.time_unit,
    mode_sizes=mm.num_fock,
    kwargs...,
)

# --- Levels.OpticalBloch: the common cooling-metrics front-end -------------------

function OpticalBloch.cooling_metrics(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    method::AdiabaticElimination;
    heating_rate=0.0u"s^-1",
    ρ=steadystate(model),
    time_unit::Unitful.Units=Unitful.µs,
)
    rates = OpticalBloch.cooling_rates(model, mcs; ρ, heating_rate, time_unit)
    [
        CoolingMetrics(method, r.nbar, r.τ_c, r.A_plus, r.A_minus, nothing, nothing) for
        r in rates
    ]
end

function OpticalBloch.cooling_metrics(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    method::LiouvillianSpectrum;
    heating_rate=0.0u"s^-1",
    ρ=nothing,
    time_unit::Unitful.Units=Unitful.µs,
)
    mm = OpticalBloch.motional_model(
        model,
        mcs;
        num_fock=method.num_fock,
        lamb_dicke_order=method.lamb_dicke_order,
        heating_rate,
        time_unit,
    )
    ρ_ss = steadystate(mm)
    nbars = OpticalBloch.fock_populations(ρ_ss, mm)
    map(eachindex(mcs)) do m
        nbar = sum((n - 1) * p for (n, p) in enumerate(nbars[m]))
        τ = OpticalBloch.cooling_time(
            mm;
            mode=m,
            eigvals=method.eigvals,
            dense_limit=method.dense_limit,
            weight_threshold=method.weight_threshold,
        )
        CoolingMetrics(method, nbar, τ, nothing, nothing, mm, nothing)
    end
end

function OpticalBloch.cooling_metrics(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    method::IntegratedTransient;
    heating_rate=0.0u"s^-1",
    ρ=steadystate(model),
    time_unit::Unitful.Units=Unitful.µs,
)
    num_modes = length(mcs)
    nbar_ini =
        method.nbar_ini isa Real ? fill(method.nbar_ini, num_modes) : method.nbar_ini
    length(nbar_ini) == num_modes ||
        throw(ArgumentError("nbar_ini must be one occupation or one per mode"))
    num_fock = isnothing(method.num_fock) ? fock_truncation.(nbar_ini) : method.num_fock
    mm = OpticalBloch.motional_model(
        model,
        mcs;
        num_fock,
        lamb_dicke_order=method.lamb_dicke_order,
        heating_rate,
        nodes=method.nodes,
        time_unit,
    )
    solver = IntegratedTransientSolver(mm; bandwidth=method.bandwidth)
    ρ0 = OpticalBloch.thermal_state(Matrix(ρ.data), mm.num_fock, nbar_ini)
    map(1:num_modes) do m
        τ = integral_relaxation_time(solver, ρ0, m) * time_unit
        CoolingMetrics(method, steady_value(solver, m), τ, nothing, nothing, mm, solver)
    end
end

OpticalBloch.cooling_metrics(
    model::LindbladModel,
    mc::MotionalCoupling,
    method;
    kwargs...,
) = only(OpticalBloch.cooling_metrics(model, [mc], method; kwargs...))

end # module
