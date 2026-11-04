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

using LinearAlgebra: Diagonal, I, diag, tr
using Unitful
using Levels: Levels, Eigenstructure, StateBasis, stateindex
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
    expectation_row,
    fock_truncation,
    integral_relaxation_time,
    joint_displacement,
    level_projector,
    mode_heating_rates,
    phonon_number_operator,
    relaxation_modes,
    require_static,
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
    QuantumObject(A::AbstractMatrix, s::Eigenstructure; time_unit = u"µs")
    QuantumObject(v::AbstractVector, s::Eigenstructure; time_unit = u"µs")

The same conversions over the basis of an [`Levels.Eigenstructure`](@ref),
for matrices and vectors expressed in its field eigenbasis.
"""
QuantumToolbox.QuantumObject(
    A::AbstractVecOrMat,
    s::Eigenstructure;
    time_unit::Unitful.Units=Unitful.µs,
) = QuantumObject(A, s.basis; time_unit)

"""
    basis(b::StateBasis, state)
    basis(b::StateBasis, level, m)
    basis(s::Eigenstructure, state)
    basis(s::Eigenstructure, level, m)

Returns the ket `QuantumObject` for the given state of the [`StateBasis`](@ref),
or of the basis of an [`Levels.Eigenstructure`](@ref) (cf.
[`Levels.stateindex`](@ref) for the accepted state forms).
"""
QuantumToolbox.fock(b::StateBasis, state) =
    QuantumToolbox.fock(length(b), stateindex(b, state) - 1)
QuantumToolbox.fock(b::StateBasis, level, m) =
    QuantumToolbox.fock(length(b), stateindex(b, level, m) - 1)
QuantumToolbox.fock(s::Eigenstructure, state) = QuantumToolbox.fock(s.basis, state)
QuantumToolbox.fock(s::Eigenstructure, level, m) =
    QuantumToolbox.fock(s.basis, level, m)

"""
    projection(b::StateBasis, ket[, bra])
    projection(s::Eigenstructure, ket[, bra])

Returns the projection operator ``|ket⟩⟨bra|`` between states of the
[`StateBasis`](@ref), or of the basis of an [`Levels.Eigenstructure`](@ref)
(cf. [`Levels.stateindex`](@ref) for the accepted state forms); `bra` defaults
to `ket`, giving the population projector.
"""
QuantumToolbox.projection(b::StateBasis, ket, bra) =
    QuantumToolbox.projection(length(b), stateindex(b, ket) - 1, stateindex(b, bra) - 1)
QuantumToolbox.projection(b::StateBasis, ket) = QuantumToolbox.projection(b, ket, ket)
QuantumToolbox.projection(s::Eigenstructure, ket, bra) =
    QuantumToolbox.projection(s.basis, ket, bra)
QuantumToolbox.projection(s::Eigenstructure, ket) =
    QuantumToolbox.projection(s.basis, ket)

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
    steadystate(model::LindbladModel; time_unit = u"µs", harmonic_cutoff = 4, solver, kwargs...)
    steadystate(mm::MotionalModel; solver, kwargs...)

Returns the steady-state internal density matrix of the
[`Levels.OpticalBloch.LindbladModel`](@ref), or the steady state of the full
internal ⊗ motional master equation of a
[`Levels.OpticalBloch.MotionalModel`](@ref). A static model is solved directly
(`QuantumToolbox.steadystate`); an internal model with beat-note harmonics at
one frequency is solved for its periodic steady state with
`steadystate_fourier` (`harmonic_cutoff` harmonics of ρ, cf. its `n_max`), of
which the period-averaged k = 0 component is returned. `solver` defaults to
the direct sparse factorisation `SteadyStateLinearSolver(; alg = nothing)`;
further keywords pass through.
"""
function QuantumToolbox.steadystate(
    model::LindbladModel;
    time_unit::Unitful.Units=Unitful.µs,
    harmonic_cutoff::Int=4,
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
        n_max=harmonic_cutoff,
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

# --- Levels.OpticalBloch: the full internal ⊗ motional model ----------------------

function OpticalBloch.motional_model(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    num_fock,
    lamb_dicke_order::Real=1,
    heating_rate=0.0u"s^-1",
    quadrature_order=nothing,
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
    num_nodes =
        isnothing(quadrature_order) ? (isinf(order) ? 12 : Int(order) + 1) :
        Int(quadrature_order)
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
            rank = emission == :isotropic ? 0 : Levels.multipole_rank(label.multipole)
            q = emission == :isotropic ? 0 : label.q
            hs, ws = emission_rule(
                rank,
                q,
                num_modes == 1 ? directions[1] : directions;
                quadrature_order=num_nodes,
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
        return (; t=times .* mm.time_unit, nbar, fock=populations)
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
    (; t=times .* mm.time_unit, nbar=Matrix(transpose(real.(sol.expect))), fock=nothing)
end

OpticalBloch.BandLiouvillian(mm::MotionalModel; bandwidth=nothing) = BandLiouvillian(
    mm.H.data,
    [c.data for c in mm.c_ops],
    mm.num_states,
    mm.num_fock;
    bandwidth,
    time_unit=mm.time_unit,
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
        Int[];
        time_unit,
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

# --- Levels.OpticalBloch: the common cooling-metrics front-end -------------------

function OpticalBloch.cooling_metrics(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling},
    method::AdiabaticElimination;
    heating_rate=0.0u"s^-1",
    ρ=steadystate(model),
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "AdiabaticElimination")
    heating_rates = mode_heating_rates(heating_rate, length(mcs))
    L = liouvillian(model; time_unit).data
    ρm = Matrix(ρ.data)
    map(zip(mcs, heating_rates)) do (mc, rate)
        h = ustrip(time_unit^-1, rate)
        A_plus, A_minus = mode_rates(L, ρm, mc, time_unit)
        net = A_minus - A_plus
        cooling = net > 0
        CoolingMetrics(
            method,
            cooling ? (A_plus + h) / net : NaN,
            (cooling ? 1 / net : NaN) * time_unit,
            A_plus * time_unit^-1,
            A_minus * time_unit^-1,
            nothing,
            nothing,
        )
    end
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
    # One factorisation serves the stationary state and the slow spectrum.
    bl = BandLiouvillian(mm; bandwidth=method.bandwidth)
    modes = relaxation_modes(
        bl;
        num_eigenpairs=method.num_eigenpairs,
        dense_limit=method.dense_limit,
    )
    nbars = OpticalBloch.fock_populations(bl, modes.steady)
    map(eachindex(mcs)) do m
        nbar = sum((n - 1) * p for (n, p) in enumerate(nbars[m]))
        # The slowest computed mode moving this mode's phonon number (the vectors
        # have unit norm); the threshold passes over the coherence modes.
        a = expectation_row(bl, phonon_number_operator(bl.num_states, bl.num_fock, m))
        weights = [abs(transpose(a) * v) for v in eachcol(modes.vectors)]
        threshold = method.weight_threshold * maximum(weights)
        slowest = findfirst(>=(threshold), weights)
        isnothing(slowest) && error(
            "No relaxation mode with phonon-number weight above the threshold found; " *
            "increase num_eigenpairs or lower weight_threshold",
        )
        τ = -1 / real(ustrip(time_unit^-1, modes.values[slowest])) * time_unit
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
    initial_nbar =
        method.initial_nbar isa Real ? fill(method.initial_nbar, num_modes) :
        method.initial_nbar
    length(initial_nbar) == num_modes ||
        throw(ArgumentError("initial_nbar must be one occupation or one per mode"))
    num_fock =
        isnothing(method.num_fock) ? fock_truncation.(initial_nbar) : method.num_fock
    mm = OpticalBloch.motional_model(
        model,
        mcs;
        num_fock,
        lamb_dicke_order=method.lamb_dicke_order,
        heating_rate,
        quadrature_order=method.quadrature_order,
        time_unit,
    )
    solver = IntegratedTransientSolver(mm; bandwidth=method.bandwidth)
    ρ0 = OpticalBloch.thermal_state(Matrix(ρ.data), mm.num_fock, initial_nbar)
    map(1:num_modes) do m
        τ = integral_relaxation_time(solver, ρ0, m)
        CoolingMetrics(method, steady_value(solver, m), τ, nothing, nothing, mm, solver)
    end
end

end # module
