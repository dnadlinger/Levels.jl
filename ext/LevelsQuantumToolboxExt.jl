# Conversions from Levels.jl objects to QuantumToolbox.jl `QuantumObject`s, so
# that the QuantumToolbox solvers (`mesolve`, `steadystate_fourier`, …) can be
# applied to models built from Levels data.
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
using Levels: StateBasis, stateindex
using Levels.PeriodicDynamics: PeriodicDynamics, DrivenTransition
using Levels.OpticalBloch:
    OpticalBloch, DecayLabel, LindbladModel, MotionalCoupling, level_projector
using QuantumToolbox:
    QuantumToolbox,
    Ket,
    Operator,
    QuantumObject,
    SteadyStateLinearSolver,
    destroy,
    eigsolve,
    liouvillian,
    mat2vec,
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

Returns the Liouvillian superoperator of the internal-state master equation of
the [`Levels.OpticalBloch.LindbladModel`](@ref) (static models only), in
inverse `time_unit`.
"""
function QuantumToolbox.liouvillian(
    model::LindbladModel;
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "liouvillian")
    liouvillian(hamiltonian_qobj(model; time_unit), jump_qobjs(model; time_unit))
end

"""
    steadystate(model::LindbladModel; time_unit = u"µs", num_harmonics = 4, solver, kwargs...)

Returns the steady-state internal density matrix of the
[`Levels.OpticalBloch.LindbladModel`](@ref). A static model is solved directly
(`QuantumToolbox.steadystate`); a model with beat-note harmonics at one frequency
is solved for its periodic steady state with `steadystate_fourier`
(`num_harmonics` harmonics of ρ, cf. its `n_max`), of which the period-averaged
k = 0
component is returned. `solver` defaults to the direct sparse factorisation
`SteadyStateLinearSolver(; alg = nothing)`; further keywords pass through.
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

# --- Levels.OpticalBloch: motional layer ------------------------------------------

# Adiabatic-elimination rates of one mode from the stripped internal Liouvillian
# `L` and steady state `ρ` (plain matrices, inverse time_unit).
function mode_rates(L, ρ, mc::MotionalCoupling, model::LindbladModel, time_unit)
    H_sb = ustrip.(time_unit^-1, mc.sideband_hamiltonian)
    ω_m = ustrip(time_unit^-1, mc.mode.frequency)
    y = mat2vec(H_sb * ρ)
    spectrum(ω) = tr(H_sb * vec2mat((-L + (im * ω) * I) \ y))
    diffusion = 0.0
    for (k, R) in zip(mc.recoil_indices, mc.recoil_operators)
        Rs = ustrip.(time_unit^(-1 // 2), R)
        diffusion += real(tr(Rs' * Rs * ρ))
    end
    A_plus = 2 * real(spectrum(ω_m)) + diffusion
    A_minus = 2 * real(spectrum(-ω_m)) + diffusion
    A_plus, A_minus
end

function OpticalBloch.cooling_rates(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    ρ=steadystate(model),
    heating_rate=0.0u"s^-1",
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "cooling_rates")
    L = liouvillian(model; time_unit).data
    ρm = Matrix(ρ.data)
    h = ustrip(time_unit^-1, heating_rate)
    map(mcs) do mc
        A_plus, A_minus = mode_rates(L, ρm, mc, model, time_unit)
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

function OpticalBloch.motional_model(
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    num_fock,
    time_unit::Unitful.Units=Unitful.µs,
)
    require_static(model, "motional_model")
    num_modes = length(mcs)
    sizes =
        num_fock isa Integer ? fill(Int(num_fock), num_modes) : collect(Int, num_fock)
    length(sizes) == num_modes ||
        throw(ArgumentError("num_fock must be one truncation or one per mode"))
    all(>(1), sizes) || throw(ArgumentError("Each Fock truncation must be at least 2"))

    # Sparse throughout: the tensor products (and hence the Liouvillian, of
    # dimension (N Π num_fock)²) would otherwise be dense.
    internal(M) = to_sparse(QuantumObject(M, model.basis; time_unit))
    eye_int = qeye(length(model.basis))
    eyes = [qeye(size) for size in sizes]
    # Operator on the full space with `op` in slot `slot` (0 = internal states).
    function embed(op, slot)
        factors = Any[eye_int; eyes...]
        factors[slot+1] = op
        tensor(factors...)
    end
    # `op_int ⊗ (a_m + a_m†)` with identities on the other modes.
    function sideband(op_int, m)
        a = destroy(sizes[m])
        factors = Any[op_int; eyes...]
        factors[m+1] = a + a'
        tensor(factors...)
    end

    H = embed(internal(model.hamiltonian), 0)
    for (m, mc) in enumerate(mcs)
        a = destroy(sizes[m])
        H += ustrip(time_unit^-1, mc.mode.frequency) * embed(a' * a, m)
        H += sideband(internal(mc.sideband_hamiltonian), m)
    end
    c_ops = QuantumObject[embed(internal(L), 0) for L in model.jump_operators]
    for (m, mc) in enumerate(mcs)
        for R in mc.recoil_operators
            push!(c_ops, sideband(internal(R), m))
        end
    end
    (; H, c_ops)
end

OpticalBloch.motional_model(model::LindbladModel, mc::MotionalCoupling; kwargs...) =
    OpticalBloch.motional_model(model, [mc]; kwargs...)

function OpticalBloch.mean_phonon_number(
    ρ::QuantumObject,
    model::LindbladModel,
    mcs::AbstractVector{<:MotionalCoupling};
    method::Symbol=:direct,
)
    method in (:direct, :thermal_ratio) ||
        throw(ArgumentError("method must be :direct or :thermal_ratio, got $method"))
    map(1:length(mcs)) do m
        p = real.(diag(ptrace(ρ, m + 1).data))
        p ./= sum(p)
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
end

OpticalBloch.mean_phonon_number(
    ρ::QuantumObject,
    model::LindbladModel,
    mc::MotionalCoupling;
    kwargs...,
) = only(OpticalBloch.mean_phonon_number(ρ, model, [mc]; kwargs...))

function OpticalBloch.cooling_time(
    H::QuantumObject,
    c_ops;
    time_unit::Unitful.Units=Unitful.µs,
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
    # Total phonon-number operator over the motional factors (slot 1 is internal).
    d = H.dims
    sizes = collect(Int, d isa Tuple ? first(d) : d)
    N_total = sum(
        tensor((k == m ? num(sizes[k]) : qeye(sizes[k]) for k in eachindex(sizes))...) for m in 2:length(sizes)
    )
    weights = [
        abs(tr(N_total.data * vec2mat(vectors[:, i]))) / norm(vectors[:, i]) for
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

end # module
