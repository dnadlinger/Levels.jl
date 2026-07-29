# Exact reference engine: one-period propagator (monodromy matrix) of the full
# time-dependent Schrödinger equation, and sideband Rabi frequencies from Floquet
# quasienergy gaps — nonperturbative in every coupling.

using SciMLBase: ODEProblem, solve
using OrdinaryDiffEqVerner: Vern8
import Optim

"""
Returns a fast unit-stripped in-place evaluator `ham!(H, t)` of the full
rotating-frame Hamiltonian of the model — the single point where the exact
engine leaves unitful quantities (angular ``s^{-1}``, time in seconds), directly
ahead of the ODE solves.
"""
function stripped_hamiltonian(dt::DrivenTransition, δ, modulation)
    Ω = ustrip(u"s^-1", dt.drive_frequency)
    L = ustrip.(u"s^-1", dt.coupling)
    amplitudes = [ustrip.(u"s^-1", drive.amplitude) for drive in dt.drives]
    phases = [drive.phase for drive in dt.drives]
    diagonal = ustrip.(u"s^-1", dt.frame)
    diagonal[dt.upper_range] .-= ustrip(u"s^-1", δ)
    phase_fn = phase_function(modulation, dt.drive_frequency)

    function ham!(H, t)
        f = 0.5 * cis(-phase_fn(t))
        @inbounds for k in axes(H, 2), i in axes(H, 1)
            H[i, k] = f * L[i, k] + conj(f * L[k, i])
        end
        for (amplitude, phase) in zip(amplitudes, phases)
            c = cos(Ω * t + phase)
            @inbounds H .+= c .* amplitude
        end
        @inbounds for i in eachindex(diagonal)
            H[i, i] += diagonal[i]
        end
        H
    end
end

"""
    monodromy(dt::DrivenTransition; δ = 0, modulation = none)

Computes the one-drive-period propagator (monodromy matrix) of
``i \\dot{ψ} = H(t) ψ`` at laser detuning `δ` from the probed transition, with
the given laser phase modulation.

Returns the (dimensionless, complex) propagator matrix over the model's basis.
"""
function monodromy(
    dt::DrivenTransition;
    δ=zero(dt.drive_frequency),
    modulation=HarmonicPhaseModulation(),
)
    period = ustrip(u"s", 2π / dt.drive_frequency)
    n = length(dt.basis)
    ham! = stripped_hamiltonian(dt, δ, modulation)
    H = zeros(ComplexF64, n, n)

    function tdse!(dU, U, _, t)
        ham!(H, t)
        mul!(dU, H, U, -1.0im, 0.0im)
    end

    problem = ODEProblem(tdse!, Matrix{ComplexF64}(I, n, n), (0.0, period))

    # High-accuracy reference propagation settings.
    solution = solve(
        problem,
        Vern8();
        reltol=1e-10,
        abstol=1e-12,
        save_everystep=false,
        save_start=false,
        dense=false,
    )

    last(solution.u)
end

"""
Quasienergies (wrapped into ``(-Ω/2, Ω/2]``, plain angular ``s^{-1}``) and
Floquet eigenvectors at ``t = 0`` from the monodromy matrix ``U(T) = e^{-iεT}``.
"""
function floquet_from_monodromy(U, Ω)
    λ, V = eigen(U)
    ε = angle.(λ) .* (-Ω / 2π)
    ε, V
end

"""
Wrapped quasienergy gap (plain angular ``s^{-1}``) between the two Floquet states
carrying the probed lower- and upper-state characters. At a sideband resonance,
the minimum of this gap over the detuning equals the effective Rabi frequency.
"""
function probe_gap(dt::DrivenTransition, δ, modulation)
    Ω = ustrip(u"s^-1", dt.drive_frequency)
    ε, V = floquet_from_monodromy(monodromy(dt; δ, modulation), Ω)
    w_lower = abs2.(V[dt.lower, :])
    w_upper = abs2.(V[dt.upper, :])
    i_lower, i_upper = argmax(w_lower), argmax(w_upper)
    if i_lower == i_upper   # strongly hybridised: take the two largest combined weights
        i_lower, i_upper = partialsortperm(w_lower .+ w_upper, 1:2; rev=true)
    end
    Δ = ε[i_lower] - ε[i_upper]
    abs(mod(Δ + Ω / 2, Ω) - Ω / 2)
end

"""
    exact_sideband(dt::DrivenTransition; sideband, modulation = none,
                   bracket_scale = 8.0, nharm = 8, ngrid = 256)

Computes the exact sideband Rabi frequency by minimising the Floquet quasienergy
gap over the laser detuning, bracketed around the dressed-state prediction of
[`sideband_amplitude`](@ref) (whose `nharm`/`ngrid` settings are passed through).

Returns `(Ω, δ_res)`. This is the gold-standard cross-check for the dressed-state
engine, nonperturbative in every coupling, so it also includes all laser-induced
shifts; it is several orders of magnitude slower.
"""
function exact_sideband(
    dt::DrivenTransition;
    sideband::Int,
    modulation=HarmonicPhaseModulation(),
    bracket_scale::Real=8.0,
    nharm::Int=8,
    ngrid::Int=256,
)
    estimate = sideband_amplitude(dress(dt; nharm, ngrid); sideband, modulation)

    δ_estimate = ustrip(u"s^-1", estimate.δ_res)
    halfwidth = max(bracket_scale * ustrip(u"s^-1", estimate.Ω), 2π * 100)
    result = Optim.optimize(
        δ -> probe_gap(dt, δ * u"s^-1", modulation),
        δ_estimate - halfwidth,
        δ_estimate + halfwidth,
        Optim.Brent();
        abs_tol=halfwidth * 1e-7,
    )

    (
        Ω=uconvert(u"ms^-1", Optim.minimum(result) * u"s^-1"),
        δ_res=uconvert(u"µs^-1", Optim.minimizer(result) * u"s^-1"),
    )
end

"""
    stroboscopic_populations(dt::DrivenTransition, num_periods::Int;
                             δ = 0, modulation = none, initial = dt.lower)

Propagates the model exactly from the basis state with index `initial`
(defaulting to the probed lower state) and returns the populations of all basis
states at the stroboscopic times ``k T`` for ``k = 0 … `` `num_periods`, as a
`length(basis) × (num_periods + 1)` matrix — a direct TDSE Rabi flop for
validation. Rows follow the basis order (cf. [`Levels.stateindex`](@ref)).
"""
function stroboscopic_populations(
    dt::DrivenTransition,
    num_periods::Int;
    δ=zero(dt.drive_frequency),
    modulation=HarmonicPhaseModulation(),
    initial::Int=dt.lower,
)
    U = monodromy(dt; δ, modulation)
    n = length(dt.basis)
    ψ = zeros(ComplexF64, n)
    ψ[initial] = 1
    populations = Matrix{Float64}(undef, n, num_periods + 1)
    populations[:, 1] .= abs2.(ψ)
    for k in 1:num_periods
        ψ = U * ψ
        populations[:, k+1] .= abs2.(ψ)
    end
    populations
end

export exact_sideband, stroboscopic_populations
public monodromy
