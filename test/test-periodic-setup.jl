# Shared scenario builder for the PeriodicDynamics tests: the ⁸⁸Sr⁺ 674 nm setup
# of TrapRFSim.jl (thesis §5.6.3) — B₀ = 0.5 mT, Ω_rf = 2π·50 MHz, probe beam
# k ⊥ B₀ with linear polarisation at 45°, so Δm = ±1, ±2 couple but Δm = 0 does
# not, matching the eight transitions of thesis fig. 5.12.
@testsnippet PeriodicSetup begin
    using LinearAlgebra
    using Unitful
    using Levels.PeriodicDynamics

    const BASIS = StateBasis(["S_1/2", "D_5/2"])
    const B0 = 0.5u"mT"
    const Ω_RF = 2π * 50.0u"MHz"
    const Ω0 = 2π * 100.0u"kHz"
    const μ_B = Unitful.q * u"ħ" / (2 * Unitful.me)

    # The eight probed transitions, sorted by increasing magnetic-field
    # sensitivity χ as in thesis fig. 5.12.
    const TRANSITIONS = sort!(
        state_pairs("S_1/2", "D_5/2"; Δm=[-2, -1, 1, 2]);
        by=t -> zeeman_sensitivity(sr88, t),
    )

    function probe_coupling(probe; Ω0=Ω0)
        n, ε = beam_vectors(π / 2, π / 4)
        rabi_normalised(
            quadrupole_couplings(BASIS, "S_1/2", "D_5/2", ε, n),
            BASIS,
            probe,
            Ω0,
        )
    end

    function make_transition(
        probe;
        drives=HarmonicDrive[],
        drive_frequency=Ω_RF,
        static_field=B0,
        coupling=probe_coupling(probe),
    )
        DrivenTransition(
            sr88,
            BASIS,
            probe;
            drive_frequency,
            static_field,
            coupling,
            drives,
        )
    end
end
