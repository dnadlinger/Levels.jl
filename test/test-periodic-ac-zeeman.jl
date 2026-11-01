# Shared scenario builder for the ac Zeeman tests: the ⁴³Ca⁺ 729 nm clock
# transition |S_1/2, F=4, m_F=+4⟩ → |D_5/2, F=4, m_F=+3⟩ at its upper
# field-insensitive point B₀ ≈ 4.96 G ([Benhelm2007]; 0.49582 mT with this
# dataset), driven by trap-rf-frequency ac magnetic fields — the scenario of the
# CaACZeemanShift.jl reference implementation whose physics validation is
# ported here.
@testsnippet Ca43ClockSetup begin
    using LinearAlgebra
    using Unitful
    using Levels.PeriodicDynamics

    const B0 = 0.49582u"mT"
    const Ω_RF = angular(30.0u"MHz")
    const S_PROBE = StateSpec("S_1/2 F=4", 4)
    const D_PROBE = StateSpec("D_5/2 F=4", 3)
    # B_⊥ x̂ + B_∥ ẑ, 50 µT each (the reference README scenario).
    const B_DRIVE = [0.05u"mT", 0.0u"mT", 0.05u"mT"]

    # The atom: both clock manifolds at the static field, solved once.
    clock_structure(b=B0) = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=b)
    const S_CLOCK = clock_structure()
end

@testitem "PT ac Zeeman shift vs reference implementation" tags=[:unit, :fast] setup=[
    Ca43ClockSetup,
] begin
    hz(x) = ustrip(u"Hz", cyclic(x))

    # The stretched S(4,+4) state is an exact product state, so a pure π
    # (B_∥) drive shifts it by exactly zero…
    @test iszero(
        ac_zeeman_shift(S_CLOCK, S_PROBE, Ω_RF, [0.0u"mT", 0.0u"mT", 0.05u"mT"]),
    )
    # …and its whole shift comes from the transverse component.
    r = ac_zeeman_shift(S_CLOCK, S_PROBE => D_PROBE, Ω_RF, B_DRIVE)
    @test ac_zeeman_shift(S_CLOCK, S_PROBE, Ω_RF, [0.05u"mT", 0.0u"mT", 0.0u"mT"]) ≈
          r.lower rtol = 1e-12
    # The state may also be given as a basis index.
    @test ac_zeeman_shift(
        S_CLOCK,
        stateindex(Levels.basis(S_CLOCK), S_PROBE),
        Ω_RF,
        B_DRIVE,
    ) == r.lower

    # No π/σ cross terms at second order: the ∥ and ⊥ shifts add exactly, and
    # the relative temporal phase of the two components is irrelevant.
    δ_par = ac_zeeman_shift(S_CLOCK, D_PROBE, Ω_RF, [0.0u"mT", 0.0u"mT", 0.05u"mT"])
    δ_perp = ac_zeeman_shift(S_CLOCK, D_PROBE, Ω_RF, [0.05u"mT", 0.0u"mT", 0.0u"mT"])
    @test δ_par + δ_perp ≈ r.upper rtol = 1e-12
    drives_dephased = [
        zeeman_drive(S_CLOCK, [0.05u"mT", 0.0u"mT", 0.0u"mT"]; phase=0.0),
        zeeman_drive(S_CLOCK, [0.0u"mT", 0.0u"mT", 0.05u"mT"]; phase=1.234),
    ]
    @test quasienergy_shift(S_CLOCK, D_PROBE, Ω_RF, drives_dephased) ≈ r.upper rtol =
        1e-12
    # A structure over the one manifold gives the same state shift.
    @test ac_zeeman_shift(Eigenstructure(ca43, "D_5/2"; B=B0), D_PROBE, Ω_RF, B_DRIVE) ≈
          r.upper rtol = 1e-12
    # The perturbation sum needs the whole manifold.
    @test_throws ArgumentError ac_zeeman_shift(
        Eigenstructure(ca43, "D_5/2 F=4"; B=B0),
        D_PROBE,
        Ω_RF,
        B_DRIVE,
    )

    # Polarisation of the transverse field: σ⁺ and σ⁻ rotating fields differ,
    # and a linear field is their mean (all at equal total |B|²).
    b_lin = [0.05u"mT", 0.0u"mT", 0.0u"mT"]
    b_plus = (0.05u"mT" / sqrt(2)) .* [1.0, -1.0im, 0.0]
    b_minus = (0.05u"mT" / sqrt(2)) .* [1.0, 1.0im, 0.0]
    δ_lin = ac_zeeman_shift(S_CLOCK, D_PROBE, Ω_RF, b_lin)
    δ_plus = ac_zeeman_shift(S_CLOCK, D_PROBE, Ω_RF, b_plus)
    δ_minus = ac_zeeman_shift(S_CLOCK, D_PROBE, Ω_RF, b_minus)
    @test !isapprox(δ_plus, δ_minus; rtol=1e-3)
    @test δ_lin ≈ (δ_plus + δ_minus) / 2 rtol = 1e-12

    # The phasor decomposition of zeeman_drives equals an explicitly built
    # two-drive representation.
    drives_manual = [
        zeeman_drive(S_CLOCK, real.(b_plus); phase=0.0),
        zeeman_drive(S_CLOCK, imag.(b_plus); phase=(-π / 2)),
    ]
    @test quasienergy_shift(S_CLOCK, D_PROBE, Ω_RF, drives_manual) ≈ δ_plus rtol = 1e-12

    # Contribution breakdown: entries carry adiabatic labels, sum to the total,
    # and are sorted by decreasing magnitude.
    contributions = quasienergy_shift_contributions(
        S_CLOCK,
        D_PROBE,
        Ω_RF,
        zeeman_drives(S_CLOCK, B_DRIVE),
    )
    @test sum(c.shift for c in contributions) ≈ r.upper rtol = 1e-12
    @test issorted([abs(c.shift) for c in contributions]; rev=true)
    @test contributions[1].state isa StateSpec{HyperfineNumberSpec}
    @test all(
        Levels.fine_structure(c.state.level) == NoHyperfineNumberSpec("D_5/2") for
        c in contributions
    )
end

@testitem "DC limit reproduces static curvature" tags=[:unit] setup=[Ca43ClockSetup] begin
    using FiniteDifferences

    # For Ω far below every internal splitting, the ac Zeeman shift of the
    # transition reduces to the static second-order response
    # δ = ¼ f″(B₀) b² of a field b cos(Ωt) along ẑ (the ⟨cos²⟩ = 1/2 average
    # of the parabolic response) — the approximation whose *failure* at trap
    # rf frequencies (D_5/2 hyperfine intervals of 5–25 MHz) motivates the
    # full treatment.
    f(b) = ustrip(
        u"Hz",
        cyclic(
            state_offset(clock_structure(b * u"mT"), D_PROBE) -
            state_offset(clock_structure(b * u"mT"), S_PROBE),
        ),
    )
    curvature = central_fdm(7, 2)(f, ustrip(u"mT", B0)) # Hz/mT²
    b_par = 0.05 # mT
    slow = ac_zeeman_shift(
        S_CLOCK,
        S_PROBE => D_PROBE,
        angular(100.0u"Hz"),
        [0.0u"mT", 0.0u"mT", b_par * u"mT"],
    )
    @test ustrip(u"Hz", cyclic(slow.shift)) ≈ curvature * b_par^2 / 4 rtol = 1e-3
end

@testitem "PT vs Floquet dressing" tags=[:unit] setup=[Ca43ClockSetup] begin
    # Second-order perturbation theory against the nonperturbative Floquet
    # dressing (dress_manifold), for drive frequencies away from the internal
    # resonances (S_1/2: intra-F Zeeman splittings ~1.7 MHz at B₀; D_5/2:
    # ΔF = ±1 hyperfine bands at ~5–25 MHz). rtol 2e-2 as in the reference
    # test suite; the residual is the next PT order.
    for (state, frequencies) in ((S_PROBE, (10.0, 30.0, 60.0)), (D_PROBE, (30.0, 60.0)))
        for f_mhz in frequencies
            ω = angular(f_mhz * u"MHz")
            pt = ac_zeeman_shift(S_CLOCK, state, ω, B_DRIVE)
            fl = ac_zeeman_shift(S_CLOCK, state, ω, B_DRIVE, Floquet())
            @test pt ≈ fl rtol = 2e-2
        end
    end

    # Harmonic-truncation convergence of the Floquet reference.
    ω = angular(30.0u"MHz")
    @test ac_zeeman_shift(S_CLOCK, D_PROBE, ω, B_DRIVE, Floquet(; harmonic_cutoff=4)) ≈
          ac_zeeman_shift(S_CLOCK, D_PROBE, ω, B_DRIVE, Floquet(; harmonic_cutoff=6)) rtol =
        1e-4
end

@testitem "Near-resonance warnings" tags=[:unit] setup=[Ca43ClockSetup] begin
    # Driving close to a D_5/2 internal splitting: the PT sum warns that its
    # result is unreliable there.
    Δ = abs(
        state_offset(S_CLOCK, StateSpec("D_5/2 F=5", 3)) -
        state_offset(S_CLOCK, D_PROBE),
    )
    ω = Δ + angular(50.0u"kHz")
    @test_logs (:warn, r"Near-resonant term") match_mode = :any ac_zeeman_shift(
        S_CLOCK,
        D_PROBE,
        ω,
        [0.05u"mT", 0.0u"mT", 0.0u"mT"],
    )

    # A strong drive on resonance spreads the bare state over several Floquet
    # modes; the dressed-state identification then warns as ambiguous.
    @test_logs (:warn, r"Ambiguous dressed-state") match_mode = :any ac_zeeman_shift(
        S_CLOCK,
        D_PROBE,
        Δ,
        [2.0u"mT", 0.0u"mT", 0.0u"mT"],
        Floquet(),
    )
end

@testitem "Hyperfine DrivenTransition" tags=[:unit] setup=[Ca43ClockSetup] begin
    using Levels.PeriodicDynamics: full_hamiltonian

    basis = Levels.basis(S_CLOCK)
    n, ε = beam_vectors(π / 2, π / 4)
    coupling = rabi_normalised(
        relative_couplings(S_CLOCK, "S_1/2" => "D_5/2", ε, n),
        basis,
        S_PROBE => D_PROBE,
        angular(10.0u"kHz"),
    )
    drives = zeeman_drives(S_CLOCK, B_DRIVE)
    dt = DrivenTransition(
        S_CLOCK,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        coupling,
        drives,
    )

    # The rotating frame vanishes on the probed pair and reproduces the
    # structure's eigen-energies (relative to the probed states) elsewhere.
    @test iszero(dt.frame[dt.lower])
    @test iszero(dt.frame[dt.upper])
    for i in dt.lower_range
        @test dt.frame[i] ≈
              state_offset(S_CLOCK, basis[i]) - state_offset(S_CLOCK, S_PROBE) atol =
            1e-9u"rad/µs"
    end

    # The model Hamiltonian is Hermitian, and the drive term seen by the model
    # is the eigenbasis Zeeman Hamiltonian of the structure.
    h = full_hamiltonian(dt, 0.01u"µs")
    @test h ≈ h'
    @test drives[1].amplitude == zeeman_hamiltonian(S_CLOCK, B_DRIVE)

    # The PT shift evaluated on the model's manifold blocks agrees exactly
    # with the structure-level wrapper (three-way consistency with the
    # reference values is covered in the PT item).
    wrapper = ac_zeeman_shift(S_CLOCK, S_PROBE => D_PROBE, Ω_RF, B_DRIVE)
    @test quasienergy_shift(dt, dt.upper) - quasienergy_shift(dt, dt.lower) ≈
          wrapper.shift rtol = 1e-12

    # The dressed carrier resonance shift agrees with the PT transition shift
    # up to the PT truncation error and the laser-induced ac Stark shift
    # (∝ Ω0², ~Hz here).
    sb = sideband_rabi(dt, DressedStates(; num_time_samples=32); sideband=0)
    @test abs(sb.δ_res - wrapper.shift) < angular(10.0u"Hz")

    # Constructor validation: same-manifold probe, wrong level kind, incomplete
    # manifolds, non-block-diagonal drives.
    @test_throws ArgumentError DrivenTransition(
        S_CLOCK,
        StateSpec("S_1/2 F=4", 4) => StateSpec("S_1/2 F=3", 3);
        drive_frequency=Ω_RF,
        coupling,
        drives,
    )
    @test_throws ArgumentError DrivenTransition(
        S_CLOCK,
        StateSpec("S_1/2", 1//2) => D_PROBE;
        drive_frequency=Ω_RF,
        coupling,
        drives,
    )
    incomplete = Eigenstructure(ca43, "S_1/2 F=4", "D_5/2"; B=B0)
    @test_throws ArgumentError DrivenTransition(
        incomplete,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        coupling=zeros(ComplexF64, 57, 57) .* u"rad/µs",
        drives=HarmonicDrive[],
    )
    bad_drive = HarmonicDrive(coupling .+ coupling', 0.0)
    @test_throws ArgumentError DrivenTransition(
        S_CLOCK,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        coupling,
        drives=[bad_drive],
    )
end

@testitem "sr88 PT cross-check vs Joshi closed form" tags=[:unit, :fast] setup=[
    PeriodicSetup,
] begin
    # Three-way check in the established ⁸⁸Sr⁺ scenario: the PT engine on the
    # DrivenTransition manifold blocks against the analytic transverse-field
    # ac Zeeman shift, Joshi et al., PRA 110, 063101 (2024), Eq. (14) — the
    # same target the dressed-state engine is validated against in
    # test-periodic-floquet.jl.
    B_perp = 5.0u"µT"
    dt = make_transition(
        TRANSITIONS[1];
        drives=[zeeman_drive(S0, [B_perp, 0.0u"µT", 0.0u"µT"])],
    )
    for (level, range) in (("S_1/2", dt.lower_range), ("D_5/2", dt.upper_range))
        g = lande_g(sr88, level)
        ω_z = g * μ_B * B0 / u"ħ"
        Ω_b = g * μ_B * B_perp / u"ħ"
        for (α, i) in enumerate(range)
            m = BASIS[i].m
            pred = -Ω_b^2 * (m / 8) * (1 / (Ω_RF - ω_z) - 1 / (Ω_RF + ω_z))
            @test quasienergy_shift(dt, i) ≈ pred rtol = 1e-6
        end
    end
end

@testitem "PT vs exact TDSE carrier resonance" tags=[:integration, :slow] setup=[
    Ca43ClockSetup,
] begin
    # End-to-end: the exact monodromy engine's carrier resonance detuning for
    # the hyperfine DrivenTransition agrees with the dressed-state engine (and
    # thus, to PT accuracy, with quasienergy_shift) — the ⁴³Ca⁺ analogue of the
    # sr88 cross-engine item in test-periodic-tdse.jl.
    n, ε = beam_vectors(π / 2, π / 4)
    coupling = rabi_normalised(
        relative_couplings(S_CLOCK, "S_1/2" => "D_5/2", ε, n),
        Levels.basis(S_CLOCK),
        S_PROBE => D_PROBE,
        angular(2.0u"kHz"),
    )
    dt = DrivenTransition(
        S_CLOCK,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        coupling,
        drives=zeeman_drives(S_CLOCK, B_DRIVE),
    )
    fl = sideband_rabi(dt, DressedStates(; num_time_samples=32); sideband=0)
    ex = sideband_rabi(
        dt,
        Monodromy(; estimate=DressedStates(; num_time_samples=32));
        sideband=0,
    )
    @test ex.Ω ≈ fl.Ω rtol = 1e-3
    # δ_res agrees between the engines up to the laser-induced ac Stark shift
    # (∝ Ω0², sub-Hz here).
    @test abs(ex.δ_res - fl.δ_res) < angular(5.0u"Hz")
end
