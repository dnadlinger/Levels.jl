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
    const Ω_RF = 2π * 30.0u"MHz"
    const S_PROBE = StateSpec("S_1/2 F=4", 4)
    const D_PROBE = StateSpec("D_5/2 F=4", 3)
    # B_⊥ x̂ + B_∥ ẑ, 50 µT each (the reference README scenario).
    const B_DRIVE = [0.05u"mT", 0.0u"mT", 0.05u"mT"]

    s_manifold(b=B0) = hyperfine_manifold(ca43, "S_1/2", b)
    d_manifold(b=B0) = hyperfine_manifold(ca43, "D_5/2", b)
end

@testitem "PT ac Zeeman shift vs reference implementation" tags=[:unit, :fast] setup=[
    Ca43ClockSetup,
] begin
    hz(x) = ustrip(u"Hz", x / 2π)

    # The stretched S(4,+4) state is an exact product state, so a pure π
    # (B_∥) drive shifts it by exactly zero…
    @test iszero(
        ac_zeeman_shift(s_manifold(), S_PROBE, Ω_RF, [0.0u"mT", 0.0u"mT", 0.05u"mT"]),
    )
    # …and its whole shift comes from the transverse component.
    r = ac_zeeman_shift(ca43, S_PROBE => D_PROBE, B0, Ω_RF, B_DRIVE)
    @test ac_zeeman_shift(
        s_manifold(),
        S_PROBE,
        Ω_RF,
        [0.05u"mT", 0.0u"mT", 0.0u"mT"],
    ) ≈ r.lower rtol = 1e-12

    # No π/σ cross terms at second order: the ∥ and ⊥ shifts add exactly, and
    # the relative temporal phase of the two components is irrelevant.
    m_d = d_manifold()
    δ_par = ac_zeeman_shift(m_d, D_PROBE, Ω_RF, [0.0u"mT", 0.0u"mT", 0.05u"mT"])
    δ_perp = ac_zeeman_shift(m_d, D_PROBE, Ω_RF, [0.05u"mT", 0.0u"mT", 0.0u"mT"])
    @test δ_par + δ_perp ≈ r.upper rtol = 1e-12
    drives_dephased = [
        zeeman_drive(ca43, m_d.basis, [0.05u"mT", 0.0u"mT", 0.0u"mT"]; phase=0.0),
        zeeman_drive(ca43, m_d.basis, [0.0u"mT", 0.0u"mT", 0.05u"mT"]; phase=1.234),
    ]
    @test ac_zeeman_shift(m_d, D_PROBE, Ω_RF, drives_dephased) ≈ r.upper rtol = 1e-12

    # Polarisation of the transverse field: σ⁺ and σ⁻ rotating fields differ,
    # and a linear field is their mean (all at equal total |B|²).
    b_lin = [0.05u"mT", 0.0u"mT", 0.0u"mT"]
    b_plus = (0.05u"mT" / sqrt(2)) .* [1.0, -1.0im, 0.0]
    b_minus = (0.05u"mT" / sqrt(2)) .* [1.0, 1.0im, 0.0]
    δ_lin = ac_zeeman_shift(m_d, D_PROBE, Ω_RF, b_lin)
    δ_plus = ac_zeeman_shift(m_d, D_PROBE, Ω_RF, b_plus)
    δ_minus = ac_zeeman_shift(m_d, D_PROBE, Ω_RF, b_minus)
    @test !isapprox(δ_plus, δ_minus; rtol=1e-3)
    @test δ_lin ≈ (δ_plus + δ_minus) / 2 rtol = 1e-12

    # The phasor decomposition of zeeman_drives equals an explicitly built
    # two-drive representation.
    drives_manual = [
        zeeman_drive(ca43, m_d.basis, real.(b_plus); phase=0.0),
        zeeman_drive(ca43, m_d.basis, imag.(b_plus); phase=(-π / 2)),
    ]
    @test ac_zeeman_shift(m_d, D_PROBE, Ω_RF, drives_manual) ≈ δ_plus rtol = 1e-12

    # Contribution breakdown: entries carry adiabatic labels, sum to the total,
    # and are sorted by decreasing magnitude.
    contributions = ac_zeeman_contributions(m_d, D_PROBE, Ω_RF, B_DRIVE)
    @test sum(c.shift for c in contributions) ≈ r.upper rtol = 1e-12
    @test issorted([abs(c.shift) for c in contributions]; rev=true)
    @test contributions[1].state isa StateSpec{HyperfineNumberSpec}
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
        (
            state_energy(d_manifold(b * u"mT"), D_PROBE) -
            state_energy(s_manifold(b * u"mT"), S_PROBE)
        ) / 2π,
    )
    curvature = central_fdm(7, 2)(f, ustrip(u"mT", B0)) # Hz/mT²
    b_par = 0.05 # mT
    slow = ac_zeeman_shift(
        ca43,
        S_PROBE => D_PROBE,
        B0,
        2π * 100.0u"Hz",
        [0.0u"mT", 0.0u"mT", b_par * u"mT"],
    )
    @test ustrip(u"Hz", slow.shift / 2π) ≈ curvature * b_par^2 / 4 rtol = 1e-3
end

@testitem "PT vs Floquet dressing" tags=[:unit] setup=[Ca43ClockSetup] begin
    # Second-order perturbation theory against the nonperturbative Floquet
    # dressing (dress_manifold), for drive frequencies away from the internal
    # resonances (S_1/2: intra-F Zeeman splittings ~1.7 MHz at B₀; D_5/2:
    # ΔF = ±1 hyperfine bands at ~5–25 MHz). rtol 2e-2 as in the reference
    # test suite; the residual is the next PT order.
    for (manifold, state, frequencies) in (
        (s_manifold(), S_PROBE, (10.0, 30.0, 60.0)),
        (d_manifold(), D_PROBE, (30.0, 60.0)),
    )
        for f_mhz in frequencies
            ω = 2π * f_mhz * u"MHz"
            pt = ac_zeeman_shift(manifold, state, ω, B_DRIVE)
            fl = floquet_zeeman_shift(manifold, state, ω, B_DRIVE)
            @test pt ≈ fl rtol = 2e-2
        end
    end

    # Harmonic-truncation convergence of the Floquet reference.
    m_d = d_manifold()
    ω = 2π * 30.0u"MHz"
    @test floquet_zeeman_shift(m_d, D_PROBE, ω, B_DRIVE; nharm=4) ≈
          floquet_zeeman_shift(m_d, D_PROBE, ω, B_DRIVE; nharm=6) rtol = 1e-4
end

@testitem "Near-resonance warnings" tags=[:unit] setup=[Ca43ClockSetup] begin
    # Driving close to a D_5/2 internal splitting: the PT sum warns that its
    # result is unreliable there.
    m_d = d_manifold()
    Δ = abs(state_energy(m_d, StateSpec("D_5/2 F=5", 3)) - state_energy(m_d, D_PROBE))
    ω = Δ + 2π * 50.0u"kHz"
    @test_logs (:warn, r"Near-resonant term") match_mode = :any ac_zeeman_shift(
        m_d,
        D_PROBE,
        ω,
        [0.05u"mT", 0.0u"mT", 0.0u"mT"],
    )

    # A strong drive on resonance spreads the bare state over several Floquet
    # modes; the dressed-state identification then warns as ambiguous.
    @test_logs (:warn, r"Ambiguous dressed-state") match_mode = :any floquet_zeeman_shift(
        m_d,
        D_PROBE,
        Δ,
        [2.0u"mT", 0.0u"mT", 0.0u"mT"],
    )
end

@testitem "Hyperfine DrivenTransition" tags=[:unit] setup=[Ca43ClockSetup] begin
    using Levels.PeriodicDynamics: full_hamiltonian

    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    n, ε = beam_vectors(π / 2, π / 4)
    coupling = rabi_normalised(
        quadrupole_couplings(ca43, basis, "S_1/2", "D_5/2", ε, n),
        basis,
        S_PROBE => D_PROBE,
        2π * 10.0u"kHz",
    )
    drives = zeeman_drives(ca43, basis, B_DRIVE)
    dt = DrivenTransition(
        ca43,
        basis,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        static_field=B0,
        coupling,
        drives,
    )

    # The rotating frame vanishes on the probed pair and reproduces the
    # manifold eigen-energies (relative to the probed states) elsewhere.
    @test iszero(dt.frame[dt.lower])
    @test iszero(dt.frame[dt.upper])
    m_s = s_manifold()
    for i in dt.lower_range
        @test dt.frame[i] ≈ state_energy(m_s, basis[i]) - state_energy(m_s, S_PROBE) atol =
            1e-9u"µs^-1"
    end

    # The model Hamiltonian is Hermitian, and the drive term seen by the model
    # equals the eigenbasis-rotated coupled-basis drive.
    h = full_hamiltonian(dt, 0.01u"µs")
    @test h ≈ h'

    # The PT shift evaluated on the model's manifold blocks agrees exactly
    # with the manifold-level wrapper (three-way consistency with the
    # reference values is covered in the PT item).
    wrapper = ac_zeeman_shift(ca43, S_PROBE => D_PROBE, B0, Ω_RF, B_DRIVE)
    @test ac_zeeman_shift(dt, dt.upper) - ac_zeeman_shift(dt, dt.lower) ≈ wrapper.shift rtol =
        1e-12

    # The dressed carrier resonance shift agrees with the PT transition shift
    # up to the PT truncation error and the laser-induced ac Stark shift
    # (∝ Ω0², ~Hz here).
    sb = sideband_rabi(dt; sideband=0, ngrid=32)
    @test abs(sb.δ_res - wrapper.shift) < 2π * 10.0u"Hz"

    # Constructor validation.
    @test_throws ArgumentError DrivenTransition(
        ca43,
        basis,
        StateSpec("S_1/2 F=4", 4) => StateSpec("S_1/2 F=3", 3);
        drive_frequency=Ω_RF,
        static_field=B0,
        coupling,
        drives,
    )
    @test_throws ArgumentError DrivenTransition(
        ca43,
        basis,
        StateSpec("S_1/2", 1//2) => D_PROBE;
        drive_frequency=Ω_RF,
        static_field=B0,
        coupling,
        drives,
    )
    incomplete = StateBasis(ca43, "S_1/2 F=4", "D_5/2")
    @test_throws ArgumentError DrivenTransition(
        ca43,
        incomplete,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        static_field=B0,
        coupling=zeros(ComplexF64, 57, 57) .* u"µs^-1",
        drives=HarmonicDrive[],
    )
    bad_drive = HarmonicDrive(coupling .+ coupling', 0.0)
    @test_throws ArgumentError DrivenTransition(
        ca43,
        basis,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        static_field=B0,
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
        drives=[zeeman_drive(sr88, BASIS, [B_perp, 0.0u"µT", 0.0u"µT"])],
    )
    for (level, range) in (("S_1/2", dt.lower_range), ("D_5/2", dt.upper_range))
        g = lande_g(sr88, level)
        ω_z = g * μ_B * B0 / u"ħ"
        Ω_b = g * μ_B * B_perp / u"ħ"
        for (α, i) in enumerate(range)
            m = BASIS[i].m
            pred = -Ω_b^2 * (m / 8) * (1 / (Ω_RF - ω_z) - 1 / (Ω_RF + ω_z))
            @test ac_zeeman_shift(dt, i) ≈ pred rtol = 1e-6
        end
    end
end

@testitem "PT vs exact TDSE carrier resonance" tags=[:integration, :slow] setup=[
    Ca43ClockSetup,
] begin
    # End-to-end: the exact monodromy engine's carrier resonance detuning for
    # the hyperfine DrivenTransition agrees with the dressed-state engine (and
    # thus, to PT accuracy, with ac_zeeman_shift) — the ⁴³Ca⁺ analogue of the
    # sr88 cross-engine item in test-periodic-tdse.jl.
    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    n, ε = beam_vectors(π / 2, π / 4)
    coupling = rabi_normalised(
        quadrupole_couplings(ca43, basis, "S_1/2", "D_5/2", ε, n),
        basis,
        S_PROBE => D_PROBE,
        2π * 2.0u"kHz",
    )
    dt = DrivenTransition(
        ca43,
        basis,
        S_PROBE => D_PROBE;
        drive_frequency=Ω_RF,
        static_field=B0,
        coupling,
        drives=zeeman_drives(ca43, basis, B_DRIVE),
    )
    fl = sideband_rabi(dt; sideband=0, ngrid=32)
    ex = exact_sideband(dt; sideband=0, ngrid=32)
    @test ex.Ω ≈ fl.Ω rtol = 1e-3
    # δ_res agrees between the engines up to the laser-induced ac Stark shift
    # (∝ Ω0², sub-Hz here).
    @test abs(ex.δ_res - fl.δ_res) < 2π * 5.0u"Hz"
end
