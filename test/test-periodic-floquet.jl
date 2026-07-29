@testitem "Pure EMM: Bessel sidebands" tags=[:unit] setup=[PeriodicSetup] begin
    using SpecialFunctions: besselj

    d = dress(make_transition(TRANSITIONS[1]))
    for (β_I, β_Q) in ((0.1, 0.0), (0.0, 0.08), (0.15, 0.06))
        β = hypot(β_I, β_Q)
        modulation = HarmonicPhaseModulation(β_I, β_Q)
        for sb in (+1, -1)
            r = sideband_amplitude(d; sideband=sb, modulation)
            @test uconvert(NoUnits, r.Ω / Ω0) ≈ abs(besselj(1, β)) rtol = 1e-10
            @test uconvert(NoUnits, r.carrier / Ω0) ≈ abs(besselj(0, β)) rtol = 1e-10
            @test r.δ_res ≈ sb * Ω_RF   # no ac fields: no resonance shift
        end
    end
end

@testitem "Longitudinal ac field: FM sidebands" tags=[:unit] setup=[PeriodicSetup] begin
    using SpecialFunctions: besselj

    B_par = 2.5u"µT"
    # Drive phase π/2 (capacitive load): the FM sidebands add coherently to
    # in-phase EMM and the sideband vanishes exactly at β_I = χ B_par / Ω_rf, for
    # both sidebands and all χ.
    drives = [zeeman_drive(sr88, BASIS, [0, 0, 1] .* B_par; phase=π / 2)]
    for tr in (TRANSITIONS[1], TRANSITIONS[5], TRANSITIONS[8]), sb in (+1, -1)
        d = dress(make_transition(tr; drives))
        β_B = uconvert(NoUnits, zeeman_sensitivity(sr88, tr) * B_par / Ω_RF)
        r0 = sideband_amplitude(
            d;
            sideband=sb,
            modulation=HarmonicPhaseModulation(β_B, 0.0),
        )
        @test uconvert(NoUnits, r0.Ω / Ω0) < 1e-9
        r1 = sideband_amplitude(d; sideband=sb)
        @test uconvert(NoUnits, r1.Ω / Ω0) ≈ abs(besselj(1, β_B)) rtol = 1e-6
    end

    # Drive phase 0: modulation in quadrature with in-phase EMM — hyperbolic floor.
    drives0 = [zeeman_drive(sr88, BASIS, [0, 0, 1] .* B_par; phase=0.0)]
    tr = TRANSITIONS[8]
    d = dress(make_transition(tr; drives=drives0))
    β_B = uconvert(NoUnits, zeeman_sensitivity(sr88, tr) * B_par / Ω_RF)
    for β_I in (-0.002, 0.0, 0.003)
        r = sideband_amplitude(
            d;
            sideband=+1,
            modulation=HarmonicPhaseModulation(β_I, 0.0),
        )
        @test uconvert(NoUnits, r.Ω / Ω0) ≈ abs(besselj(1, hypot(β_I, β_B))) rtol = 1e-4
    end
end

@testitem "Transverse ac field: ac Zeeman shifts" tags=[:unit] setup=[PeriodicSetup] begin
    # Quasienergy shifts of the dressed levels against Joshi et al. Eq. (14).
    B_perp = 2.0u"µT"
    drives = [zeeman_drive(sr88, BASIS, [1, 0, 0] .* B_perp; phase=π / 2)]
    dt = make_transition(StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2); drives)
    d = dress(dt; ngrid=16)

    for (level, r, ε) in
        (("S_1/2", dt.lower_range, d.ε_lower), ("D_5/2", dt.upper_range, d.ε_upper))
        g = lande_g(sr88, level)
        ω_z = g * μ_B * B0 / u"ħ"
        Ω_b = g * μ_B * B_perp / u"ħ"
        for (α, i) in enumerate(r)
            m = BASIS[i].m
            pred = -Ω_b^2 * (m / 8) * (1 / (Ω_RF - ω_z) - 1 / (Ω_RF + ω_z))
            @test ε[α] - dt.frame[i] ≈ pred rtol = 1e-3
        end
    end
end

@testitem "Transverse ac field: two-photon sidebands" tags=[:unit] setup=[PeriodicSetup] begin
    # Joshi et al. Eq. (15)/Table I structure, in the Ca⁺ QSIM-trap configuration:
    # the Δm = −3 two-photon transition |S,+1/2⟩ ↔ |D,−5/2⟩ has no direct E2
    # path, only mediated coupling. Each path amplitude and their destructive
    # combination must match Eq. (15), for both sidebands (whose mediator
    # detunings Ω_rf ± ω_z differ). Both paths involve one Δm = −1 rf photon, so
    # the transverse-field azimuth θ enters as a common phase only: the isolated
    # two-photon Rabi frequency is θ-independent.
    Ω_rf = 2π * 29.687u"MHz"
    b_static = 0.417u"mT"
    B_perp = 27.6u"µT"
    Ω_probe = 2π * 131.7u"kHz"
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", -5//2)
    C1, C2 = 0.515 / sqrt(3), 0.515   # laser coupling factors (Table I, RSB row)

    path1 = (-3//2, 1//2, C1)    # laser Δm = −2, then rf photon within D
    path2 = (-5//2, -1//2, C2)   # rf photon within S, then laser Δm = −2
    function coupling_for(paths...)
        L = zeros(ComplexF64, length(BASIS), length(BASIS)) .* u"kHz"
        for (m_upper, m_lower, c) in paths
            L[
                stateindex(BASIS, "D_5/2", m_upper),
                stateindex(BASIS, "S_1/2", m_lower),
            ] = c * Ω_probe
        end
        L
    end

    function twophoton(coupling, sb; θ=0.0)
        drives = [zeeman_drive(sr88, BASIS, [cos(θ), sin(θ), 0] .* B_perp; phase=π / 2)]
        dt = DrivenTransition(
            sr88,
            BASIS,
            probe;
            drive_frequency=Ω_rf,
            static_field=b_static,
            coupling,
            drives,
        )
        sideband_amplitude(dress(dt); sideband=sb).Ω
    end

    ω_zS = lande_g(sr88, "S_1/2") * μ_B * b_static / u"ħ"
    ω_zD = lande_g(sr88, "D_5/2") * μ_B * b_static / u"ħ"
    Ω_b1 = lande_g(sr88, "D_5/2") * μ_B * B_perp / u"ħ" * sqrt(5) / 2  # ⟨−5/2|J₋|−3/2⟩ = √5
    Ω_b2 = lande_g(sr88, "S_1/2") * μ_B * B_perp / u"ħ" / 2            # ⟨−1/2|J₋|+1/2⟩ = 1
    for (sb, s) in ((-1, +1), (+1, -1))   # RSB: mediator detuning Ω_rf + ω_z; BSB: −
        t1 = uconvert(u"kHz", C1 * Ω_probe * Ω_b1 / (2 * (Ω_rf + s * ω_zD)))
        t2 = uconvert(u"kHz", C2 * Ω_probe * Ω_b2 / (2 * (Ω_rf + s * ω_zS)))
        @test twophoton(coupling_for(path1), sb) ≈ t1 rtol = 1e-3
        @test twophoton(coupling_for(path2), sb) ≈ t2 rtol = 1e-3
        both = twophoton(coupling_for(path1, path2), sb)
        @test both ≈ abs(t1 - t2) rtol = 1e-2
        @test twophoton(coupling_for(path1, path2), sb; θ=1.1) ≈ both rtol = 1e-6
    end
end
