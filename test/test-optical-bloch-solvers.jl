# Solver-side tests for Levels.OpticalBloch through the LevelsQuantumToolboxExt
# extension: internal steady states (two-level saturation, Λ dark resonance,
# hyperfine ≡ fine-structure ⊗ nuclear spin), adiabatic-elimination cooling
# rates against the two-level sideband-cooling formula and against the full
# internal ⊗ motional model, and the phonon-number observables.

@testitem "Two-level steady state and populations" tags=[:integration] setup=[
    OpticalBlochSetup,
] begin
    using QuantumToolbox: liouvillian, steadystate

    # |S₁/₂, −½⟩ ↔ |P₁/₂, +½⟩ under σ⁺ light along ẑ, the other S state left out
    # of the basis, is a two-level system with Γ = (2/3) A (the σ⁺ decay
    # channel) — its steady state is the saturation formula exactly.
    basis = StateBasis([StateSpec("S_1/2", -1//2), StateSpec("P_1/2", 1//2)])
    n_z, ε_σ = beam_vectors(0.0, π / 4, π / 2)
    Γ = 2 / 3 * einstein_a(sr88, "S_1/2", "P_1/2")
    I_sat = saturation_intensity(sr88, "S_1/2", "P_1/2")
    for (s, δ_over_Γ) in ((0.5, 0.0), (2.0, 1.5), (7.0, -3.0))
        beam = LaserBeam("S_1/2" => "P_1/2", δ_over_Γ * Γ, s * I_sat, ε_σ, n_z)
        scheme =
            LaserScheme(basis, [beam]; static_field=0.0u"mT", frame_reference="S_1/2")
        model = lindblad_model(sr88, scheme; open_channels=:drop)
        ρ = steadystate(model)
        Ω = 2 * abs(model.hamiltonian[2, 1])
        # In the frame of the lower level, the upper state sits at minus the
        # laser detuning.
        δ = real(model.hamiltonian[1, 1] - model.hamiltonian[2, 2])
        @test δ ≈ δ_over_Γ * Γ
        s_eff = 2 * Ω^2 / Γ^2
        @test populations(ρ, model, "P_1/2") ≈ (s_eff / 2) / (1 + s_eff + (2δ / Γ)^2) rtol =
            1e-9
        @test sum(populations(ρ, model)) ≈ 1 rtol = 1e-12
        @test populations(ρ, model, "S_1/2") + populations(ρ, model, "P_1/2") ≈ 1 rtol =
            1e-12
        @test size(liouvillian(model).data) == (4, 4)
    end
    # A density matrix of the wrong size is rejected.
    beam = LaserBeam("S_1/2" => "P_1/2", 0.0u"µs^-1", I_sat, ε_σ, n_z)
    two = lindblad_model(
        sr88,
        LaserScheme(basis, [beam]; static_field=0.0u"mT", frame_reference="S_1/2");
        open_channels=:drop,
    )
    eight =
        lindblad_model(CA40_NONREL_G, lambda_scheme(2π * 300.0u"MHz", 2π * 290.0u"MHz"))
    @test_throws ArgumentError populations(steadystate(two), eight)
end

@testitem "Λ dark resonance" tags=[:integration] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: steadystate

    # Golden-section minimiser over a bracket (the test project carries no
    # optimisation package).
    function golden_min(f, a, b; tol=1e-5)
        φ = (sqrt(5) - 1) / 2
        c, d = b - φ * (b - a), a + φ * (b - a)
        fc, fd = f(c), f(d)
        while b - a > tol
            if fc < fd
                b, d, fd = d, c, fc
                c = b - φ * (b - a)
                fc = f(c)
            else
                a, c, fc = c, d, fd
                d = a + φ * (b - a)
                fd = f(d)
            end
        end
        (a + b) / 2
    end

    # In the large-detuning limit the excited population dips where the
    # two-photon (Raman) resonance condition of one Λ subsystem is met, i.e. for
    # |S, m_s⟩–|P, m_p⟩–|D, m_d⟩ at Δ866 = Δ397 + Δω_s − Δω_d up to the light
    # shift of the excited state by the 397 nm beam. The lowest-lying dip is the
    # |S, −½⟩–|P, +½⟩–|D, +3/2⟩ system.
    Δ397 = 2π * 300.0u"MHz"
    u = uconvert(u"µs^-1", Levels.BOHR_MAGNETON * B_LAMBDA / u"ħ")
    function p_population(Δ866)
        model = lindblad_model(CA40_NONREL_G, lambda_scheme(Δ397, Δ866))
        populations(steadystate(model), model, "P_1/2")
    end
    analytic = Δ397 + 2 * (-1 / 2) * u - 4 / 5 * (3 / 2) * u   # Zeeman shifts g m u
    @test analytic / 2π ≈ 287.923u"MHz" atol = 0.002u"MHz"
    dip = golden_min(x -> p_population(2π * x * u"MHz"), 287.5, 288.4)
    # The 397 nm light shift of the excited state moves the resonance up by
    # ~0.1 MHz from the bare two-photon condition.
    @test 287.92 < dip < 288.1
    @test dip > ustrip(u"MHz", analytic / 2π)
    # Contrast: the dip is deep against the off-resonant background.
    background = p_population(2π * 285.0u"MHz")
    @test p_population(2π * dip * u"MHz") < 0.7 * background
    @test background ≈ 5.0e-4 rtol = 0.15
end

@testitem "Hyperfine steady state reduces to the fine-structure one" tags=[:integration] setup=[
    OpticalBlochSetup,
] begin
    using QuantumToolbox: steadystate

    # Zero hyperfine constants: the two m_I sectors are decoupled copies of the
    # fine-structure dynamics (shifted by the nuclear Zeeman energy, a constant
    # within each sector), so the level populations coincide.
    B = 0.5u"mT"
    n, ε = beam_vectors(deg2rad(60.0), 0.4, 0.2)
    I0 = 3.0u"W/m^2"
    Δ = 2π * 5.0u"MHz"
    hf_basis = StateBasis(SR88_TOY_HF, "S_1/2", "P_1/2", "D_3/2")
    fs_basis = StateBasis(["S_1/2", "P_1/2", "D_3/2"])
    hf = lindblad_model(
        SR88_TOY_HF,
        LaserScheme(
            hf_basis,
            [
                LaserBeam("S_1/2 F=1" => "P_1/2 F=1", Δ, I0, ε, n),
                LaserBeam("D_3/2 F=2" => "P_1/2 F=1", -Δ, 4I0, ε, n),
            ];
            static_field=B,
        ),
    )
    fs = lindblad_model(
        sr88,
        LaserScheme(
            fs_basis,
            [
                LaserBeam("S_1/2" => "P_1/2", Δ, I0, ε, n),
                LaserBeam("D_3/2" => "P_1/2", -Δ, 4I0, ε, n),
            ];
            static_field=B,
        ),
    )
    ρ_hf = steadystate(hf)
    ρ_fs = steadystate(fs)
    for level in ("S_1/2", "P_1/2", "D_3/2")
        @test populations(ρ_hf, hf, level) ≈ populations(ρ_fs, fs, level) rtol = 1e-8
    end
end

@testitem "Adiabatic-elimination cooling rates: two-level sideband cooling" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: steadystate, ptrace

    # Resolved-sideband regime of the two-level toy (Γ = 2A/3 for the σ⁺ channel
    # kept): the rates equal the weak-coupling scattering rates on the red and
    # blue sidebands, η² Γ Ω² / (Γ² + 4 (δ ± ω_m)²), plus the recoil diffusion
    # α η₀² Γ ρ_ee on both — red detuning cools, blue heats — and the full
    # motional model reproduces n̄ and the relaxation time.
    basis = StateBasis([StateSpec("S_1/2", -1//2), StateSpec("P_1/2", 1//2)])
    n_z, ε_σ = beam_vectors(0.0, π / 4, π / 2)
    Γ = 2 / 3 * einstein_a(sr88, "S_1/2", "P_1/2")
    I_sat = saturation_intensity(sr88, "S_1/2", "P_1/2")
    ω_m = 2π * 20.0u"MHz"
    mode = MotionalMode(ω_m, [0, 0, 1.0])
    for (s, δ_over_ω) in ((0.01, -1.0), (0.01, -0.7), (0.05, -1.0), (0.01, 1.0))
        δ = δ_over_ω * ω_m
        beam = LaserBeam("S_1/2" => "P_1/2", δ, s * I_sat, ε_σ, n_z)
        scheme =
            LaserScheme(basis, [beam]; static_field=0.0u"mT", frame_reference="S_1/2")
        model = lindblad_model(sr88, scheme; open_channels=:drop)
        mc = motional_coupling(sr88, scheme, model, mode)
        ρ = steadystate(model)
        rates = cooling_rates(model, mc; ρ)
        Ω = 2 * abs(model.hamiltonian[2, 1])
        η, η0, α =
            mc.projected_lamb_dicke[1], mc.recoil_lamb_dicke[1], mc.recoil_moments[1]
        R(d) = Γ * Ω^2 / (Γ^2 + 4d^2)
        D = α * η0^2 * Γ * populations(ρ, model, "P_1/2")
        A_minus = η^2 * R(δ + ω_m) + D
        A_plus = η^2 * R(δ - ω_m) + D
        @test rates.A_minus ≈ A_minus rtol = 3e-2 * (1 + 10s)
        @test rates.A_plus ≈ A_plus rtol = 3e-2 * (1 + 10s)
        if δ < zero(δ)
            @test rates.A_minus > rates.A_plus
            @test rates.nbar ≈ A_plus / (A_minus - A_plus) rtol = 3e-2
            @test rates.τ_c ≈ 1 / (A_minus - A_plus) rtol = 3e-2
            # Background heating adds to the steady-state occupation.
            heated = cooling_rates(model, mc; ρ, heating_rate=0.1 * rates.A_minus)
            @test heated.nbar ≈
                  rates.nbar + 0.1 * rates.A_minus / (rates.A_minus - rates.A_plus) rtol =
                1e-9
            @test heated.τ_c == rates.τ_c
        else
            @test rates.A_plus > rates.A_minus
            @test isnan(rates.nbar) && isnan(ustrip(rates.τ_c))
        end
        if s == 0.01 && δ_over_ω == -1.0
            for n_max in (4, 7)
                full = motional_model(model, mc; n_max)
                ρ_full = steadystate(full.H, full.c_ops)
                @test mean_phonon_number(ρ_full, model, mc) ≈ rates.nbar rtol = 5e-3
                @test mean_phonon_number(ρ_full, model, mc; method=:thermal_ratio) ≈
                      rates.nbar rtol = 2e-3
                @test populations(ptrace(ρ_full, 1), model, "P_1/2") ≈
                      populations(ρ, model, "P_1/2") rtol = 1e-3
                # The slowest phonon-number relaxation is the rate-equation τ_c,
                # not the twice slower coherence decay.
                @test cooling_time(full.H, full.c_ops) ≈ rates.τ_c rtol = 2e-2
            end
        end
    end
    # Vector form: one steady state shared by several modes.
    beam = LaserBeam("S_1/2" => "P_1/2", -ω_m, 0.01I_sat, ε_σ, n_z)
    scheme = LaserScheme(basis, [beam]; static_field=0.0u"mT", frame_reference="S_1/2")
    model = lindblad_model(sr88, scheme; open_channels=:drop)
    modes = [mode, MotionalMode(2π * 5.0u"MHz", [1.0, 0, 1.0])]
    mcs = motional_coupling(sr88, scheme, model, modes)
    many = cooling_rates(model, mcs)
    @test length(many) == 2
    @test many[1] == cooling_rates(model, mcs[1])
    @test many[2].A_minus ≈ cooling_rates(model, mcs[2]).A_minus
end

@testitem "Λ cooling: adiabatic elimination vs full model" tags=[:integration, :slow] setup=[
    OpticalBlochSetup,
] begin
    using QuantumToolbox: steadystate, ptrace

    # Dark-resonance cooling on the |S, −½⟩–|P, +½⟩–|D, +3/2⟩ Λ system at an
    # 866 nm detuning one trap frequency above its dark resonance, at reduced
    # intensity so that η Ω ≲ ω_m/4 keeps the adiabatic elimination roughly
    # valid. Even so, the full model's phonon distribution is not quite thermal:
    # p_{n+1}/p_n rises from 0.73 to 0.77 over the first dozen Fock states (rate
    # model: 0.715 throughout), and the converged ⟨n⟩ = 3.05 (n_max = 32) lies
    # 22 % above the rate model's n̄ = 2.51. The truncation at 8 perturbs only the
    # top ratio, so the thermal-ratio estimate (2.77; 2.86 converged) is compared
    # instead, while the direct ⟨n⟩ (2.04) is truncation-limited. The internal
    # populations agree to ~10 %. (The relaxation time is not compared: at
    # n̄ ≈ 3 the slowest mode of a Fock space truncated at 8 is set by the
    # truncation; the two-level item above pins cooling_time.)
    ω_m = 2π * 1.1u"MHz"
    mode = MotionalMode(ω_m, [0, 0, 1.0])
    scheme =
        lambda_scheme(2π * 300.0u"MHz", 2π * 289.0u"MHz"; P397=1.0u"µW", P866=5.0u"µW")
    model = lindblad_model(CA40_NONREL_G, scheme)
    mc = motional_coupling(CA40_NONREL_G, scheme, model, mode)
    ρ = steadystate(model)
    rates = cooling_rates(model, mc; ρ)
    @test rates.A_minus > rates.A_plus
    @test abs(mc.projected_lamb_dicke[2]) *
          2 *
          maximum(abs, model.hamiltonian[7:8, 3:6]) < 0.3ω_m
    full = motional_model(model, mc; n_max=8)
    ρ_full = steadystate(full.H, full.c_ops)
    nbar_full = mean_phonon_number(ρ_full, model, mc; method=:thermal_ratio)
    @test nbar_full ≈ rates.nbar rtol = 0.15
    @test mean_phonon_number(ρ_full, model, mc) < 0.8nbar_full
    @test_throws ArgumentError mean_phonon_number(ρ_full, model, mc; method=:thermal)
    @test populations(ptrace(ρ_full, 1), model, "P_1/2") ≈
          populations(ρ, model, "P_1/2") rtol = 0.15
    @test cooling_time(full.H, full.c_ops) > 1000u"µs"
    # Two modes in one full model: each recovers its own occupation (the
    # second, off-resonant mode is heated, i.e. its occupation is set by the
    # truncation and the rate model reports NaN).
    modes = [mode, MotionalMode(2π * 2.3u"MHz", [0, 0, 1.0])]
    mcs = motional_coupling(CA40_NONREL_G, scheme, model, modes)
    two = motional_model(model, mcs; n_max=[6, 3])
    ρ_two = steadystate(two.H, two.c_ops)
    nbars = mean_phonon_number(ρ_two, model, mcs; method=:thermal_ratio)
    both = cooling_rates(model, mcs; ρ)
    @test length(nbars) == 2
    @test nbars[1] ≈ both[1].nbar rtol = 0.15
    @test isnan(both[2].nbar) && both[2].A_plus > both[2].A_minus
    @test_throws ArgumentError motional_model(model, mcs; n_max=[6])
    @test_throws ArgumentError motional_model(model, mc; n_max=1)
end
