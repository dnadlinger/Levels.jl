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
            for num_fock in (4, 7)
                full = motional_model(model, mc; num_fock)
                @test full isa MotionalModel
                @test full.num_fock == [num_fock] && full.lamb_dicke_order == 1
                ρ_full = steadystate(full)
                @test mean_phonon_number(ρ_full, full) ≈ rates.nbar rtol = 5e-3
                @test mean_phonon_number(ρ_full, full; method=:thermal_ratio) ≈
                      rates.nbar rtol = 2e-3
                @test populations(ptrace(ρ_full, 1), model, "P_1/2") ≈
                      populations(ρ, model, "P_1/2") rtol = 1e-3
                # The slowest phonon-number relaxation is the rate-equation τ_c,
                # not the twice slower coherence decay.
                @test cooling_time(full) ≈ rates.τ_c rtol = 2e-2
                @test cooling_time(full.H, full.c_ops) == cooling_time(full)
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
    # model: 0.715 throughout), and the converged ⟨n⟩ = 3.05 (32 Fock states) lies
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
    full = motional_model(model, mc; num_fock=8)
    ρ_full = steadystate(full)
    nbar_full = mean_phonon_number(ρ_full, full; method=:thermal_ratio)
    @test nbar_full ≈ rates.nbar rtol = 0.15
    @test mean_phonon_number(ρ_full, full) < 0.8nbar_full
    @test_throws ArgumentError mean_phonon_number(ρ_full, full; method=:thermal)
    @test populations(ptrace(ρ_full, 1), model, "P_1/2") ≈
          populations(ρ, model, "P_1/2") rtol = 0.15
    @test cooling_time(full) > 1000u"µs"
    # Two modes in one full model: each recovers its own occupation (the
    # second, off-resonant mode is heated, i.e. its occupation is set by the
    # truncation and the rate model reports NaN).
    modes = [mode, MotionalMode(2π * 2.3u"MHz", [0, 0, 1.0])]
    mcs = motional_coupling(CA40_NONREL_G, scheme, model, modes)
    two = motional_model(model, mcs; num_fock=[6, 3])
    ρ_two = steadystate(two)
    nbars = mean_phonon_number(ρ_two, two; method=:thermal_ratio)
    both = cooling_rates(model, mcs; ρ)
    @test length(nbars) == 2
    @test nbars[1] ≈ both[1].nbar rtol = 0.15
    @test isnan(both[2].nbar) && both[2].A_plus > both[2].A_minus
    @test_throws ArgumentError motional_model(model, mcs; num_fock=[6])
    @test_throws ArgumentError motional_model(model, mc; num_fock=1)
    @test_throws ArgumentError motional_model(model, mc; num_fock=4, lamb_dicke_order=0)
    @test_throws ArgumentError motional_model(
        model,
        mc;
        num_fock=4,
        lamb_dicke_order=1.5,
    )
end

@testitem "Motional model beyond first order in the Lamb–Dicke parameters" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: liouvillian, steadystate, destroy, expect
    using SparseArrays

    ω = 2π * 0.6u"MHz"
    γ = 2π * 50.0u"kHz"
    Ω = 2π * 200.0u"kHz"
    t = sideband_toy(; γ, ω, η=0.1, Ω, δ=(-ω))
    @test abs(2 * t.model.hamiltonian[2, 1]) ≈ Ω
    @test t.mc.emission == :isotropic && t.mc.recoil_moments == [1 / 3]
    num_fock = 10
    first_order = motional_model(t.model, t.mc; num_fock)
    exact = motional_model(t.model, t.mc; num_fock, lamb_dicke_order=Inf)
    second = motional_model(t.model, t.mc; num_fock, lamb_dicke_order=2)
    @test first_order.lamb_dicke_order == 1 && isinf(exact.lamb_dicke_order)
    @test first_order.num_states == 2 && first_order.num_fock == [num_fock]
    @test first_order.modes == [t.mode]
    # Labels: order 1 carries L ⊗ 1 and one recoil operator per decay and mode;
    # the exact model the quadrature set (12 nodes by default, k + 1 at order k).
    @test first_order.c_labels[1] isa DecayLabel
    @test first_order.c_labels[2] == RecoilLabel(first_order.c_labels[1], 1, 0)
    @test length(exact.c_ops) == 12 &&
          all(l -> l isa RecoilLabel && l.mode == 0, exact.c_labels)
    @test [l.node for l in exact.c_labels] == 1:12
    @test length(second.c_ops) == 3
    @test length(
        motional_model(t.model, t.mc; num_fock, lamb_dicke_order=Inf, nodes=5).c_ops,
    ) == 5

    # First order: H = H_int ⊗ 1 + ω a†a + H_sb ⊗ (a + a†), the previous model.
    a = Matrix(destroy(num_fock).data)
    x = a + a'
    H_int = strip_h(t.model.hamiltonian)
    H_sb = strip_h(t.mc.sideband_hamiltonian)
    expected =
        kron(H_int, I(num_fock)) +
        ustrip(u"µs^-1", ω) * kron(I(2), a' * a) +
        kron(H_sb, x)
    @test Matrix(first_order.H.data) ≈ expected atol = 1e-12
    R = ustrip.(u"µs^(-1/2)", t.mc.recoil_operators[1])
    @test Matrix(first_order.c_ops[2].data) ≈ kron(R, x) atol = 1e-12
    # Exact: the beam coupling carries the displacement elements, each
    # quadrature operator a displaced decay operator with the node's weight.
    C = strip_h(t.model.couplings[1])
    η = t.mc.projected_lamb_dicke[1]
    D = displacement_elements(η, num_fock)
    T = kron(C, D) / 2
    @test Matrix(exact.H.data) ≈
          kron(Diagonal(diag(H_int)), I(num_fock)) +
          ustrip(u"µs^-1", ω) * kron(I(2), a' * a) +
          T +
          T' atol = 1e-12
    h, w = emission_rule(0, 0, Z_AXIS)
    L = ustrip.(u"µs^(-1/2)", t.model.jump_operators[1])
    η0 = t.mc.recoil_lamb_dicke[1]
    @test η0 ≈ abs(η)
    for j in (1, 7)
        @test Matrix(exact.c_ops[j].data) ≈
              sqrt(w[j]) * kron(L, displacement_elements(-η0 * h[j], num_fock)) atol =
            1e-12
    end
    # Σ_j L_j†L_j = L†L ⊗ 1 (the displacements are unitary), so the total decay
    # rate is untouched — up to the leakage out of the truncated ladder, which
    # only the top Fock states see (η₀√n ≈ 0.3 there); the recoil heats by
    # α η₀² per emission at every order, as in the first-order model.
    G = sum(Matrix(c.data)' * Matrix(c.data) for c in exact.c_ops)
    low = [i * num_fock + n for i in 0:1 for n in 1:6]
    @test G[low, low] ≈ kron(L' * L, I(num_fock))[low, low] atol = 1e-8
    @test real(G[2num_fock, 2num_fock]) < 0.99 * real((L'*L)[2, 2])
    N_op = kron(I(2), a' * a)
    ρ_e0 = kron(
        [0 0; 0 1.0],
        [n == m == 1 ? 1.0 : 0.0 for n in 1:num_fock, m in 1:num_fock],
    )
    function heating(mm)
        ops = [Matrix(c.data) for c in mm.c_ops]
        dissipated =
            sum(c * ρ_e0 * c' - (c' * c * ρ_e0 + ρ_e0 * c' * c) / 2 for c in ops)
        real(tr(N_op * dissipated))
    end
    @test heating(exact) ≈ η0^2 / 3 * ustrip(u"µs^-1", γ) rtol = 1e-6
    @test heating(first_order) ≈ η0^2 / 3 * ustrip(u"µs^-1", γ) rtol = 1e-6
    # A finite order k ≥ 2 gets it only to O(η^{k+2}): the truncated series
    # leaves spurious terms (here +0.6 %), which the exact operator cancels.
    @test heating(second) ≈ η0^2 / 3 * ustrip(u"µs^-1", γ) rtol = 1e-2
    @test !isapprox(heating(second), η0^2 / 3 * ustrip(u"µs^-1", γ); rtol=1e-3)

    # The orders converge: first order agrees with the exact model at O(η²),
    # second order at O(η³), in the Liouvillian.
    for η in (1e-3, 2e-3)
        s = sideband_toy(; γ, ω, η, Ω, δ=(-ω))
        ex = motional_model(s.model, s.mc; num_fock, lamb_dicke_order=Inf)
        ld = motional_model(s.model, s.mc; num_fock)
        sq = motional_model(s.model, s.mc; num_fock, lamb_dicke_order=2)
        scale = ustrip(u"µs^-1", Ω) * num_fock
        @test maximum(abs, ex.H.data - ld.H.data) < 2 * η^2 * scale
        @test maximum(abs, ex.H.data - sq.H.data) < 2 * η^3 * scale
        Lex = liouvillian(ex).data
        @test maximum(abs, Lex - liouvillian(ld).data) <
              5 * η^2 * ustrip(u"µs^-1", γ) * num_fock^2
        @test maximum(abs, Lex - liouvillian(sq).data) <
              5 * η^3 * ustrip(u"µs^-1", γ) * num_fock^2
        @test maximum(abs, Lex - liouvillian(ld).data) > η^2 * ustrip(u"µs^-1", γ) / 10
    end
    # A constant recoil moment does not define the emission distribution.
    c = sideband_toy(; γ, ω, η=0.1, Ω, δ=(-ω), recoil_moment=2 / 5)
    @test c.mc.emission == :constant
    @test motional_model(c.model, c.mc; num_fock=4).c_labels[2] isa RecoilLabel
    @test_throws ArgumentError motional_model(
        c.model,
        c.mc;
        num_fock=4,
        lamb_dicke_order=2,
    )
    @test_throws ArgumentError motional_model(
        c.model,
        c.mc;
        num_fock=4,
        lamb_dicke_order=Inf,
    )

    # Band Liouvillian from the model: the full band is the QuantumToolbox
    # Liouvillian up to the index permutation.
    bl = BandLiouvillian(exact)
    Lf = liouvillian(exact).data
    perm = [(c - 1) * 2num_fock + r for (r, c) in bl.entries]
    @test maximum(abs, bl.L - Lf[perm, perm]) < 1e-12
    @test BandLiouvillian(exact; bandwidth=3).bandwidth == 3

    # Heating bath: two operators per mode, raising the phonon number at the
    # bath rate for any state.
    rate = 2.0u"ms^-1"
    heated = motional_model(t.model, t.mc; num_fock, heating_rate=rate)
    @test count(l -> l isa HeatingLabel, heated.c_labels) == 2
    @test heated.c_labels[end] == HeatingLabel(1, true)
    ρ_th = thermal_state(heated, steadystate(t.model), 1.5)
    @test tr(ρ_th.data) ≈ 1
    bath = [c for (c, l) in zip(heated.c_ops, heated.c_labels) if l isa HeatingLabel]
    dN = sum(
        Matrix(c.data) * ρ_th.data * Matrix(c.data)' -
        (
            Matrix(c.data)' * Matrix(c.data) * ρ_th.data +
            ρ_th.data * Matrix(c.data)' * Matrix(c.data)
        ) / 2 for c in bath
    )
    # (the a† operator loses the kick out of the top Fock state, so compare below the truncation)
    @test real(tr(N_op * dN)) ≈
          ustrip(u"µs^-1", rate) *
          (1 - num_fock * thermal_populations(1.5, num_fock)[end]) rtol = 1e-9
    @test_throws ArgumentError motional_model(
        t.model,
        t.mc;
        num_fock,
        heating_rate=-1.0u"s^-1",
    )
    @test_throws ArgumentError motional_model(
        t.model,
        t.mc;
        num_fock,
        heating_rate=1.0u"m",
    )

    # Thermal states, Fock populations and phonon numbers on QuantumToolbox
    # objects; the two propagators (mesolve on the full model, Verner on the
    # band vector) agree.
    @test mean_phonon_number(ρ_th, first_order) ≈
          sum((0:(num_fock-1)) .* thermal_populations(1.5, num_fock))
    @test fock_populations(ρ_th, first_order) ==
          [fock_populations(ρ_th, first_order)[1]]
    @test fock_populations(ρ_th, first_order)[1] ≈ thermal_populations(1.5, num_fock)
    @test fock_populations(bl, band_vector(bl, Matrix(ρ_th.data)))[1] ≈
          thermal_populations(1.5, num_fock)
    ts = range(0.0, 40.0; length=21)
    full = cooling_curve(exact, ρ_th, ts)
    band = cooling_curve(
        BandLiouvillian(exact; bandwidth=4),
        Matrix(ρ_th.data),
        ts;
        fock=true,
    )
    @test full.t == collect(ts) && size(full.nbar) == (21, 1) && isnothing(full.fock)
    @test band.nbar ≈ full.nbar rtol = 1e-4
    @test band.fock[1][:, 1] ≈ thermal_populations(1.5, num_fock) atol = 1e-8
    unitful = cooling_curve(exact, ρ_th, ts .* u"µs"; fock=true)
    @test unitful.nbar ≈ full.nbar
    @test unitful.fock[1] ≈ band.fock[1] rtol = 1e-4
    @test steadystate(exact) isa typeof(ρ_th)
end

@testitem "Cooling metrics: the three methods agree in the rate-equation regime" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: steadystate

    # Weak drive, small η: the full model relaxes like the adiabatic-elimination
    # rate equation, whatever the Lamb–Dicke order.
    ω = 2π * 0.6u"MHz"
    γ = 2π * 50.0u"kHz"
    s = sideband_toy(; γ, ω, η=0.02, Ω=2π * 20.0u"kHz", δ=(-ω))
    ae = cooling_metrics(s.model, s.mc, AdiabaticElimination())
    @test ae isa CoolingMetrics
    @test ae.method == AdiabaticElimination() &&
          isnothing(ae.model) &&
          isnothing(ae.solver)
    rates = cooling_rates(s.model, s.mc)
    @test ae.nbar == rates.nbar && ae.τ_c == rates.τ_c
    @test ae.A_plus == rates.A_plus && ae.A_minus == rates.A_minus
    # the low-intensity limit 7/48 (γ/ω)² of isotropic emission
    @test ae.nbar ≈ 7 / 48 * ustrip(NoUnits, γ / ω)^2 rtol = 0.02
    it =
        cooling_metrics(s.model, s.mc, IntegratedTransient(; nbar_ini=1.0, num_fock=20))
    @test it.method.num_fock == 20 && isinf(it.method.lamb_dicke_order)
    @test it.model isa MotionalModel && it.solver isa IntegratedTransientSolver
    @test isnothing(it.A_plus)
    @test it.nbar ≈ ae.nbar rtol = 0.02
    @test it.τ_c ≈ ae.τ_c rtol = 0.02
    it1 = cooling_metrics(
        s.model,
        s.mc,
        IntegratedTransient(;
            nbar_ini=1.0,
            num_fock=20,
            lamb_dicke_order=1,
            bandwidth=nothing,
        ),
    )
    @test it1.model.lamb_dicke_order == 1
    @test it1.τ_c ≈ ae.τ_c rtol = 0.02
    sm = cooling_metrics(s.model, s.mc, LiouvillianSpectrum(; num_fock=8))
    @test sm.model isa MotionalModel && isnothing(sm.solver)
    @test sm.nbar ≈ ae.nbar rtol = 0.02
    @test sm.τ_c ≈ ae.τ_c rtol = 0.03
    # Other initial states from the solver carried in the result.
    ρ0 = thermal_state(Matrix(steadystate(s.model).data), it.model.num_fock, 3.0)
    @test integral_relaxation_time(it.solver, ρ0) * u"µs" ≈ ae.τ_c rtol = 0.02
    # The default truncation follows the thermal tail.
    auto = IntegratedTransient(; nbar_ini=2.0)
    @test isnothing(auto.num_fock) && auto.bandwidth == 4 && isnothing(auto.nodes)
    @test cooling_metrics(s.model, s.mc, auto).model.num_fock == [fock_truncation(2.0)]
    @test_throws ArgumentError IntegratedTransient(; nbar_ini=-1.0)
    @test_throws ArgumentError IntegratedTransient(; nbar_ini=1.0, lamb_dicke_order=0)
    @test_throws ArgumentError LiouvillianSpectrum(; num_fock=8, lamb_dicke_order=2.5)
    # A heating bath shifts the occupation alike in the rate equations and the
    # full model.
    h = 1.0u"s^-1"
    ae_h = cooling_metrics(s.model, s.mc, AdiabaticElimination(); heating_rate=h)
    it_h = cooling_metrics(
        s.model,
        s.mc,
        IntegratedTransient(; nbar_ini=1.0, num_fock=20);
        heating_rate=h,
    )
    @test ae_h.nbar > 1.5ae.nbar
    @test it_h.nbar ≈ ae_h.nbar rtol = 0.03
    @test count(l -> l isa HeatingLabel, it_h.model.c_labels) == 2
    # Several modes: one joint model, per-mode metrics (the second, far
    # off-resonant mode barely couples, so its full-model occupation is set by
    # its tiny rates and compares loosely).
    modes = [s.mode, MotionalMode(2π * 2.5u"MHz", Z_AXIS)]
    mcs =
        motional_coupling(s.species, s.scheme, s.model, modes; recoil_moment=:isotropic)
    many_ae = cooling_metrics(s.model, mcs, AdiabaticElimination())
    many_it = cooling_metrics(
        s.model,
        mcs,
        IntegratedTransient(; nbar_ini=[1.0, 0.5], num_fock=[12, 4]),
    )
    @test length(many_ae) == 2 && length(many_it) == 2
    @test many_ae[1].nbar == ae.nbar
    @test many_it[1].model === many_it[2].model &&
          many_it[1].solver === many_it[2].solver
    @test many_it[1].model.num_fock == [12, 4]
    @test many_it[1].nbar ≈ ae.nbar rtol = 0.03
    @test many_it[1].τ_c ≈ ae.τ_c rtol = 0.03
    @test_throws ArgumentError cooling_metrics(
        s.model,
        mcs,
        IntegratedTransient(; nbar_ini=[1.0], num_fock=4),
    )
    # The internal problem alone: the integrated transient of a level population
    # — the quench-rate construction — on a LindbladModel.
    decaying =
        sideband_toy(; γ=2π * 1.0u"MHz", ω, η=0.05, Ω=1e-9u"µs^-1", δ=0.0u"µs^-1")
    solver = IntegratedTransientSolver(
        decaying.model,
        level_projector(decaying.model, "P_3/2"),
    )
    @test integral_relaxation_time(solver, [0 0; 0 1.0]) ≈
          1 / ustrip(u"µs^-1", 2π * 1.0u"MHz") rtol = 1e-6
    @test steady_value(solver) ≈ 0 atol = 1e-9
end

@testitem "Integrated transient on a long Fock ladder" tags=[:integration, :slow] setup=[
    OpticalBlochSetup,
] begin
    using QuantumToolbox: steadystate

    # A Doppler-cooled ion at η = 0.1: the steady state spans > 100 orders of
    # magnitude along the ladder, where UMFPACK's default threshold pivoting
    # silently returned τ ≈ −7e9; strict pivoting gives the converged value
    # (Kulosa et al. 2023, Fig. 1(b) scenario; regression from the
    # quenched-sideband-cooling example).
    ω = 2π * 0.6u"MHz"
    γ = 2π * 50.0u"kHz"
    s = sideband_toy(; γ, ω, η=0.1, Ω=2π * 500.0u"kHz", δ=(-ω))
    num_fock = fock_truncation(20.0)
    @test num_fock == 189
    m = cooling_metrics(s.model, s.mc, IntegratedTransient(; nbar_ini=20.0))
    @test m.model.num_fock == [num_fock]
    @test m.τ_c * γ ≈ 81.5 rtol = 2e-3
    @test 0.05 < m.nbar < 0.2
    # Relaxation from a colder start is faster, and the integral time is a
    # linear function of the initial occupation (one phonon per ≈ 2 lifetimes).
    ρ_int = Matrix(steadystate(s.model).data)
    τ5 = integral_relaxation_time(m.solver, thermal_state(ρ_int, num_fock, 5.0))
    τ10 = integral_relaxation_time(m.solver, thermal_state(ρ_int, num_fock, 10.0))
    @test τ5 < τ10 < ustrip(u"µs", m.τ_c)
    @test (τ10 - τ5) * ustrip(u"µs^-1", γ) / 5 ≈ 2 rtol = 0.25
end
