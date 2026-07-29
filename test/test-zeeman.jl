@testitem "Landé g-factors" tags=[:unit, :fast] begin
    # LS-coupling values with the reduced-mass-corrected g_L for ⁸⁸Sr⁺, as used in
    # hoa2-common (free-electron g_s, M from AME2020).
    @test lande_g(sr88, "S_1/2") ≈ 2.0023193043618 rtol = 1e-10
    @test lande_g(sr88, "D_5/2") ≈ 1.2004588684586 rtol = 1e-10

    # Uncorrected Landé values (g_s = 2, g_l = 1) as coarse sanity checks.
    @test lande_g(sr88, "P_1/2") ≈ 2 / 3 rtol = 2e-3
    @test lande_g(sr88, "P_3/2") ≈ 4 / 3 rtol = 2e-3
    @test lande_g(sr88, "D_3/2") ≈ 4 / 5 rtol = 2e-3
end

@testitem "Zeeman shifts and sensitivities" tags=[:unit, :fast] begin
    using Unitful

    μ_B = Unitful.q * u"ħ" / (2 * Unitful.me)
    s_up = StateSpec("S_1/2", 1//2)
    d_up = StateSpec("D_5/2", 5//2)

    shift = zeeman_shift(sr88, s_up, 0.5u"mT")
    @test shift ≈ uconvert(u"MHz", lande_g(sr88, "S_1/2") / 2 * μ_B * 0.5u"mT" / u"ħ")
    @test zeeman_shift(sr88, StateSpec("S_1/2", -1//2), 0.5u"mT") ≈ -shift

    # χ of the strongest S ↔ D transition, in the (non-angular) units of thesis
    # fig. 5.12: (g_D·5/2 − g_S·1/2)·μ_B/h = 27.9923 MHz/mT.
    χ = zeeman_sensitivity(sr88, s_up, d_up)
    @test uconvert(u"MHz/mT", χ / 2π) ≈ 27.9923u"MHz/mT" rtol = 1e-5
    @test zeeman_sensitivity(sr88, s_up => d_up) == χ

    # The eight quadrupole components sorted by sensitivity: ±{11.2, 22.4, 28.0,
    # 39.2} MHz/mT, antisymmetric under m → −m.
    pairs = sort!(
        state_pairs("S_1/2", "D_5/2"; Δm=[-2, -1, 1, 2]);
        by=t -> zeeman_sensitivity(sr88, t),
    )
    χs = [zeeman_sensitivity(sr88, t) for t in pairs]
    @test first(pairs) == (StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", -3//2))
    @test last(pairs) == (StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", 3//2))
    @test uconvert(u"MHz/mT", last(χs) / 2π) ≈ 39.215u"MHz/mT" rtol = 1e-4
    @test all(χs[i] ≈ -χs[9-i] for i in 1:8)
end

@testitem "Zeeman Hamiltonian" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful
    using Levels: jx_matrix, jy_matrix, jz_matrix

    basis = StateBasis(["S_1/2", "D_5/2"])

    # Field along the quantisation axis: diagonal of first-order shifts.
    H = zeeman_hamiltonian(sr88, basis, [0.0, 0.0, 0.5]u"mT")
    @test all(iszero, H - Diagonal(diag(H)))
    for (i, state) in enumerate(basis)
        @test real(H[i, i]) ≈ zeeman_shift(sr88, state, 0.5u"mT")
    end

    # General field direction: Hermitian, block-diagonal in the levels, and each
    # block matches g μ_B/ħ B·J on the spin-j matrices.
    B = [0.3, -0.2, 0.7]u"mT"
    H = zeeman_hamiltonian(sr88, basis, B)
    @test H ≈ H'
    s_range = staterange(basis, "S_1/2")
    d_range = staterange(basis, "D_5/2")
    @test all(iszero, H[s_range, d_range])
    @test all(iszero, H[d_range, s_range])
    μ_B = Unitful.q * u"ħ" / (2 * Unitful.me)
    for (level, j, r) in (("S_1/2", 1//2, s_range), ("D_5/2", 5//2, d_range))
        bj = B[1] * jx_matrix(j) + B[2] * jy_matrix(j) + B[3] * Matrix(jz_matrix(j))
        expected = uconvert.(u"MHz", lande_g(sr88, level) * μ_B / u"ħ" * bj)
        @test H[r, r] ≈ expected
    end

    # Elements are constructed per-state, so scrambled/partial bases work too.
    scrambled = StateBasis([
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", -1//2),
        StateSpec("S_1/2", -1//2),
    ])
    Hs = zeeman_hamiltonian(sr88, scrambled, B)
    @test Hs ≈ Hs'
    @test Hs[1, 3] ≈ H[2, 1]   # ⟨S,+1/2| · |S,−1/2⟩
    @test Hs[3, 1] ≈ H[1, 2]
    i = stateindex(basis, "D_5/2", -1//2)
    @test Hs[2, 2] ≈ H[i, i]
end
