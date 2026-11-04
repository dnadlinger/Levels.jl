@testitem "Level lifetimes" tags = [:unit, :fast] begin
    using Unitful

    # 88Sr+ P_3/2 decays to S_1/2, D_3/2 and D_5/2 (the last rate being the
    # 8.06 µs⁻¹ implied by the [Jiang2009] matrix element, cf. species_data).
    @test lifetime(sr88, "P_3/2") ≈ 1 / ((141 + 1.0 + 8.056)u"µs^-1") rtol = 1e-4

    # The ground state does not decay.
    @test isnothing(lifetime(sr88, "S_1/2"))
end

@testitem "ca40 data" tags = [:unit, :fast] begin
    using Unitful

    # The rate, g-factor and polarisability tables are shared between the Ca⁺
    # isotopes (all ⁴⁰Ca⁺ measurements); only masses and level energies differ.
    for (lo, hi) in keys(ca43.einstein_as)
        @test einstein_a(ca40, lo => hi) == einstein_a(ca43, lo => hi)
    end
    @test ca40.lande_g_overrides == ca43.lande_g_overrides
    @test ca40.mass < ca43.mass
    @test ca40.mass ≈ 39.9626u"u" rtol = 1e-4

    # Lifetimes reproduce the entered measurements: τ(P_1/2) [Hettrich2015],
    # τ(P_3/2) [Meir2020], τ(D_5/2) [Shao2017] and τ(D_3/2) = 1.0257 τ(D_5/2)
    # [Shao2018] (consistent with the direct 1.195(8) s of [Shao2016]).
    @test lifetime(ca40, "P_1/2") ≈ 6.904u"ns" rtol = 1e-6
    # (Branching fractions sum to 1.00001, hence the tolerance.)
    @test lifetime(ca40, "P_3/2") ≈ 6.639u"ns" rtol = 2e-5
    @test lifetime(ca40, "D_5/2") ≈ 1.1649u"s" rtol = 1e-6
    @test lifetime(ca40, "D_3/2") ≈ 1.195u"s" rtol = 1e-3

    # Measured g-factors for S_1/2 and D_5/2, LS coupling elsewhere.
    @test lande_g(ca40, "S_1/2") == 2.00225664
    @test lande_g(ca40, "D_5/2") ≈ 1.20033046 atol = 1e-8
    @test lande_g(ca40, "P_1/2") ≈ 2 // 3 rtol = 2e-3
    @test lande_g(ca40, "D_3/2") ≈ 4 // 5 rtol = 2e-3

    # Ring closure of the level energies (729 nm [Zhang2023] + fine-structure
    # interval [Solaro2018], 397 nm [Gebert2015], 393 nm [Shi2017]) against the
    # directly measured 3d → 4p intervals of [Muller2020], Table III.
    interval(lo, hi) =
        uconvert(u"MHz", cyclic(Levels.transition_frequency(ca40, lo => hi)))
    @test interval("D_3/2", "P_1/2") ≈ 346_000_235.13u"MHz" atol = 0.2u"MHz"
    @test interval("D_3/2", "P_3/2") ≈ 352_682_481.93u"MHz" atol = 0.2u"MHz"
    @test interval("D_5/2", "P_3/2") ≈ 350_862_882.63u"MHz" atol = 0.3u"MHz"
    @test interval("S_1/2", "D_5/2") ≈ 411_042_129.77640026u"MHz" rtol = 1e-12
end

@testitem "ba138 data" tags = [:unit, :fast] begin
    using Unitful
    using WignerSymbols: wigner6j

    cyclic(pair) = uconvert(u"MHz", Levels.transition_frequency(ba138, pair) / 2π)

    @test ba138.mass ≈ 137.9047u"u" rtol = 1e-6

    # Level energies (clock [Arnold2020] + D fine structure [Whitford1994], 493 and
    # 455 nm [Imgram2019]) against the independent direct measurements of [Dijck2015]
    # and the line centre of [Jayjong2026b].
    @test cyclic("D_3/2" => "P_1/2") ≈ 461_311_878.5u"MHz" atol = 0.2u"MHz"
    @test cyclic("S_1/2" => "D_3/2") ≈ 146_114_384.0u"MHz" atol = 0.2u"MHz"
    @test cyclic("S_1/2" => "P_1/2") ≈ 607_426_262.6455u"MHz" atol = 0.13u"MHz"

    # The dipoles of [Jayjong2026b] reproduce the lifetimes it states (Eqs. 16–17)…
    @test lifetime(ba138, "P_1/2") ≈ 7.8663u"ns" rtol = 7e-4
    @test lifetime(ba138, "P_3/2") ≈ 6.2900u"ns" rtol = 7.5e-4
    # … and the independently measured branching of [Zhang2020] (ratios p and q).
    a(pair) = einstein_a(ba138, pair)
    s, d3, d5 = a("S_1/2" => "P_3/2"), a("D_3/2" => "P_3/2"), a("D_5/2" => "P_3/2")
    @test s / (s + d5) ≈ 0.763107 atol = 3 * 65e-6
    @test d3 / (d3 + d5) ≈ 0.108527 atol = 3 * 77e-6
    # P_1/2 → D_3/2 fraction, as adopted by [Jayjong2026b].
    d1, s1 = a("D_3/2" => "P_1/2"), a("S_1/2" => "P_1/2")
    @test d1 / (d1 + s1) ≈ 0.268167 atol = 3 * 47e-6

    # D_5/2 decays by E2 to S_1/2 and by M1 to D_3/2, at the measured total rate.
    @test length(keys(ba138.einstein_as)) == 8
    @test lifetime(ba138, "D_5/2") ≈ 30.14u"s" rtol = 1e-12
    @test Levels.multipole(ba138, "D_3/2" => "D_5/2") == :M1
    @test Levels.multipole(ba138, "D_5/2" => "D_3/2") == :M1
    @test Levels.multipole(ba138, "S_1/2" => "D_5/2") == :E2
    @test Levels.multipole(ba138, "S_1/2" => "P_3/2") == :E1
    # The M1 rate against its LS-coupling value: between fine-structure partners,
    # μ = −μ_B (L + g_s S) reduces to −μ_B (g_s − 1) S, whose reduced element follows
    # from recoupling l and s (g_l = 1 up to the reduced-mass correction).
    l, s, j_lo, j_hi = 2, 1//2, 3//2, 5//2
    s_reduced =
        sqrt((2j_lo + 1) * (2j_hi + 1)) *
        wigner6j(s, j_hi, l, j_lo, s, 1) *
        sqrt(s * (s + 1) * (2s + 1))
    @test s_reduced^2 ≈ 12 / 5
    μ = (2.00231930436 - 1) * abs(s_reduced) * u"μB"
    ω = Levels.transition_frequency(ba138, "D_3/2" => "D_5/2")
    a_ls = Levels.einstein_a_from_magnetic_dipole(μ, ω, j_hi)
    @test einstein_a(ba138, "D_3/2" => "D_5/2") ≈ a_ls rtol = 0.02
    @test lifetime(ba138, "D_3/2") ≈ 82.4u"s" rtol = 1e-12

    # Measured g-factors, including the field-independent ratio of [Knoll1996].
    @test lande_g(ba138, "S_1/2") == 2.00249192
    @test lande_g(ba138, "D_5/2") == 1.20036739
    @test lande_g(ba138, "S_1/2") / lande_g(ba138, "D_3/2") ≈ 2.505220 atol = 3e-6
end
