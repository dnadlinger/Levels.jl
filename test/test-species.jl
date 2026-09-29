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
        @test einstein_a(ca40, lo, hi) == einstein_a(ca43, lo, hi)
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
    interval(lo, hi) = uconvert(u"MHz", Levels.transition_frequency(ca40, lo, hi) / 2π)
    @test interval("D_3/2", "P_1/2") ≈ 346_000_235.13u"MHz" atol = 0.2u"MHz"
    @test interval("D_3/2", "P_3/2") ≈ 352_682_481.93u"MHz" atol = 0.2u"MHz"
    @test interval("D_5/2", "P_3/2") ≈ 350_862_882.63u"MHz" atol = 0.3u"MHz"
    @test interval("S_1/2", "D_5/2") ≈ 411_042_129.77640026u"MHz" rtol = 1e-12
end
