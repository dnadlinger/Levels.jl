@testitem "Species input and nuclear-moment conventions" tags=[:unit, :fast] begin
    using Unitful
    species = NoHyperfineOneElectronSpecies(;
        mass=40u"u",
        energies=Dict("S_1/2" => 0u"J", "P_1/2" => 3u"eV"),
        einstein_as=Dict(("S_1/2", "P_1/2") => ReducedDipole(3 * Levels.DIPOLE_AU)),
    )
    @test einstein_a(species, "S_1/2" => "P_1/2") > 0u"s^-1"
    @test ca43.nuclear_g * ca43.nuclear_spin ≈ -1.315350
    @test_throws ArgumentError HyperfineConstants(; a=1u"m")
    @test_throws ArgumentError StateSpec("S_1/2", 3//2)
    @test_throws ArgumentError StateSpec("S_1/2", 0)
    @test_throws ArgumentError NoHyperfineOneElectronSpecies(;
        mass=-1u"u",
        energies=species.energies,
        einstein_as=species.einstein_as,
    )
end

@testitem "Pair queries, angular conversion and physical energy" tags=[:unit, :fast] begin
    using Unitful
    pair = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -1//2)
    B = 0.5u"mT"
    @test transition_frequency(ca40, pair; B) ≈
          transition_frequency(ca40, "S_1/2" => "D_5/2") +
          zeeman_shift(ca40, pair.second; B) - zeeman_shift(ca40, pair.first; B)
    @test zeeman_sensitivity(ca40, pair; B) == zeeman_sensitivity(ca40, pair)
    @test angular(1u"MHz") ≈ 2π * 1u"µs^-1"
    @test cyclic(angular(1u"MHz")) ≈ 1u"MHz"
    @test level_energy(ca40, "S_1/2") == 0u"J"
    @test collect(StateBasis(ca40, "S_1/2")) == collect(StateBasis("S_1/2"))
    @test all(p -> p isa Pair, keys(ca40.einstein_as))
    m = hyperfine_manifold(ca43, "S_1/2"; B)
    @test state_energy(m, m.basis[1]) ≈ u"ħ" * state_offset(m, m.basis[1])
end
