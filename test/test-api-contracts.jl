@testitem "Species and state input validation" tags=[:unit, :fast] begin
    using Unitful
    species = NoHyperfineOneElectronSpecies(;
        mass=40u"u",
        energies=Dict("S_1/2" => 0u"J", "P_1/2" => 3u"eV"),
        einstein_as=Dict(("S_1/2", "P_1/2") => ReducedDipole(3 * Levels.DIPOLE_AU)),
    )
    @test einstein_a(species, "S_1/2", "P_1/2") > 0u"s^-1"
    @test_throws ArgumentError HyperfineConstants(; a=1u"m")
    @test_throws ArgumentError StateSpec("S_1/2", 3//2)
    @test_throws ArgumentError StateSpec("S_1/2", 0)
    @test_throws ArgumentError NoHyperfineOneElectronSpecies(;
        mass=-1u"u",
        energies=species.energies,
        einstein_as=species.einstein_as,
    )
end
