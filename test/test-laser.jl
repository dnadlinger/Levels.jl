@testitem "RelativeFrequency" tags=[:unit, :fast] begin
    using Unitful

    # The reference levels must be a distinct, same-kind pair, and the offset
    # an angular frequency.
    @test_throws ArgumentError RelativeFrequency("S_1/2" => "S_1/2", 0.0u"s^-1")
    @test_throws ArgumentError RelativeFrequency("S_1/2" => "D_5/2 F=4", 0.0u"s^-1")
    @test_throws ArgumentError RelativeFrequency("S_1/2" => "D_5/2", 5.0u"nm")

    # Resolving to an absolute frequency needs the atomic data — and inherits
    # its precision, which is the point of keeping the offset separate.
    rel = RelativeFrequency("S_1/2" => "D_5/2", 2π * 1.0u"MHz")
    @test_throws ArgumentError Levels.photon_energy(rel)
    ω0 = Levels.transition_frequency(sr88, "S_1/2", "D_5/2")
    @test Levels.photon_energy(sr88, rel) ≈ uconvert(u"J", u"ħ" * (ω0 + 2π * 1.0u"MHz"))

    # The wavelength and angular-frequency forms agree with each other (and
    # ignore a species argument).
    λ = 674.0u"nm"
    @test Levels.photon_energy(λ) ≈
          Levels.photon_energy(uconvert(u"s^-1", 2π * u"c" / λ))
    @test Levels.photon_energy(sr88, λ) == Levels.photon_energy(λ)
end
