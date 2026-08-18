@testitem "Level lifetimes" tags = [:unit, :fast] begin
    using Unitful

    # 88Sr+ P_3/2 decays to S_1/2, D_3/2 and D_5/2 (the last rate being the
    # 8.06 µs⁻¹ implied by the [Jiang2009] matrix element, cf. species_data).
    @test lifetime(sr88, "P_3/2") ≈ 1 / ((141 + 1.0 + 8.056)u"µs^-1") rtol = 1e-4

    # The ground state does not decay.
    @test isnothing(lifetime(sr88, "S_1/2"))
end
