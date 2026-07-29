@testitem "Level lifetimes" tags = [:unit, :fast] begin
    using Unitful

    # 88Sr+ P_3/2 decays to S_1/2, D_3/2 and D_5/2.
    @test lifetime(sr88, "P_3/2") ≈ 1 / ((141 + 1.0 + 8.7)u"µs^-1")

    # The ground state does not decay.
    @test isnothing(lifetime(sr88, "S_1/2"))
end
