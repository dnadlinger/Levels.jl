using Levels
using Base.Test

@testset "Spectroscopic notation" begin
    is_lj(val, l, j) = NoHyperfineNumberSpec(val) == NoHyperfineNumberSpec(l, j)
    @test is_lj("S_{1/2}", 0, 1//2)
    @test is_lj("P_{3/2}", 1, 3//2)
    @test is_lj("D_{5/2}", 2, 5//2)

    @test is_lj("D_{5//2}", 2, 5//2)
    @test is_lj("D_5/2", 2, 5//2)
    @test is_lj("D_5//2", 2, 5//2)
    @test is_lj("D5/2", 2, 5//2)
    @test is_lj("D5//2", 2, 5//2)

    @test_throws ArgumentError NoHyperfineNumberSpec("D_5//2a")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_52")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_a")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_1/2}")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_{1/2")
end
