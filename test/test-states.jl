@testitem "Spectroscopic notation parser" tags=[:unit, :fast] begin
    is_lj(val, l, j) = NoHyperfineNumberSpec(val) == NoHyperfineNumberSpec(l, j)
    @test is_lj("S_{1/2}", 0, 1//2)
    @test is_lj("P_{3/2}", 1, 3//2)
    @test is_lj("D_{5/2}", 2, 5//2)

    @test is_lj("D_{5//2}", 2, 5//2)
    @test is_lj("D_5/2", 2, 5//2)
    @test is_lj("D_5//2", 2, 5//2)
    @test is_lj("D5/2", 2, 5//2)
    @test is_lj("D5//2", 2, 5//2)

    # Not just `String`s; any `AbstractString` will do.
    @test is_lj(split("D_5/2 (metastable)")[1], 2, 5//2)

    @test_throws ArgumentError NoHyperfineNumberSpec("X_1/2")
    @test_throws ArgumentError NoHyperfineNumberSpec("")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_5//2a")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_52")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_a")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_1/2}")
    @test_throws ArgumentError NoHyperfineNumberSpec("D_{1/2")
end

@testitem "Hyperfine level parsing" tags=[:unit, :fast] begin
    using Levels: momentum, parse_level

    is_ljf(val, l, j, f) = HyperfineNumberSpec(val) == HyperfineNumberSpec(l, j, f)
    @test is_ljf("S_1/2 F=4", 0, 1//2, 4//1)
    @test is_ljf("S_{1/2} F=4", 0, 1//2, 4//1)
    @test is_ljf("S_1/2 F = 4", 0, 1//2, 4//1)
    @test is_ljf("D_5/2 F=7/2", 2, 5//2, 7//2)
    @test is_ljf("D5//2 F=7//2", 2, 5//2, 7//2)

    # The 4f fine-structure levels remain parseable; the suffix needs whitespace.
    @test NoHyperfineNumberSpec("F_7/2") == NoHyperfineNumberSpec(3, 7//2)

    # Kind routing: hyperfine iff the F suffix is present.
    @test parse_level("S_1/2") == NoHyperfineNumberSpec(0, 1//2)
    @test parse_level("S_1/2 F=4") == HyperfineNumberSpec(0, 1//2, 4//1)
    @test parse_level(HyperfineNumberSpec(0, 1//2, 4)) ==
          HyperfineNumberSpec(0, 1//2, 4)

    @test fine_structure(HyperfineNumberSpec(2, 5//2, 3)) ==
          NoHyperfineNumberSpec(2, 5//2)
    @test fine_structure("D_5/2 F=3") == NoHyperfineNumberSpec(2, 5//2)
    @test fine_structure("D_5/2") == NoHyperfineNumberSpec(2, 5//2)
    @test momentum(NoHyperfineNumberSpec(2, 5//2)) == 5//2
    @test momentum("D_5/2 F=3") == 3//1

    # A fine-structure parse must not silently drop the F suffix, and vice versa.
    @test_throws ArgumentError NoHyperfineNumberSpec("S_1/2 F=4")
    @test_throws ArgumentError HyperfineNumberSpec("S_1/2")
    @test_throws ArgumentError HyperfineNumberSpec("S_1/2 F=a")
    @test_throws ArgumentError HyperfineNumberSpec("S_1/2 F=4/")
    @test_throws ArgumentError HyperfineNumberSpec("X_1/2 F=4")

    @test StateSpec("S_1/2 F=4", 2) == StateSpec(HyperfineNumberSpec(0, 1//2, 4), 2//1)
end

@testitem "StateSpec construction" tags=[:unit, :fast] begin
    @test StateSpec("D_5/2", 3//2) == StateSpec(NoHyperfineNumberSpec(2, 5//2), 3//2)
    @test StateSpec("P_1/2", 1//2).m == 1//2

    # Integer projections convert to the rational field type.
    @test StateSpec("P_1/2", 0).m == 0//1

    # Spectroscopic level parts canonicalise via convert.
    spectro = StateSpec(SpectroscopicSpec("D_5/2"), 3//2)
    @test convert(StateSpec{NoHyperfineNumberSpec}, spectro) == StateSpec("D_5/2", 3//2)
end

@testitem "state_pairs" tags=[:unit, :fast] begin
    # The eight quadrupole components observable with Δm = 0 suppressed.
    e2 = state_pairs("S_1/2", "D_5/2"; Δm=[-2, -1, 1, 2])
    @test length(e2) == 8
    @test all(p -> abs(p.second.m - p.first.m) in 1:2, e2)
    @test (StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -5//2)) in e2
    @test !((StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -1//2)) in e2)

    @test length(state_pairs("S_1/2", "D_5/2"; Δm=-2:2)) == 10

    # E1-type enumeration, and a Δm = 0 (π-transition) singleton collection.
    @test length(state_pairs("S_1/2", "P_1/2"; Δm=-1:1)) == 4
    @test state_pairs("S_1/2", "P_1/2"; Δm=0) == [
        StateSpec("S_1/2", -1//2) => StateSpec("P_1/2", -1//2),
        StateSpec("S_1/2", 1//2) => StateSpec("P_1/2", 1//2),
    ]

    # Δm is deliberately without default.
    @test_throws UndefKeywordError state_pairs("S_1/2", "D_5/2")

    # Hyperfine levels enumerate m_F; kinds must not mix.
    hf = state_pairs("S_1/2 F=4", "D_5/2 F=3"; Δm=0)
    @test length(hf) == 7
    @test hf[1] == (StateSpec("S_1/2 F=4", -3) => StateSpec("D_5/2 F=3", -3))
    # All 9 lower m_F reach 5 upper states each within |m_F| ≤ 6.
    @test length(state_pairs("S_1/2 F=4", "D_5/2 F=6"; Δm=-2:2)) == 9 * 5
    @test_throws ArgumentError state_pairs("S_1/2", "D_5/2 F=3"; Δm=0)
end
