@testitem "StateBasis from levels" tags=[:unit, :fast] begin
    basis = StateBasis(["S_1/2", "D_5/2"])
    @test length(basis) == 8

    # Blocks in the given order, each ordered by increasing m.
    @test stateindex(basis, "S_1/2", -1//2) == 1
    @test stateindex(basis, "S_1/2", 1//2) == 2
    @test stateindex(basis, "D_5/2", -5//2) == 3
    @test stateindex(basis, "D_5/2", 5//2) == 8
    @test stateindex(basis, StateSpec("D_5/2", 3//2)) == 7
    @test staterange(basis, "S_1/2") == 1:2
    @test staterange(basis, "D_5/2") == 3:8

    # Collection interface.
    @test basis[3] == StateSpec("D_5/2", -5//2)
    @test basis[end] == StateSpec("D_5/2", 5//2)
    @test collect(basis) == [basis[i] for i in 1:8]
    @test eltype(basis) == StateSpec{NoHyperfineNumberSpec}
    @test StateSpec("D_5/2", 1//2) in basis
    @test !(StateSpec("D_3/2", 1//2) in basis)

    # Vararg convenience and level order.
    basis2 = StateBasis("D_5/2", "S_1/2")
    @test staterange(basis2, "D_5/2") == 1:6
    @test staterange(basis2, "S_1/2") == 7:8

    @test_throws ArgumentError StateBasis(["S_1/2", "S_{1/2}"])
    @test_throws ArgumentError stateindex(basis, "P_1/2", 1//2)
    @test_throws ArgumentError stateindex(basis, "S_1/2", 3//2)
    @test_throws ArgumentError staterange(basis, "P_1/2")
end

@testitem "Hyperfine StateBasis" tags=[:unit, :fast] begin
    basis = StateBasis(["S_1/2 F=3", "S_1/2 F=4"])
    @test length(basis) == 7 + 9
    @test eltype(basis) == StateSpec{HyperfineNumberSpec}
    @test stateindex(basis, "S_1/2 F=3", -3) == 1
    @test stateindex(basis, "S_1/2 F=4", 4) == 16
    @test staterange(basis, "S_1/2 F=3") == 1:7
    @test staterange(basis, "S_1/2 F=4") == 8:16

    # A fine-structure query returns the whole-manifold range.
    @test staterange(basis, "S_1/2") == 1:16
    @test staterange(basis, NoHyperfineNumberSpec(0, 1//2)) == 1:16
    @test_throws ArgumentError staterange(basis, "D_5/2")

    # m_F outside the level, and one of incompatible parity.
    @test_throws ArgumentError StateBasis([StateSpec("S_1/2 F=4", 5)])
    @test_throws ArgumentError StateBasis([StateSpec("S_1/2 F=4", 1//2)])

    # Level kinds must not mix, whether given as levels or states.
    @test_throws ArgumentError StateBasis(["S_1/2 F=4", "D_5/2"])
    @test_throws ArgumentError StateBasis([
        StateSpec("S_1/2 F=4", 0),
        StateSpec("D_5/2", 1//2),
    ])

    # Manifold queries require contiguous F blocks.
    scattered = StateBasis([
        StateSpec("S_1/2 F=3", 0),
        StateSpec("D_5/2 F=2", 0),
        StateSpec("S_1/2 F=4", 0),
    ])
    @test staterange(scattered, "D_5/2") == 2:2
    @test_throws ArgumentError staterange(scattered, "S_1/2")

    # Hyperfine queries on a no-hyperfine basis fail loudly rather than
    # silently reinterpreting m_F as m_J.
    nh = StateBasis(["S_1/2"])
    @test_throws ArgumentError staterange(nh, "S_1/2 F=4")
    @test_throws ArgumentError staterange(nh, HyperfineNumberSpec(0, 1//2, 4))
end

@testitem "StateBasis from explicit states" tags=[:unit, :fast] begin
    basis = StateBasis([
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", -1//2),
        StateSpec("S_1/2", -1//2),
    ])
    @test length(basis) == 3
    @test stateindex(basis, "D_5/2", -1//2) == 2

    # S states sit at indices 1 and 3, so no contiguous range exists for them.
    @test staterange(basis, "D_5/2") == 2:2
    @test_throws ArgumentError staterange(basis, "S_1/2")

    @test_throws ArgumentError StateBasis(StateSpec{NoHyperfineNumberSpec}[])
    @test_throws ArgumentError StateBasis([
        StateSpec("S_1/2", 1//2),
        StateSpec("S_1/2", 1//2),
    ])
    # m out of range, and integer m for a half-integer j level.
    @test_throws ArgumentError StateBasis([StateSpec("S_1/2", 3//2)])
    @test_throws ArgumentError StateBasis([StateSpec("D_5/2", 1)])
end
