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

@testitem "Physical wrappers retain time units" tags=[:unit, :fast] begin
    using Levels.OpticalBloch
    using Unitful, SparseArrays, LinearAlgebra
    num_fock = 4
    a = spdiagm(1 => sqrt.(1.0:(num_fock-1)))
    function damping(rate, time_unit)
        BandLiouvillian(
            spzeros(num_fock, num_fock),
            [sqrt(rate)*a],
            1,
            num_fock;
            time_unit,
        )
    end
    microseconds = damping(0.3, u"µs")
    seconds = damping(3e5, u"s")
    observable = phonon_number_operator(1, num_fock)
    ρ0 = thermal_state(ones(1, 1), num_fock, 0.5)
    first = IntegratedTransientSolver(microseconds, observable)
    second = IntegratedTransientSolver(seconds, observable)
    @test integral_relaxation_time(first, ρ0) ≈ integral_relaxation_time(second, ρ0)
    @test integral_relaxation_time(first, ρ0) ≈ (1/0.3)u"µs"
    times = [0.0, 1.0, 2.0]u"µs"
    left = cooling_curve(microseconds, ρ0, times)
    right = cooling_curve(seconds, ρ0, times)
    @test left.t ≈ right.t
    @test left.nbar ≈ right.nbar
end

@testitem "Omitted decay states require an explicit policy" tags=[:unit, :fast] begin
    using Levels.OpticalBloch
    using Unitful
    basis = StateBasis(
        vcat([StateSpec("S_1/2", -1//2)], collect(StateBasis("D_3/2", "P_1/2"))),
    )
    scheme = LaserScheme(basis, LaserBeam[]; B=0.5u"mT")
    @test_throws ArgumentError lindblad_model(ca40, scheme)
    dropped = lindblad_model(ca40, scheme; open_channels=:drop)
    @test length(dropped.jump_operators)>0
end

@testitem "Relative couplings are normalised over the geometry" tags=[:unit, :fast] begin
    using Unitful
    direction, ε = beam_vectors(π / 2, π / 4)
    basis = StateBasis("S_1/2", "P_1/2")
    E1 = relative_couplings(ca40, basis, "S_1/2" => "P_1/2", ε, direction)
    @test any(!iszero, E1)
    @test E1 ≈ relative_couplings(ca40, basis, "S_1/2" => "P_1/2", 2ε, 3direction)
end
