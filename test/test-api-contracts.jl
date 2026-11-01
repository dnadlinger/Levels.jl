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
    @test angular(1u"MHz") ≈ 2π * 1u"rad/µs"
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

@testitem "Angular frequencies carry rad, rates do not" tags=[:unit, :fast] begin
    using Levels.PeriodicDynamics
    using Levels.OpticalBloch
    using Unitful
    rad_per_µs = u"rad/µs"
    B = 0.5u"mT"
    direction, ε = beam_vectors(π / 2, π / 4)
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    @test unit(angular(1u"MHz")) == rad_per_µs
    @test unit(cyclic(angular(1u"MHz"))) == u"Hz"
    @test unit(transition_frequency(sr88, "S_1/2" => "D_5/2")) == rad_per_µs
    @test unit(transition_frequency(ca43, "S_1/2 F=4" => "D_5/2 F=4")) == rad_per_µs
    @test unit(zeeman_shift(sr88, probe.first; B)) == rad_per_µs
    @test unit(zeeman_sensitivity(sr88, probe)) == u"rad/µs/mT"
    @test unit(Levels.hyperfine_shift(ca43, "S_1/2 F=4")) == rad_per_µs
    manifold = hyperfine_manifold(ca43, "S_1/2"; B)
    @test unit(state_offset(manifold, manifold.basis[1])) == rad_per_µs
    @test unit(state_energy(manifold, manifold.basis[1])) == u"J"
    hf_pair = StateSpec("S_1/2 F=4", 0) => StateSpec("S_1/2 F=3", 0)
    @test unit(zeeman_sensitivity(ca43, hf_pair; B)) == u"rad/µs/mT"
    @test unit(rabi_frequency(sr88, probe, 1.0u"W/m^2", ε, direction)) == rad_per_µs
    basis = StateBasis(sr88, "S_1/2", "D_5/2")
    C = coupling_matrix(sr88, basis, "S_1/2" => "D_5/2", 1.0u"W/m^2", ε, direction)
    @test unit(eltype(C)) == rad_per_µs
    @test unit(RelativeFrequency(sr88, probe, angular(1u"kHz"); B).offset) == rad_per_µs
    c = LightShiftCoefficients(sr88, basis, 674.0u"nm"; B)
    @test unit(driven_light_shift(c, probe, 1.0u"W/m^2", ε; direction)) == rad_per_µs
    @test unit(light_shift(c, probe.first, 1.0u"W/m^2", ε)) == rad_per_µs
    parked = LightShiftCoefficients(
        sr88,
        basis,
        RelativeFrequency("S_1/2" => "D_5/2", angular(1u"kHz"));
        B,
    )
    @test unit(light_shift(parked, probe.first, 1.0u"W/m^2", ε; direction)) ==
          rad_per_µs
    @test unit(
        light_shift(parked, probe.first, 1.0u"W/m^2", ε; direction, parts=:resonant),
    ) == rad_per_µs
    @test unit(polarisability(sr88, probe.first, 0u"rad/s", ε)) == u"C*m^2/V"
    # Rates are plain inverse times.
    @test unit(einstein_a(sr88, "S_1/2" => "P_1/2")) == u"µs^-1"
    @test unit(lifetime(sr88, "P_1/2")) == u"s"
    @test unit(saturation_intensity(sr88, "S_1/2" => "P_1/2")) == u"W/m^2"
    # Periodic dynamics.
    dt = DrivenTransition(
        sr88,
        basis,
        probe;
        B,
        drive_frequency=angular(50u"MHz"),
        coupling=rabi_normalised(C, basis, probe, angular(100u"kHz")),
        drives=zeeman_drives(sr88, basis, [1.0, 0.0, 0.0]u"µT"),
    )
    @test unit(eltype(dt.frame)) == rad_per_µs
    @test unit(eltype(dt.coupling)) == rad_per_µs
    sideband = sideband_rabi(dt; sideband=1)
    @test unit(sideband.Ω) == rad_per_µs && unit(sideband.δ_res) == rad_per_µs
    shift = ac_zeeman_shift(sr88, probe, angular(50u"MHz"), [1.0, 0.0, 0.0]u"µT"; B)
    @test unit(shift.shift) == rad_per_µs
    # Optical Bloch: angular Hamiltonian, √rate jump operators.
    beam = LaserBeam("S_1/2" => "P_1/2", angular(-10u"MHz"), 1.0u"W/m^2", ε, direction)
    scheme = LaserScheme(StateBasis(sr88, "S_1/2", "P_1/2"), [beam]; B)
    model = lindblad_model(sr88, scheme; open_channels=:drop)
    @test unit(eltype(model.hamiltonian)) == rad_per_µs
    @test unit(eltype(first(model.jump_operators))) == u"µs^(-1/2)"
end
