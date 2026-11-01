@testitem "Species input and nuclear-moment conventions" tags=[:unit, :fast] begin
    using Unitful
    toy = NoHyperfineOneElectronSpecies(;
        mass=40u"u",
        energies=Dict("S_1/2" => 0u"J", "P_1/2" => 3u"eV"),
        einstein_as=Dict(("S_1/2", "P_1/2") => ReducedDipole(3 * Levels.DIPOLE_AU)),
    )
    @test einstein_a(toy, "S_1/2" => "P_1/2") > 0u"s^-1"
    @test ca43.nuclear_g * ca43.nuclear_spin ≈ -1.315350
    @test_throws ArgumentError HyperfineConstants(; a=1u"m")
    @test_throws ArgumentError StateSpec("S_1/2", 3//2)
    @test_throws ArgumentError StateSpec("S_1/2", 0)
    @test_throws ArgumentError NoHyperfineOneElectronSpecies(;
        mass=-1u"u",
        energies=toy.energies,
        einstein_as=toy.einstein_as,
    )
end

@testitem "Pair queries, angular conversion and physical energy" tags=[:unit, :fast] begin
    using Unitful
    pair = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -1//2)
    B = 0.5u"mT"
    s = Eigenstructure(ca40, "S_1/2", "D_5/2"; B)
    @test transition_frequency(s, pair) ≈
          transition_frequency(ca40, "S_1/2" => "D_5/2") +
          zeeman_shift(s, pair.second) - zeeman_shift(s, pair.first)
    @test zeeman_sensitivity(s, pair) ≈ zeeman_sensitivity(ca40, pair)
    @test angular(1u"MHz") ≈ 2π * 1u"rad/µs"
    @test cyclic(angular(1u"MHz")) ≈ 1u"MHz"
    @test level_energy(ca40, "S_1/2") == 0u"J"
    @test collect(StateBasis(ca40, "S_1/2")) == collect(StateBasis("S_1/2"))
    @test all(p -> p isa Pair, keys(ca40.einstein_as))
    s43 = Eigenstructure(ca43, "S_1/2"; B)
    st = first(Levels.basis(s43))
    @test state_energy(s43, st) ≈
          level_energy(ca43, "S_1/2") + u"ħ" * state_offset(s43, st)
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
    s = Eigenstructure(ca40, basis; B=0.5u"mT")
    scheme = LaserScheme(LaserBeam[])
    @test_throws ArgumentError lindblad_model(s, scheme)
    dropped = lindblad_model(s, scheme; open_channels=:drop)
    @test length(dropped.jump_operators)>0
end

@testitem "The static field enters only through Eigenstructure" tags=[:unit, :fast] begin
    using Levels.PeriodicDynamics
    using Levels.OpticalBloch

    # The working field is specified in exactly one place: no exported function
    # of any of the three modules takes a `B` keyword, apart from the structure
    # constructors, and the representation of a matrix is never a tag.
    offenders = String[]
    for mod in (Levels, Levels.PeriodicDynamics, Levels.OpticalBloch),
        name in names(mod)

        f = getfield(mod, name)
        (f isa Function || f isa Type) || continue
        for m in methods(f)
            kws = Base.kwarg_decl(m)
            if :B in kws && f !== Eigenstructure
                push!(offenders, "$name: B")
            end
            if :coupling_basis in kws || :drive_basis in kws
                push!(offenders, "$name: basis tag")
            end
        end
    end
    @test isempty(offenders)
    @test any(m -> :B in Base.kwarg_decl(m), methods(Eigenstructure))
end

@testitem "Operator matrices use the structure's eigenbasis" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Levels.PeriodicDynamics
    using Unitful

    B = 0.5u"mT"
    s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B)
    probe = StateSpec("S_1/2 F=4", 0) => StateSpec("D_5/2 F=4", 1)
    direction, ε = beam_vectors(π / 2, π / 4)

    # Species forms give coupled-basis operators over the complete coupled
    # basis; structure forms the same operators rotated into the eigenbasis.
    cb = coupled_basis(s)
    C = coupling_matrix(s, "S_1/2" => "D_5/2", 1e6u"W/m^2", ε, direction)
    @test C ≈ to_eigenbasis(
        s,
        coupling_matrix(ca43, cb, "S_1/2" => "D_5/2", 1e6u"W/m^2", ε, direction),
    )
    b = [10.0, 0.0, 0.0]u"µT"
    @test zeeman_hamiltonian(s, b) ≈ to_eigenbasis(s, zeeman_hamiltonian(ca43, cb, b))
    sp = Eigenstructure(ca43, "S_1/2", "P_1/2"; B)
    E1 = relative_couplings(sp, "S_1/2" => "P_1/2", ε, direction)
    @test E1 ≈ to_eigenbasis(
        sp,
        relative_couplings(ca43, coupled_basis(sp), "S_1/2" => "P_1/2", ε, direction),
    )
    @test any(!iszero, E1)

    # A model takes eigenbasis matrices as they are: a hand-written
    # coupled-basis operator enters through to_eigenbasis.
    X = zeros(ComplexF64, length(cb), length(cb)) .* u"rad/µs"
    X[stateindex(cb, probe.second), stateindex(cb, probe.first)] = angular(100u"kHz")
    dt = DrivenTransition(
        s,
        probe;
        drive_frequency=angular(50u"MHz"),
        coupling=to_eigenbasis(s, X),
        drives=zeeman_drives(s, b),
    )
    @test dt.coupling ≈ to_eigenbasis(s, X)
    @test abs(
        dt.coupling[
            stateindex(Levels.basis(s), probe.second),
            stateindex(Levels.basis(s), probe.first),
        ],
    ) < angular(100u"kHz")
    dt_e = DrivenTransition(
        s,
        probe;
        drive_frequency=angular(50u"MHz"),
        coupling=C,
        drives=zeeman_drives(s, b),
    )
    @test dt_e.coupling ≈ C
end

@testitem "Angular frequencies carry rad, rates do not" tags=[:unit, :fast] begin
    using Levels.PeriodicDynamics
    using Levels.OpticalBloch
    using Unitful
    rad_per_µs = u"rad/µs"
    B = 0.5u"mT"
    direction, ε = beam_vectors(π / 2, π / 4)
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    s88 = Eigenstructure(sr88, "S_1/2", "D_5/2"; B)
    s43 = Eigenstructure(ca43, "S_1/2"; B)
    @test unit(angular(1u"MHz")) == rad_per_µs
    @test unit(cyclic(angular(1u"MHz"))) == u"Hz"
    @test unit(transition_frequency(sr88, "S_1/2" => "D_5/2")) == rad_per_µs
    @test unit(transition_frequency(ca43, "S_1/2 F=4" => "D_5/2 F=4")) == rad_per_µs
    @test unit(transition_frequency(s88, probe)) == rad_per_µs
    @test unit(zeeman_shift(s88, probe.first)) == rad_per_µs
    @test unit(zeeman_sensitivity(sr88, probe)) == u"rad/µs/mT"
    @test unit(zeeman_sensitivity(s88, probe)) == u"rad/µs/mT"
    @test unit(Levels.hyperfine_shift(ca43, "S_1/2 F=4")) == rad_per_µs
    @test unit(eltype(offsets(s43))) == rad_per_µs
    @test unit(state_offset(s43, first(Levels.basis(s43)))) == rad_per_µs
    @test unit(state_energy(s43, first(Levels.basis(s43)))) == u"J"
    @test unit(eltype(coupled_hamiltonian(s43))) == rad_per_µs
    @test unit(field(s43)) == u"mT"
    hf_pair = StateSpec("S_1/2 F=4", 0) => StateSpec("S_1/2 F=3", 0)
    @test unit(zeeman_sensitivity(s43, hf_pair)) == u"rad/µs/mT"
    @test unit(rabi_frequency(sr88, probe, 1.0u"W/m^2", ε, direction)) == rad_per_µs
    @test unit(rabi_frequency(s88, probe, 1.0u"W/m^2", ε, direction)) == rad_per_µs
    basis = StateBasis(sr88, "S_1/2", "D_5/2")
    C = coupling_matrix(sr88, basis, "S_1/2" => "D_5/2", 1.0u"W/m^2", ε, direction)
    @test unit(eltype(C)) == rad_per_µs
    @test unit(
        eltype(coupling_matrix(s88, "S_1/2" => "D_5/2", 1.0u"W/m^2", ε, direction)),
    ) == rad_per_µs
    @test unit(eltype(zeeman_hamiltonian(s88, [1.0, 0.0, 0.0]u"µT"))) == rad_per_µs
    @test unit(
        eltype(
            dipole_components(
                Eigenstructure(sr88, "S_1/2", "P_1/2"; B),
                "S_1/2" => "P_1/2",
            )[2],
        ),
    ) == u"C*m"
    @test unit(RelativeFrequency(s88, probe, angular(1u"kHz")).offset) == rad_per_µs
    c = LightShiftCoefficients(s88, 674.0u"nm")
    @test unit(driven_light_shift(c, probe, 1.0u"W/m^2", ε; direction)) == rad_per_µs
    @test unit(light_shift(c, probe.first, 1.0u"W/m^2", ε)) == rad_per_µs
    parked = LightShiftCoefficients(
        s88,
        RelativeFrequency("S_1/2" => "D_5/2", angular(1u"kHz")),
    )
    @test unit(light_shift(parked, probe.first, 1.0u"W/m^2", ε; direction)) ==
          rad_per_µs
    @test unit(
        light_shift(parked, probe.first, 1.0u"W/m^2", ε; direction, parts=:resonant),
    ) == rad_per_µs
    @test unit(polarisability(sr88, probe.first, 0u"rad/s", ε)) == u"C*m^2/V"
    @test unit(polarisability(s88, probe.first, 0u"rad/s", ε)) == u"C*m^2/V"
    # Rates are plain inverse times.
    @test unit(einstein_a(sr88, "S_1/2" => "P_1/2")) == u"µs^-1"
    @test unit(lifetime(sr88, "P_1/2")) == u"s"
    @test unit(saturation_intensity(sr88, "S_1/2" => "P_1/2")) == u"W/m^2"
    # Periodic dynamics.
    dt = DrivenTransition(
        s88,
        probe;
        drive_frequency=angular(50u"MHz"),
        coupling=rabi_normalised(C, basis, probe, angular(100u"kHz")),
        drives=zeeman_drives(s88, [1.0, 0.0, 0.0]u"µT"),
    )
    @test unit(eltype(dt.frame)) == rad_per_µs
    @test unit(eltype(dt.coupling)) == rad_per_µs
    sideband = sideband_rabi(dt; sideband=1)
    @test unit(sideband.Ω) == rad_per_µs && unit(sideband.δ_res) == rad_per_µs
    shift = quasienergy_shift(s88, probe, dt.drive_frequency, dt.drives)
    @test unit(shift.shift) == rad_per_µs
    @test unit(
        ac_zeeman_shift(s88, probe.first, dt.drive_frequency, [1.0, 0.0, 0.0]u"µT"),
    ) == rad_per_µs
    # Optical Bloch: angular Hamiltonian, √rate jump operators.
    beam = LaserBeam("S_1/2" => "P_1/2", angular(-10u"MHz"), 1.0u"W/m^2", ε, direction)
    model = lindblad_model(
        Eigenstructure(sr88, "S_1/2", "P_1/2"; B),
        LaserScheme([beam]);
        open_channels=:drop,
    )
    @test unit(eltype(model.hamiltonian)) == rad_per_µs
    @test unit(eltype(first(model.jump_operators))) == u"µs^(-1/2)"
end
