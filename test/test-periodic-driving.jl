@testitem "DrivenTransition assembly" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful
    using Levels.PeriodicDynamics
    using Levels.PeriodicDynamics: full_hamiltonian

    basis = StateBasis(["S_1/2", "D_5/2"])
    B0 = 0.5u"mT"
    s = Eigenstructure(sr88, basis; B=B0)
    Ω_rf = angular(50.0u"MHz")
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    Ω0 = angular(100.0u"kHz")

    n, ε = beam_vectors(π / 2, π / 4)
    L = rabi_normalised(
        relative_couplings(s, "S_1/2" => "D_5/2", ε, n),
        basis,
        probe,
        Ω0,
    )
    drive = zeeman_drive(s, [0.0, 0.0, 2.5]u"µT"; phase=π / 2)

    dt = DrivenTransition(s, probe; drive_frequency=Ω_rf, coupling=L, drives=[drive])
    @test dt.basis === Levels.basis(s)

    # Frame diagonal vanishes on the probed pair and reproduces relative Zeeman
    # shifts elsewhere.
    @test iszero(dt.frame[stateindex(basis, probe.first)])
    @test iszero(dt.frame[stateindex(basis, probe.second)])
    i = stateindex(basis, "S_1/2", -1//2)
    @test dt.frame[i] ≈
          zeeman_shift(s, StateSpec("S_1/2", -1//2)) - zeeman_shift(s, probe.first)

    # H(t) is Hermitian, and the probed off-diagonal element is f(t)·L.
    t = 3.7u"ns"
    modulation = HarmonicPhaseModulation(0.05, -0.02)
    H = full_hamiltonian(dt, t; δ=angular(1.0u"kHz"), modulation)
    @test H ≈ H'
    φ_mod =
        0.05 * cos(uconvert(NoUnits, Ω_rf * t)) -
        0.02 * sin(uconvert(NoUnits, Ω_rf * t))
    f = cis(-φ_mod) / 2
    iu = stateindex(basis, probe.second)
    il = stateindex(basis, probe.first)
    @test H[iu, il] ≈ f * L[iu, il]

    # The detuning acts on the upper manifold only.
    H0 = full_hamiltonian(dt, t; modulation)
    δ = angular(1.0u"kHz")
    Hδ = full_hamiltonian(dt, t; δ, modulation)
    for k in 1:length(basis)
        expected = k in staterange(basis, "D_5/2") ? -δ : 0.0u"rad/ms"
        # (atol set by float cancellation against the µs^-1-scale diagonal)
        @test Hδ[k, k] - H0[k, k] ≈ expected atol = 1e-9u"rad/ms"
    end

    # A callable phase modulation reproduces the harmonic one.
    same = full_hamiltonian(
        dt,
        t;
        modulation=τ ->
            0.05 * cos(uconvert(NoUnits, Ω_rf * τ)) -
            0.02 * sin(uconvert(NoUnits, Ω_rf * τ)),
    )
    @test same ≈ full_hamiltonian(dt, t; modulation)

    # At t = 0 with φ = π/2 the harmonic drive contributes nothing.
    Hstart = full_hamiltonian(dt, 0.0u"s")
    @test Hstart[il, il] ≈ 0.0u"rad/µs" atol = 1e-12u"rad/µs"
    @test Hstart[iu, il] ≈ L[iu, il] / 2

    # The drives are the Zeeman Hamiltonian of the structure.
    @test drive.amplitude == zeeman_hamiltonian(s, [0.0, 0.0, 2.5]u"µT")
    b = [12.0 + 5.0im, 5.0 - 3.0im, 8.0im] .* u"µT"
    drives = zeeman_drives(s, b)
    @test length(drives) == 2
    @test drives[1].amplitude ≈ zeeman_hamiltonian(s, real.(b)) &&
          drives[1].phase == 0.0
    @test drives[2].amplitude ≈ zeeman_hamiltonian(s, imag.(b)) &&
          drives[2].phase ≈ -π / 2
    @test length(zeeman_drives(s, [1.0, 0.0, 0.0]u"µT")) == 1
end

@testitem "DrivenTransition validation" tags=[:unit, :fast] begin
    using Unitful
    using Levels.PeriodicDynamics

    basis = StateBasis(["S_1/2", "D_5/2"])
    s = Eigenstructure(sr88, basis; B=0.5u"mT")
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    n, ε = beam_vectors(π / 2, π / 4)
    L = rabi_normalised(
        relative_couplings(s, "S_1/2" => "D_5/2", ε, n),
        basis,
        probe,
        angular(100.0u"kHz"),
    )
    kwargs = (drive_frequency=angular(50.0u"MHz"), coupling=L)

    # The structure must consist of exactly the probed manifolds…
    wide = Eigenstructure(sr88, ["S_1/2", "D_3/2", "D_5/2"]; B=0.5u"mT")
    Lwide = rabi_normalised(
        relative_couplings(wide, "S_1/2" => "D_5/2", ε, n),
        Levels.basis(wide),
        probe,
        angular(100.0u"kHz"),
    )
    @test_throws ArgumentError DrivenTransition(wide, probe; kwargs..., coupling=Lwide)
    # …and contain them completely.
    partial = Eigenstructure(
        sr88,
        [StateSpec("S_1/2", 1//2), StateSpec("D_5/2", 5//2)];
        B=0.5u"mT",
    )
    @test_throws ArgumentError DrivenTransition(
        partial,
        probe;
        drive_frequency=angular(50.0u"MHz"),
        coupling=[0.0 0.0; 1.0 0.0] .* u"rad/µs",
    )

    # Probed states must be in different levels.
    @test_throws ArgumentError DrivenTransition(
        s,
        StateSpec("D_5/2", 1//2) => StateSpec("D_5/2", 5//2);
        kwargs...,
    )

    # Couplings outside the upper⟨row|lower⟩⟨col| block are rejected.
    bad = copy(L)
    bad[stateindex(basis, probe.first), stateindex(basis, probe.second)] = 1.0u"rad/ms"
    @test_throws ArgumentError DrivenTransition(s, probe; kwargs..., coupling=bad)

    # Harmonic drives must be block-diagonal in the levels.
    amplitude = zeros(ComplexF64, 8, 8) * u"rad/µs"
    amplitude[3, 1] = amplitude[1, 3] = 1.0u"rad/µs"
    @test_throws ArgumentError DrivenTransition(
        s,
        probe;
        kwargs...,
        drives=[HarmonicDrive(amplitude, 0.0)],
    )

    # A probe of the wrong level kind, and the static field as a keyword, are
    # rejected.
    @test_throws ArgumentError DrivenTransition(
        s,
        StateSpec("S_1/2 F=4", 4) => StateSpec("D_5/2 F=4", 3);
        kwargs...,
    )
    @test_throws MethodError DrivenTransition(s, probe; kwargs..., B=0.5u"mT")
end
