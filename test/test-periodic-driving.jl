@testitem "DrivenTransition assembly" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful
    using Levels.PeriodicDynamics
    using Levels.PeriodicDynamics: full_hamiltonian

    basis = StateBasis(["S_1/2", "D_5/2"])
    B0 = 0.5u"mT"
    Ω_rf = 2π * 50.0u"MHz"
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    Ω0 = 2π * 100.0u"kHz"

    n, ε = beam_vectors(π / 2, π / 4)
    L = rabi_normalised(
        quadrupole_couplings(basis, "S_1/2", "D_5/2", ε, n),
        basis,
        probe,
        Ω0,
    )
    drive = zeeman_drive(sr88, basis, [0.0, 0.0, 2.5]u"µT"; phase=π / 2)

    dt = DrivenTransition(
        sr88,
        basis,
        probe;
        drive_frequency=Ω_rf,
        static_field=B0,
        coupling=L,
        drives=[drive],
    )

    # Frame diagonal vanishes on the probed pair and reproduces relative Zeeman
    # shifts elsewhere.
    @test iszero(dt.frame[stateindex(basis, probe.first)])
    @test iszero(dt.frame[stateindex(basis, probe.second)])
    i = stateindex(basis, "S_1/2", -1//2)
    @test dt.frame[i] ≈
          zeeman_shift(sr88, StateSpec("S_1/2", -1//2), B0) -
          zeeman_shift(sr88, probe.first, B0)

    # H(t) is Hermitian, and the probed off-diagonal element is f(t)·L.
    t = 3.7u"ns"
    modulation = HarmonicPhaseModulation(0.05, -0.02)
    H = full_hamiltonian(dt, t; δ=2π * 1.0u"kHz", modulation)
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
    δ = 2π * 1.0u"kHz"
    Hδ = full_hamiltonian(dt, t; δ, modulation)
    for k in 1:length(basis)
        expected = k in staterange(basis, "D_5/2") ? -δ : 0.0u"kHz"
        # (atol set by float cancellation against the MHz-scale diagonal)
        @test Hδ[k, k] - H0[k, k] ≈ expected atol = 1e-9u"kHz"
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
    @test Hstart[il, il] ≈ 0.0u"MHz" atol = 1e-12u"MHz"
    @test Hstart[iu, il] ≈ L[iu, il] / 2
end

@testitem "DrivenTransition validation" tags=[:unit, :fast] begin
    using Unitful
    using Levels.PeriodicDynamics

    basis = StateBasis(["S_1/2", "D_5/2"])
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    n, ε = beam_vectors(π / 2, π / 4)
    L = rabi_normalised(
        quadrupole_couplings(basis, "S_1/2", "D_5/2", ε, n),
        basis,
        probe,
        2π * 100.0u"kHz",
    )
    kwargs = (drive_frequency=2π * 50.0u"MHz", static_field=0.5u"mT", coupling=L)

    # The basis must consist of exactly the probed levels.
    wide = StateBasis(["S_1/2", "D_3/2", "D_5/2"])
    n2, ε2 = n, ε
    Lwide = rabi_normalised(
        quadrupole_couplings(wide, "S_1/2", "D_5/2", ε2, n2),
        wide,
        probe,
        2π * 100.0u"kHz",
    )
    @test_throws ArgumentError DrivenTransition(
        sr88,
        wide,
        probe;
        kwargs...,
        coupling=Lwide,
    )

    # Probed states must be in different levels.
    @test_throws ArgumentError DrivenTransition(
        sr88,
        basis,
        StateSpec("D_5/2", 1//2) => StateSpec("D_5/2", 5//2);
        kwargs...,
    )

    # Couplings outside the upper⟨row|lower⟩⟨col| block are rejected.
    bad = copy(L)
    bad[stateindex(basis, probe.first), stateindex(basis, probe.second)] = 1.0u"kHz"
    @test_throws ArgumentError DrivenTransition(
        sr88,
        basis,
        probe;
        kwargs...,
        coupling=bad,
    )

    # Harmonic drives must be block-diagonal in the levels.
    amplitude = zeros(ComplexF64, 8, 8) * u"MHz"
    amplitude[3, 1] = amplitude[1, 3] = 1.0u"MHz"
    @test_throws ArgumentError DrivenTransition(
        sr88,
        basis,
        probe;
        kwargs...,
        drives=[HarmonicDrive(amplitude, 0.0)],
    )
end
