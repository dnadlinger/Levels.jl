# Tests for the LevelsQuantumToolboxExt package extension: QuantumObject/ket
# conversions over a StateBasis (with the 𝐓ᵖ → time_unitᵖ stripping rule) and
# the fourier_hamiltonians decomposition, checked against full_hamiltonian and
# end-to-end by solving a σ⁻ optical-pumping steady state with
# steadystate_fourier against mesolve time propagation.

@testitem "QuantumToolbox conversions" tags = [:unit] begin
    using LinearAlgebra
    using QuantumToolbox: QuantumObject, basis, expect, projection
    using Unitful

    sb = StateBasis(["S_1/2", "P_1/2"])
    n = length(sb)
    A = randn(ComplexF64, n, n)

    # Operators: quantities of time dimension 𝐓ᵖ are stripped in time_unitᵖ
    # (µsᵖ by default) — angular frequencies as µs⁻¹, Lindblad jump operators
    # as µs^(-1/2), dimensionless matrices unchanged.
    @test QuantumObject(A .* u"ms^-1", sb).data ≈ A ./ 1e3
    @test QuantumObject(A .* u"ms^-1", sb; time_unit=u"ms").data ≈ A
    @test QuantumObject(A .* u"ms^-1", sb; time_unit=u"s").data ≈ 1e3 .* A
    @test QuantumObject(A .* sqrt(1.0u"ms^-1"), sb).data ≈ A ./ sqrt(1e3)
    @test QuantumObject(A, sb).data == A
    @test size(QuantumObject(A, sb)) == (n, n)

    # Kets, both from coefficient vectors and as basis states.
    v = normalize(randn(ComplexF64, n))
    @test QuantumObject(v, sb).data == v
    for (i, s) in enumerate(sb)
        @test basis(sb, s).data == ComplexF64.(1:n .== i)
    end
    @test basis(sb, "P_1/2", 1 // 2) == QuantumObject(ComplexF64.(1:n .== n), sb)

    # Projections address states through stateindex.
    p = projection(sb, StateSpec("S_1/2", 1 // 2), StateSpec("P_1/2", -1 // 2))
    @test p.data == [i == 2 && k == 3 for i in 1:n, k in 1:n]
    @test projection(sb, StateSpec("S_1/2", 1 // 2)) ==
          projection(sb, StateSpec("S_1/2", 1 // 2), StateSpec("S_1/2", 1 // 2))
    ket = basis(sb, "S_1/2", 1 // 2)
    @test expect(projection(sb, StateSpec("S_1/2", 1 // 2)), ket) ≈ 1
    @test abs(expect(projection(sb, StateSpec("S_1/2", -1 // 2)), ket)) < 1e-15

    # Error paths: size mismatches, non-time dimensions, bad reference units.
    @test_throws ArgumentError QuantumObject(zeros(ComplexF64, n + 1, n + 1), sb)
    @test_throws ArgumentError QuantumObject(zeros(ComplexF64, n + 1), sb)
    @test_throws ArgumentError QuantumObject(A .* u"m", sb)
    @test_throws ArgumentError QuantumObject(A, sb; time_unit=u"m")
    @test_throws ArgumentError QuantumObject(A .* u"ms^-1", sb; time_unit=u"µm")
end

@testitem "fourier_hamiltonians vs full_hamiltonian" tags = [:unit] setup =
    [PeriodicSetup] begin
    using QuantumToolbox

    dt = make_transition(
        TRANSITIONS[3];
        drives=zeeman_drives(sr88, BASIS, [12.0 + 5.0im, 5.0 - 3.0im, 8.0im] .* u"µT"),
    )
    δ = 2π * 1.23u"kHz"
    fh = fourier_hamiltonians(dt; δ)

    @test fh.ωd ≈ ustrip(u"µs^-1", Ω_RF)
    @test fh.H_m.data == adjoint(fh.H_p.data)
    @test ishermitian(fh.H_0.data)

    # The three components reassemble the full Hamiltonian at any time.
    for θ in 0.0:0.4:6.0
        t = θ / dt.drive_frequency
        H_t = ustrip.(u"µs^-1", Levels.PeriodicDynamics.full_hamiltonian(dt, t; δ))
        @test fh.H_0.data .+ cis(θ) .* fh.H_p.data .+ cis(-θ) .* fh.H_m.data ≈ H_t
    end

    # The reference time scale is configurable; the µs default matches.
    @test fourier_hamiltonians(dt; δ, time_unit=u"ns").ωd ≈ 1e-3 * fh.ωd
    @test fourier_hamiltonians(dt; δ, time_unit=u"ns").H_0.data ≈ 1e-3 .* fh.H_0.data
    @test fourier_hamiltonians(dt; δ, time_unit=u"µs") == fh
end

@testitem "steadystate_fourier σ⁻ optical pumping" tags = [:integration, :slow] begin
    using LinearAlgebra
    using QuantumToolbox:
        QobjEvo,
        QuantumObject,
        SteadyStateLinearSolver,
        basis,
        expect,
        mesolve,
        projection,
        steadystate_fourier
    using Unitful
    using Levels.PeriodicDynamics

    # The ⁸⁸Sr⁺ σ⁻ 422 nm state-preparation scenario of
    # examples/sr88-ac-zeeman-state-prep, at s = 1: with the trap rf off,
    # |S_1/2, m = −1/2⟩ is exactly dark; the transverse ac Zeeman coupling
    # makes it leaky, so the periodic steady state carries a small error.
    sb = StateBasis(["S_1/2", "P_1/2"])
    dark = StateSpec("S_1/2", -1 // 2)
    bright = StateSpec("S_1/2", 1 // 2)
    excited = StateSpec("P_1/2", -1 // 2)

    # σ⁻ beam along ẑ at saturation intensity: the only coupling out of the
    # S manifold is bright → excited.
    ε = [1.0, -1.0im, 0.0] ./ sqrt(2)
    beam = [0.0, 0.0, 1.0]
    coupling = zeros(ComplexF64, 4, 4) .* u"µs^-1"
    for ml in (-1 // 2, 1 // 2), mu in (-1 // 2, 1 // 2)
        lo, up = StateSpec("S_1/2", ml), StateSpec("P_1/2", mu)
        coupling[stateindex(sb, up), stateindex(sb, lo)] = rabi_frequency(
            sr88,
            lo,
            up,
            saturation_intensity(sr88, "S_1/2", "P_1/2"),
            ε,
            beam,
        )
    end

    # Polarisation-resolved spontaneous emission P_1/2 → S_1/2.
    a_sp = einstein_a(sr88, "S_1/2", "P_1/2")
    c_ops = map(-1:1) do q
        C = zeros(ComplexF64, 4, 4)
        for m in (-1 // 2, 1 // 2)
            abs(m + q) <= 1 // 2 || continue
            C[stateindex(sb, "S_1/2", m), stateindex(sb, "P_1/2", m + q)] =
                Levels.clebsch_gordan(StateSpec("S_1/2", m), StateSpec("P_1/2", m + q))
        end
        QuantumObject(sqrt(a_sp) .* C, sb)
    end

    model(b) = DrivenTransition(
        sr88,
        sb,
        bright => excited;
        drive_frequency=2π * 50.48559u"MHz",
        static_field=0.4955125u"mT",
        coupling,
        drives=zeeman_drives(sr88, sb, b),
    )
    # Direct sparse factorisation instead of the preconditioned-GMRES default:
    # the system is tiny and the cross-checks below resolve ≲1e-10.
    solver = SteadyStateLinearSolver(; alg=nothing)
    p_dark = ρ -> real(expect(projection(sb, dark), ρ))

    # Trap rf off: the dark state is the exact steady state.
    ρs = steadystate_fourier(
        fourier_hamiltonians(model(zeros(3) .* u"µT"))...,
        c_ops;
        n_max=2,
        solver,
    )
    @test p_dark(ρs[1]) ≈ 1 atol = 1e-10

    # Trap rf on: the k = 0 harmonic agrees with mesolve time propagation
    # into the steady state, averaged over one drive period.
    fh = fourier_hamiltonians(model([30.0 + 10.0im, 5.0 - 12.0im, 46.0] .* u"µT"))
    ρs = steadystate_fourier(fh..., c_ops; n_max=6, solver)
    error_fourier = 1 - p_dark(ρs[1])
    @test error_fourier > 1e-5

    H_t = QobjEvo((
        fh.H_0,
        (fh.H_p, (p, t) -> cis(fh.ωd * t)),
        (fh.H_m, (p, t) -> cis(-fh.ωd * t)),
    ))
    period = 2π / fh.ωd
    ts = 10.0 .+ period .* (0:63) ./ 64   # ≫ the ~0.03 µs pumping time constant
    sol = mesolve(
        H_t,
        basis(sb, bright),
        [0.0; ts],
        c_ops;
        e_ops=[projection(sb, dark)],
        progress_bar=Val(false),
        reltol=1e-10,
        abstol=1e-12,
    )
    error_mesolve = 1 - sum(real, sol.expect[1, 2:end]) / 64
    @test error_fourier ≈ error_mesolve rtol = 1e-4
end
