# Tests for the semiclassical cooling rates (semiclassical.jl): the Lamb–Dicke limit,
# the agreement of its two engines and of the Doppler gauge with the explicit beam
# phases, the classical limit of the Fock ladder, and the thermal sampling.

@testitem "Semiclassical: participation scales the Lamb–Dicke factors" tags=[
    :unit,
    :fast,
] setup=[OpticalBlochSetup] begin
    ω = angular(1.0u"MHz")
    direction = [1.0, 1.0, 1.0]
    k = [0.0, 0.0, 1.0] .* (2π / 411u"nm")
    plain = MotionalMode(ω, direction)
    crystal = MotionalMode(ω, direction; participation=0.7)
    @test plain.participation == 1
    for projected in (true, false)
        @test lamb_dicke(sr88, crystal, k; projected) ≈
              0.7 * lamb_dicke(sr88, plain, k; projected)
    end
    @test_throws ArgumentError MotionalMode(ω, direction; participation=0)
    @test_throws ArgumentError MotionalMode(ω, direction; participation=-0.5)
    @test_throws ArgumentError MotionalMode(ω, direction; participation=Inf)

    # The recoil factors of a mode scale alike.
    model = lambda_model(angular(-20.0u"MHz"), angular(20.0u"MHz"))
    mc_plain = motional_coupling(model, plain)
    mc_crystal = motional_coupling(model, crystal)
    @test mc_crystal.projected_lamb_dicke ≈ 0.7 .* mc_plain.projected_lamb_dicke
    @test mc_crystal.recoil_lamb_dicke ≈ 0.7 .* mc_plain.recoil_lamb_dicke
end

@testitem "Semiclassical: the limit of small occupation is adiabatic elimination" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: QuantumToolbox

    # Doppler cooling (Lamb–Dicke regime, A₊ comparable to A₋) …
    ω = angular(1.0u"MHz")
    s = sideband_toy(;
        γ=angular(20.0u"MHz"),
        ω,
        η=0.05,
        Ω=angular(5.0u"MHz"),
        δ=angular(-10.0u"MHz"),
    )
    ae = cooling_metrics(s.model, s.mc, AdiabaticElimination())
    W = only(semiclassical_cooling_rates(s.model, [s.mc], [1e-4]))
    @test W ≈ ae.A_minus - ae.A_plus rtol = 1e-3

    # … a multi-level dark-resonance scheme with both beam projections …
    model = lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"))
    mode = MotionalMode(angular(1.2u"MHz"), [1.0, 0.0, 1.0])
    mc = motional_coupling(model, mode)
    ae = cooling_metrics(model, mc, AdiabaticElimination())
    W = only(semiclassical_cooling_rates(model, [mc], [1e-4]))
    @test W ≈ ae.A_minus - ae.A_plus rtol = 1e-3

    # … and through the cooling_metrics front-end, without a steady occupation.
    metrics = cooling_metrics(model, mc, SemiclassicalRates(; occupation=1e-4))
    @test 1 / metrics.τ_c ≈ W rtol = 1e-12
    @test isnan(metrics.nbar)
    @test_throws ArgumentError cooling_metrics(
        model,
        mc,
        SemiclassicalRates(; occupation=1.0);
        heating_rate=1.0u"s^-1",
    )

    # Independent modes at low occupation: each mode at its own limit.
    other = MotionalMode(angular(2.03u"MHz"), [0.0, 1.0, 1.0])
    mcs = motional_coupling(model, [mode, other])
    limits = [
        (m=cooling_metrics(model, c, AdiabaticElimination()); m.A_minus - m.A_plus)
        for c in mcs
    ]
    W = semiclassical_cooling_rates(
        model,
        mcs,
        [1e-2, 1e-2];
        engine=TimeAverage(; averaging_time=40.0u"µs"),
    )
    @test W[1] ≈ limits[1] rtol = 0.02
    @test W[2] ≈ limits[2] rtol = 0.02
end

@testitem "Semiclassical: the engines agree, for hot motion too" tags=[:integration] setup=[
    OpticalBlochSetup,
] begin
    using Levels.OpticalBloch: SemiclassicalProblem, trajectory_rates

    model = lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"))
    mode = MotionalMode(angular(1.2u"MHz"), [1.0, 0.0, 1.0])
    mc = motional_coupling(model, mode)
    p = SemiclassicalProblem(model, [mc])
    # (A single mode needs no phase average, and an integer number of periods no window.)
    periods = 40 * 2π / p.frequencies[1]
    engine = TimeAverage(;
        settling_time=10.0u"µs",
        averaging_time=periods * u"µs",
        num_phases=1,
        window=:rectangular,
        reltol=1e-9,
        abstol=1e-11,
    )
    for n in (0.5, 30.0, 3000.0)
        R_balance = only(semiclassical_phonon_rates(model, [mc], [n]))
        R_average = only(trajectory_rates(p, [n], [0.7], engine)) * u"µs^-1"
        @test R_balance ≈ R_average rtol = 2e-4
    end
    # A larger cutoff margin changes nothing.
    @test semiclassical_phonon_rates(model, [mc], [3000.0]) ≈
          semiclassical_phonon_rates(
        model,
        [mc],
        [3000.0];
        engine=HarmonicBalance(; margin=48),
    ) rtol = 1e-8
    @test_throws ArgumentError semiclassical_phonon_rates(
        model,
        motional_coupling(model, [mode, mode]),
        [1.0, 1.0];
        engine=HarmonicBalance(),
    )
end

@testitem "Semiclassical: the Doppler gauge matches the explicit beam phases" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using LinearAlgebra
    using Levels.OpticalBloch: ODEProblem, Vern7, solve

    # Direct integration of Eqs. (2) and (6) of Bartolotta et al.: the frame diagonal,
    # each beam's C e^{iφ_b(t)}, and the force operator with the same phases.
    model = lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"))
    mode = MotionalMode(angular(1.2u"MHz"), [1.0, 0.0, 1.0])
    mc = motional_coupling(model, mode)
    n = 400.0
    ω = ustrip(u"µs^-1", mode.frequency)
    N = length(model.basis)
    H0 = Diagonal(diag(ComplexF64.(ustrip.(u"µs^-1", model.hamiltonian))))
    Cs = [ComplexF64.(ustrip.(u"µs^-1", C)) for C in model.couplings]
    ηs = mc.projected_lamb_dicke
    Ls = [ComplexF64.(ustrip.(u"µs^(-1/2)", L)) for L in model.jump_operators]
    decay = sum(L' * L for L in Ls)
    function rhs!(du, u, _, t)
        ρ = reshape(view(u, 1:(N^2)), N, N)
        H = Matrix{ComplexF64}(H0)
        force = zeros(ComplexF64, N, N)
        φ = 2 * sqrt(n) * cos(ω * t)
        for (C, η) in zip(Cs, ηs)
            X = C .* cis(η * φ)
            H .+= (X .+ X') ./ 2
            force .+= η .* (im .* X .+ (im .* X)') ./ 2
        end
        dρ = -im .* (H * ρ .- ρ * H) .- (decay * ρ .+ ρ * decay) ./ 2
        for L in Ls
            dρ .+= L * ρ * L'
        end
        du[1:(N^2)] .= vec(dρ)
        du[N^2+1] = 2 * sqrt(n) * sin(ω * t) * real(tr(force * ρ))
    end
    u0 = zeros(ComplexF64, N^2 + 1)
    u0[1:(N+1):(N^2)] .= 1 / N
    settle = 60 * 2π / ω
    averaging = 40 * 2π / ω
    settled = solve(
        ODEProblem(rhs!, u0, (0.0, settle)),
        Vern7();
        reltol=1e-9,
        abstol=1e-11,
        save_everystep=false,
    )
    u1 = copy(last(settled.u))
    u1[end] = 0
    final = solve(
        ODEProblem(rhs!, u1, (settle, settle + averaging)),
        Vern7();
        reltol=1e-9,
        abstol=1e-11,
        save_everystep=false,
    )
    R_direct = real(last(final.u)[end]) / averaging
    R = ustrip(u"µs^-1", only(semiclassical_phonon_rates(model, [mc], [n])))
    @test R ≈ R_direct rtol = 5e-4

    # The rates are independent of the frame reference and of a mode sign.
    reference =
        lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"); frame_reference="S_1/2")
    @test semiclassical_phonon_rates(
        reference,
        [motional_coupling(reference, mode)],
        [n],
    ) ≈ semiclassical_phonon_rates(model, [mc], [n]) rtol = 1e-8
    flipped = MotionalMode(mode.frequency, -mode.direction)
    @test semiclassical_phonon_rates(model, [motional_coupling(model, flipped)], [n]) ≈
          semiclassical_phonon_rates(model, [mc], [n]) rtol = 1e-8
end

@testitem "Semiclassical: the classical limit of the Fock ladder" tags=[:unit, :fast] setup=[
    OpticalBlochSetup,
] begin
    using Levels.OpticalBloch: RungLadder, ladder_net_rate

    # At large n the ladder's elements ⟨n-k|D(η)|n⟩ tend to iᵏJ_k(2η√n), the harmonics
    # of the classical phase, so that its net rate Γ↑ - Γ↓ becomes R(n).
    model = lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"))
    mode = MotionalMode(angular(1.2u"MHz"), [1.0, 0.0, 1.0])
    mc = motional_coupling(model, mode)
    # This needs bandwidth ≳ β = 2 max|η| √n harmonics (3.2 and 10.2 here).
    for (n, bandwidth) in ((100, 10), (1000, 20))
        ladder = RungLadder(model, mc; lamb_dicke_order=Inf, bandwidth, time_unit=u"µs")
        R = ustrip(u"µs^-1", only(semiclassical_phonon_rates(model, [mc], [n - 0.5])))
        @test ladder_net_rate(ladder, n) ≈ R rtol = 0.01
    end
end

@testitem "Semiclassical: thermal sampling" tags=[:unit, :fast] setup=[
    OpticalBlochSetup,
] begin
    using Levels.OpticalBloch: gausslegendre
    using Random: Xoshiro

    model = lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"))
    mode = MotionalMode(angular(1.2u"MHz"), [1.0, 0.0, 1.0])
    mc = motional_coupling(model, mode)
    nbar = 200.0
    # W_th = -⟨R⟩/n̄ over the exponential P function, with n = -n̄ ln u for uniform u.
    nodes, weights = gausslegendre(60)
    R = [
        ustrip(
            u"µs^-1",
            only(semiclassical_phonon_rates(model, [mc], [-nbar * log((x + 1) / 2)])),
        ) for x in nodes
    ]
    W_exact = -sum(weights .* R) / 2 / nbar
    sampled =
        thermal_cooling_rates(model, [mc], [nbar]; num_samples=400, rng=Xoshiro(1))
    W = ustrip(u"µs^-1", only(sampled.rates))
    σ = ustrip(u"µs^-1", only(sampled.standard_errors))
    @test abs(W - W_exact) < 4σ
    @test σ < 0.2 * abs(W_exact)
    # Reproducible given the generator.
    again = thermal_cooling_rates(model, [mc], [nbar]; num_samples=400, rng=Xoshiro(1))
    @test again.rates == sampled.rates
end
