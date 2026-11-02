# Tests for the saturated phonon-number rate equation (rate_equation.jl): its
# weak-sideband limit, its agreement with the full internal ⊗ Fock model where
# adiabatic elimination fails, and the evolution of its populations.

@testitem "Fock displacement element" tags=[:unit, :fast] begin
    using Levels.OpticalBloch: displacement_element, displacement_elements

    for η in (0.3, -0.3)
        elements = displacement_elements(η, 8)
        @test elements ≈ transpose(elements)
        @test [displacement_element(η, n, k) for n in 7:7, k in 0:7] ≈ permutedims(elements[8, 8:-1:1])
    end
    # The first-order expansion of the first sideband.
    @test displacement_element(1e-4, 5, 1) ≈ im * 1e-4 * sqrt(5) rtol = 1e-6
end

@testitem "Phonon rate equation: recovers adiabatic elimination for weak sidebands" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: QuantumToolbox

    ω = angular(0.6u"MHz")
    s = sideband_toy(; γ=2π * 50.0u"kHz", ω, η=0.02, Ω=angular(20.0u"kHz"), δ=(-ω))
    ae = cooling_metrics(s.model, s.mc, AdiabaticElimination())
    rungs = 1:7
    for order in (1, Inf)
        r = PhononRateEquation(s.model, s.mc; num_fock=8, lamb_dicke_order=order)
        @test r.lowering_rates ≈ rungs .* ae.A_minus rtol = 5e-3
        @test r.raising_rates[1] ≈ ae.A_plus rtol = 5e-3
        @test steady_value(r) ≈ ae.nbar rtol = 5e-3
        @test integral_relaxation_time(r, 1.0) ≈ ae.τ_c rtol = 5e-3
    end

    # The ratio departs from A₊/A₋ along the ladder at order η⁴, in the full model too.
    r = PhononRateEquation(s.model, s.mc; num_fock=8, lamb_dicke_order=1)
    full = cooling_metrics(
        s.model,
        s.mc,
        IntegratedTransient(; initial_nbar=0.5, num_fock=14, lamb_dicke_order=1),
    )
    bl = BandLiouvillian(full.model; bandwidth=4)
    populations = fock_populations(bl, full.solver.steady, 1)
    @test r.raising_rates[1:4] ./ r.lowering_rates[1:4] ≈
          populations[2:5] ./ populations[1:4] rtol = 0.02
    @test r.raising_rates[4] / r.lowering_rates[4] > 1.1ae.A_plus / ae.A_minus
end

@testitem "Phonon rate equation: saturated cooling agrees with the full model" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: steadystate

    # Two-level sideband cooling at η Ω √n̄₀ ~ γ.
    ω = angular(0.6u"MHz")
    s = sideband_toy(; γ=2π * 50.0u"kHz", ω, η=0.1, Ω=angular(200.0u"kHz"), δ=(-ω))
    ae = cooling_metrics(s.model, s.mc, AdiabaticElimination())
    for order in (1, Inf)
        saturated = cooling_metrics(
            s.model,
            s.mc,
            SaturatedRates(; initial_nbar=5.0, lamb_dicke_order=order),
        )
        exact = cooling_metrics(
            s.model,
            s.mc,
            IntegratedTransient(; initial_nbar=5.0, lamb_dicke_order=order),
        )
        @test saturated.τ_c ≈ exact.τ_c rtol = 0.05
        @test saturated.nbar ≈ exact.nbar rtol = 0.03
        @test exact.τ_c > 1.5ae.τ_c
    end

    # A multi-level dark resonance, whose occupation adiabatic elimination misses by 20 %.
    mode = MotionalMode(angular(1.1u"MHz"), [0, 0, 1.0])
    scheme = lambda_scheme(
        angular(300.0u"MHz"),
        angular(289.0u"MHz");
        P397=1.0u"µW",
        P866=5.0u"µW",
    )
    model = lindblad_model(LAMBDA, scheme)
    mc = motional_coupling(model, mode)
    saturated = cooling_metrics(
        model,
        mc,
        SaturatedRates(; initial_nbar=1.0, num_fock=32, lamb_dicke_order=1),
    )
    exact = cooling_metrics(
        model,
        mc,
        IntegratedTransient(;
            initial_nbar=1.0,
            num_fock=32,
            lamb_dicke_order=1,
            bandwidth=2,
        ),
    )
    @test saturated.nbar ≈ exact.nbar rtol = 0.02
    @test saturated.τ_c ≈ exact.τ_c rtol = 0.05
    @test exact.nbar > 1.15cooling_metrics(model, mc, AdiabaticElimination()).nbar
end

@testitem "Phonon rate equation: evolution from an initial state" tags=[:integration] setup=[
    OpticalBlochSetup,
] begin
    using QuantumToolbox: steadystate

    ω = angular(0.6u"MHz")
    s = sideband_toy(; γ=2π * 50.0u"kHz", ω, η=0.1, Ω=angular(200.0u"kHz"), δ=(-ω))
    num_fock = fock_truncation(5.0)
    r = PhononRateEquation(s.model, s.mc; num_fock)
    @test sum(fock_populations(r)) ≈ 1
    @test steady_value(r) ≈ sum((0:(num_fock-1)) .* fock_populations(r))

    start = thermal_populations(5.0, num_fock)
    ts = range(0.0, 400.0; length=401)
    curve = cooling_curve(r, 5.0, ts; fock=true)
    @test curve.t == collect(ts) .* u"µs"
    @test size(curve.nbar) == (401, 1)
    @test curve.fock[1][:, 1] ≈ start
    @test vec(cooling_curve(r, 5.0, [0.0, 1e5]; fock=true).fock[1][:, end]) ≈
          fock_populations(r) atol = 1e-6
    @test cooling_curve(r, start, ts .* u"µs").nbar ≈ curve.nbar
    @test isnothing(cooling_curve(r, 5.0, ts).fock)

    # The relaxation time is the normalised area of the curve.
    nbar = vec(curve.nbar)
    area =
        sum(nbar[1:(end-1)] .+ nbar[2:end]) / 2 * step(ts) - steady_value(r) * last(ts)
    @test integral_relaxation_time(r, 5.0) ≈ area / (nbar[1] - steady_value(r)) * u"µs" rtol =
        1e-2

    # The time to a target occupation is where the curve crosses it.
    crossing = time_to_nbar(r, 5.0, 1.0)
    k = findfirst(<(1.0), nbar)
    @test ustrip(u"µs", crossing) ≈ ts[k] atol = step(ts)
    @test time_to_nbar(r, start, 1.0) ≈ crossing
    @test iszero(time_to_nbar(r, 0.5, 1.0))
    @test isinf(time_to_nbar(r, 5.0, 0.1steady_value(r)))

    # … and follows the full model in saturation.
    full = BandLiouvillian(
        motional_model(s.model, s.mc; num_fock, lamb_dicke_order=Inf);
        bandwidth=4,
    )
    ρ_internal = Matrix(steadystate(s.model).data)
    exact = cooling_curve(full, thermal_state(ρ_internal, num_fock, 5.0), ts)
    @test maximum(abs, nbar .- vec(exact.nbar)) < 0.15
end

@testitem "Phonon rate equation: dark Fock states" tags=[:integration] setup=[
    OpticalBlochSetup,
] begin
    using Levels.OpticalBloch: RungLadder, fock_coupling_zeros, rung_nodes

    # The first sideband vanishes near n ≈ 3.67/η² ≈ 41.
    ω = angular(0.6u"MHz")
    η = 0.3
    s = sideband_toy(; γ=2π * 50.0u"kHz", ω, η, Ω=angular(200.0u"kHz"), δ=(-ω))
    @test any(in(39:42), fock_coupling_zeros(η, 100))
    @test isempty(fock_coupling_zeros(0.01, 100))

    num_fock = fock_truncation(20.0)
    resolved = PhononRateEquation(s.model, s.mc; num_fock, node_ratio=1)
    coarse = PhononRateEquation(s.model, s.mc; num_fock)
    rate_per_quantum(r) = r.lowering_rates ./ (1:(num_fock-1))
    @test argmin(rate_per_quantum(resolved)[20:60]) + 19 in 38:44
    @test argmin(rate_per_quantum(coarse)[20:60]) + 19 in 38:44

    ladder =
        RungLadder(s.model, s.mc; lamb_dicke_order=Inf, bandwidth=2, time_unit=u"µs")
    @test length(rung_nodes(ladder, num_fock - 1, 2)) < num_fock - 1
    @test integral_relaxation_time(coarse, 20.0) ≈
          integral_relaxation_time(resolved, 20.0) rtol = 0.08
end

@testitem "Phonon rate equation: heating, several modes and validation" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    ω = angular(0.6u"MHz")
    s = sideband_toy(; γ=2π * 50.0u"kHz", ω, η=0.05, Ω=angular(100.0u"kHz"), δ=(-ω))
    method = SaturatedRates(; initial_nbar=1.0, num_fock=12)
    result = cooling_metrics(s.model, s.mc, method)
    @test result isa CoolingMetrics && result.method == method
    @test result.solver isa PhononRateEquation
    @test isnothing(result.model) &&
          isnothing(result.A_plus) &&
          isnothing(result.A_minus)
    @test result.nbar == steady_value(result.solver)
    @test result.τ_c == integral_relaxation_time(result.solver, 1.0)
    @test SaturatedRates(; initial_nbar=2.0).num_fock === nothing
    @test cooling_metrics(s.model, s.mc, SaturatedRates(; initial_nbar=2.0)).solver.lowering_rates |>
          length == fock_truncation(2.0) - 1

    # Heating adds to both rates in proportion to n; other time units are equivalent.
    heating = 1000.0u"s^-1"
    heated = PhononRateEquation(s.model, s.mc; num_fock=12, heating_rate=heating)
    @test heated.lowering_rates ≈ result.solver.lowering_rates .+ (1:11) .* heating
    @test heated.raising_rates ≈ result.solver.raising_rates .+ (1:11) .* heating
    @test steady_value(heated) > result.nbar
    in_ms = cooling_metrics(s.model, s.mc, method; time_unit=u"ms")
    @test in_ms.τ_c ≈ result.τ_c rtol = 1e-8
    @test unit(in_ms.τ_c) == u"ms"

    # Modes are independent.
    modes = [s.mode, MotionalMode(angular(1.2u"MHz"), Z_AXIS)]
    mcs = motional_coupling(s.model, modes; recoil_moment=:isotropic)
    many = cooling_metrics(
        s.model,
        mcs,
        SaturatedRates(; initial_nbar=[1.0, 0.5], num_fock=[12, 6]);
        heating_rate=[0.0u"s^-1", heating],
    )
    @test many[1].τ_c == result.τ_c && many[1].nbar == result.nbar
    @test length(many[2].solver.lowering_rates) == 5
    @test_throws ArgumentError cooling_metrics(
        s.model,
        mcs,
        SaturatedRates(; initial_nbar=[1.0], num_fock=4),
    )

    @test_throws ArgumentError SaturatedRates(; initial_nbar=-1.0)
    @test_throws ArgumentError SaturatedRates(; initial_nbar=1.0, lamb_dicke_order=2)
    @test_throws ArgumentError SaturatedRates(; initial_nbar=1.0, bandwidth=0)
    @test_throws ArgumentError SaturatedRates(; initial_nbar=1.0, node_ratio=0.5)
    @test_throws ArgumentError PhononRateEquation(s.model, s.mc; num_fock=1)

    # A constant recoil moment only suffices at first order.
    constant = motional_coupling(s.model, s.mode; recoil_moment=0.4)
    @test_throws ArgumentError PhononRateEquation(s.model, constant; num_fock=4)
    @test PhononRateEquation(s.model, constant; num_fock=4, lamb_dicke_order=1) isa
          PhononRateEquation

    # Beat notes have no static rate equation.
    two_tones = LaserScheme(
        [
            LaserBeam(
                "S_1/2" => "P_1/2",
                Δ,
                peak_intensity(10u"µW", 75u"µm"),
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ) for Δ in (angular(100.0u"MHz"), angular(120.0u"MHz"))
        ];
        frame_reference="S_1/2",
    )
    beating = lindblad_model(LAMBDA, two_tones; open_channels=:drop)
    @test !isempty(beating.harmonics)
    @test_throws ArgumentError PhononRateEquation(
        beating,
        motional_coupling(beating, s.mode);
        num_fock=4,
    )
end
