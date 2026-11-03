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

@testitem "Phonon rate equation: closure diagnostics" tags=[:unit, :fast] begin
    using Levels.OpticalBloch
    using Levels.OpticalBloch: close_rung_rates
    using Unitful

    diffusion = close_rung_rates(-1.0, 4.0)
    @test diffusion.rates == (2.5, 1.5)
    @test diffusion.source == :diffusion && !diffusion.clipped
    @test diffusion.window_status == :not_requested
    @test isnothing(diffusion.window_ratio)

    window = close_rung_rates(-1.0, 4.0, 0.2)
    @test last(window.rates) / first(window.rates) ≈ 0.2
    @test window.source == :window && window.window_status == :accepted
    @test !window.clipped && window.total_rate != window.diffusion

    rejected = close_rung_rates(-1.0, 4.0, 2.0)
    @test rejected.rates == diffusion.rates
    @test rejected.source == :diffusion && rejected.window_status == :incompatible
    @test rejected.window_ratio == 2.0

    clipped = close_rung_rates(-1.0, 0.3, 2.0)
    @test clipped.clipped && clipped.window_status == :incompatible
    @test all(>(0), clipped.rates)
    @test last(clipped.rates) - first(clipped.rates) ≈ clipped.net_rate
    @test clipped.total_rate > abs(clipped.net_rate) > clipped.diffusion
    @test close_rung_rates(-1.0, 4.0, Inf).window_status == :invalid
    @test_throws ArgumentError close_rung_rates(NaN, 4.0)

    explicit = PhononRateEquation([2.0]u"µs^-1", [1.0]u"µs^-1", u"µs")
    @test isempty(explicit.diagnostics)
end

@testitem "Phonon rate equation: recoil diffusion retains the motional phase" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using SparseArrays
    using Levels.OpticalBloch:
        RungLadder, ladder_cumulants, ladder_hamiltonians, commutator_superoperator

    ω = angular(0.6u"MHz")
    γ = 2π * 50.0u"kHz"
    s = sideband_toy(; γ, ω, η=0.2, Ω=angular(200.0u"kHz"), δ=(-angular(25.0u"kHz")))
    n = 40
    K = 8
    ladder =
        RungLadder(s.model, s.mc; lamb_dicke_order=Inf, bandwidth=K, time_unit=u"µs")
    hamiltonians = ladder_hamiltonians(ladder, n)
    num_states = ladder.num_states
    num_orders = 2K + 1
    state_identity = sparse(ComplexF64(1) * I, num_states, num_states)
    same_order = sparse(1.0I, num_orders, num_orders)
    shift(k) = spdiagm(k => ones(num_orders - abs(k)))

    # Independently tilt ket and bra motion by ±χ/2. The exact frozen-amplitude
    # recoil factor per kick is exp(-4i η√n sin(χ/2) sin(θ)); symmetric emission
    # makes it real. Its Fourier coefficients retain emission/phase correlations.
    function tilted_generator(χ; phase_averaged=false)
        L =
            kron(
                same_order,
                ladder.dissipator - im * commutator_superoperator(first(hamiltonians)),
            ) + kron(
                spdiagm(0 => -im * ladder.mode_frequency .* ((-K):K)),
                sparse(I, num_states^2, num_states^2),
            )
        for k in 1:K
            H = sparse(hamiltonians[k+1])
            left = -im * kron(state_identity, H)
            right = im * kron(sparse(transpose(H)), state_identity)
            L += kron(shift(k), cis(-χ * k / 2) * left + cis(χ * k / 2) * right)
            L += kron(shift(-k), cis(χ * k / 2) * left + cis(-χ * k / 2) * right)
        end
        phases = 2π .* (0:127) ./ 128
        for (c, kicks, weights) in
            zip(ladder.jump_operators, ladder.recoil_kicks, ladder.recoil_weights)

            isempty(kicks) && continue
            factors = [
                sum(weights .* cos.(4sqrt(n) * sin(χ / 2) * sin(θ) .* kicks)) for
                θ in phases
            ]
            phase_averaged && fill!(factors, sum(factors) / length(factors))
            for k in (-2K):2K
                coefficient = sum(factors .* cis.(-k .* phases)) / length(phases)
                k == 0 && (coefficient -= 1)
                L += coefficient * kron(shift(k), kron(conj(c), c))
            end
        end
        L
    end
    χ = 1e-3
    function cumulants_from_eigenvalue(; kwargs...)
        values = eigvals(Matrix(tilted_generator(χ; kwargs...)))
        λ = values[argmin(abs.(values))]
        (; net_rate=imag(λ) / χ, diffusion=-2real(λ) / χ^2)
    end
    counted = cumulants_from_eigenvalue()
    averaged = cumulants_from_eigenvalue(; phase_averaged=true)
    cumulants = ladder_cumulants(ladder, n)
    @test cumulants.net_rate ≈ counted.net_rate rtol = 1e-4
    @test cumulants.diffusion ≈ counted.diffusion rtol = 1e-4
    @test abs(counted.diffusion - averaged.diffusion) > 0.05counted.diffusion
end

@testitem "Phonon rate equation: recovers adiabatic elimination for weak sidebands" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: QuantumToolbox

    using Levels.OpticalBloch: RungLadder, ladder_cumulants

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
        @test all(
            d -> d.source == :window && d.window_status == :accepted,
            r.diagnostics,
        )
        @test all(
            d -> r.lowering_rates[d.n] + r.raising_rates[d.n] ≈ d.total_rate,
            r.diagnostics,
        )
        # The ladder's current and diffusion are those of Poissonian steps n A₊, n A₋.
        ladder = RungLadder(
            s.model,
            s.mc;
            lamb_dicke_order=order,
            bandwidth=2,
            time_unit=u"µs",
        )
        for n in (1, 5)
            c = ladder_cumulants(ladder, n)
            @test c.net_rate * u"µs^-1" ≈ n * (ae.A_plus - ae.A_minus) rtol = 5e-3
            @test c.diffusion * u"µs^-1" ≈ n * (ae.A_plus + ae.A_minus) rtol = 5e-3
        end
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

@testitem "Phonon rate equation: the spread of a hot state follows the full model" tags=[
    :integration,
    :slow,
] setup=[OpticalBlochSetup] begin
    using QuantumToolbox: steadystate
    using Levels.OpticalBloch: motional_model

    # Beyond a sideband reach of one (n > 6 here), the total rate of a rung is the
    # ladder's diffusion. From a Fock state, the mean and the variance of the
    # populations then follow the full model, where the window's ratio would
    # overstate the variance by 20 %.
    ω = angular(0.6u"MHz")
    s = sideband_toy(; γ=2π * 50.0u"kHz", ω, η=0.2, Ω=angular(200.0u"kHz"), δ=(-ω))
    num_fock = 70
    r = PhononRateEquation(s.model, s.mc; num_fock, bandwidth=6, node_ratio=1.25)
    @test any(d -> d.source == :diffusion && isnothing(d.window_ratio), r.diagnostics)
    moments(p) = (
        nbar=sum((0:(num_fock-1)) .* p);
        (nbar, sum((0:(num_fock-1)) .^ 2 .* p) - nbar^2)
    )
    start = zeros(num_fock)
    start[41] = 1
    ts = [0.0, 63.0, 126.0, 190.0]
    curve = cooling_curve(r, start, ts; fock=true)
    bl = BandLiouvillian(
        motional_model(s.model, s.mc; num_fock, lamb_dicke_order=Inf);
        bandwidth=6,
    )
    ρ_internal = Matrix(steadystate(s.model).data)
    exact = cooling_curve(bl, kron(ρ_internal, Diagonal(start)), ts; fock=true)
    for j in 2:4
        nbar, variance = moments(curve.fock[1][:, j])
        nbar_exact, variance_exact = moments(exact.fock[1][:, j])
        @test nbar ≈ nbar_exact rtol = 0.03
        @test variance ≈ variance_exact rtol = 0.08
    end
end

@testitem "Phonon rate equation: hot rungs are smooth and interpolate" tags=[
    :integration,
] setup=[OpticalBlochSetup] begin
    # A Doppler-cooled Λ system up to n = 1000 (reach 2|η|√n ≈ 10): the net rates tend to
    # the semiclassical drift, the total per quantum has no spikes along the ladder, and
    # the default node spacing interpolates both to a few percent.
    model = lambda_model(angular(-20.0u"MHz"), angular(-5.0u"MHz"))
    mode = MotionalMode(angular(1.2u"MHz"), [1.0, 0.0, 1.0])
    mc = motional_coupling(model, mode)
    resolved =
        PhononRateEquation(model, mc; num_fock=1001, bandwidth=20, node_ratio=1.05)
    coarse = PhononRateEquation(model, mc; num_fock=1001, bandwidth=20)
    net(r) = r.raising_rates .- r.lowering_rates
    total(r) = r.raising_rates .+ r.lowering_rates
    hot = 100:1000
    semiclassical =
        only.(semiclassical_phonon_rates(model, [mc], [n - 0.5]) for n in (100, 1000))
    @test net(resolved)[100] ≈ semiclassical[1] rtol = 0.02
    @test net(resolved)[1000] ≈ semiclassical[2] rtol = 0.02
    @test maximum(abs.(net(coarse)[hot] ./ net(resolved)[hot] .- 1)) < 0.05
    @test maximum(abs.(total(coarse)[hot] ./ total(resolved)[hot] .- 1)) < 0.05
    per_quantum = total(resolved) ./ (1:1000)
    @test maximum(per_quantum[2:end] ./ per_quantum[1:(end-1)]) < 1.1
    @test minimum(per_quantum[2:end] ./ per_quantum[1:(end-1)]) > 0.9
    @test all(total(resolved) .> abs.(net(resolved)))
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
    @test getproperty.(heated.diagnostics, :total_rate) ≈
          getproperty.(result.solver.diagnostics, :total_rate)
    in_ms = cooling_metrics(s.model, s.mc, method; time_unit=u"ms")
    @test in_ms.τ_c ≈ result.τ_c rtol = 1e-8
    @test unit(in_ms.τ_c) == u"ms"
    @test getproperty.(in_ms.solver.diagnostics, :diffusion) ≈
          getproperty.(result.solver.diagnostics, :diffusion)
    @test all(d -> unit(d.diffusion) == u"ms^-1", in_ms.solver.diagnostics)

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
