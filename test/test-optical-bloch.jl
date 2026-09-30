# Tests for Levels.OpticalBloch: specification validation, the rotating frame
# over the beam graph, the Lindblad model (frame diagonal, couplings, decay and
# dephasing operators) for fine-structure and hyperfine species, and the
# motional layer (Lamb–Dicke factors, sideband Hamiltonian, recoil operators).
# The solver-side checks (needing QuantumToolbox) are in
# test-optical-bloch-solvers.jl.

@testsnippet OpticalBlochSetup begin
    using LinearAlgebra
    using Unitful
    using Levels.OpticalBloch

    # ⁴⁰Ca⁺ with the nonrelativistic Landé g-factors (g_s = 2, g_L = 1: 2, 2/3,
    # 4/5), for closed-form Zeeman bookkeeping in the tests.
    const CA40_NONREL_G = NoHyperfineOneElectronSpecies(;
        mass=ca40.mass,
        energies=ca40.energies,
        einstein_as=ca40.einstein_as,
        lande_g_overrides=Dict(
            convert(NoHyperfineNumberSpec, k) => v for
            (k, v) in ["S_1/2" => 2.0, "P_1/2" => 2 / 3, "D_3/2" => 4 / 5]
        ),
    )
    const LAMBDA_BASIS = StateBasis(["S_1/2", "D_3/2", "P_1/2"])
    const B_LAMBDA = 3.922e-4u"T"

    # Peak intensity of a Gaussian beam of power P and waist w.
    peak_intensity(P, w) = uconvert(u"W/m^2", 2P / (π * w^2))

    # A 397 nm beam at 140° to the field (mixed π/σ± polarisation) plus an 866 nm
    # beam along the field (σ± only), as in a typical Λ dark-resonance scheme.
    function lambda_scheme(
        Δ397,
        Δ866;
        P397=10u"µW",
        P866=50u"µW",
        linewidth=0.0u"µs^-1",
        kwargs...,
    )
        n397, ε397 = beam_vectors(deg2rad(-140.0), 0.0)
        n866, ε866 = beam_vectors(0.0, 0.0)
        LaserScheme(
            LAMBDA_BASIS,
            [
                LaserBeam(
                    "S_1/2" => "P_1/2",
                    Δ397,
                    peak_intensity(P397, 75u"µm"),
                    ε397,
                    n397;
                    linewidth,
                ),
                LaserBeam(
                    "D_3/2" => "P_1/2",
                    Δ866,
                    peak_intensity(P866, 150u"µm"),
                    ε866,
                    n866;
                    linewidth,
                ),
            ];
            static_field=B_LAMBDA,
            kwargs...,
        )
    end

    # Toy hyperfine twin of ⁸⁸Sr⁺ with vanishing hyperfine constants: its
    # at-field model must be the fine-structure model tensored with the nuclear
    # spin, the small nuclear g-factor (which keeps the adiabatic labels well
    # defined by lifting the accidental degeneracies) only adding a nuclear
    # Zeeman shift constant within each m_I sector.
    const SR88_TOY_HF = HyperfineOneElectronSpecies(;
        mass=sr88.mass,
        nuclear_spin=1//2,
        nuclear_g=2e-4,
        energies=sr88.energies,
        lande_g_overrides=sr88.lande_g_overrides,
        hyperfine=Dict(
            convert(NoHyperfineNumberSpec, k) => HyperfineConstants(; a=0.0u"J") for
            k in ["S_1/2", "P_1/2", "D_3/2"]
        ),
        einstein_as=sr88.einstein_as,
    )

    strip_h(M) = ustrip.(u"µs^-1", M)
end

@testitem "LaserBeam, LaserScheme and MotionalMode validation" tags=[:unit, :fast] setup=[
    OpticalBlochSetup,
] begin
    using StaticArrays: SVector

    n, ε = beam_vectors(0.4, 0.3)
    beam = LaserBeam("S_1/2" => "P_1/2", 2π * 10.0u"MHz", 3.0u"W/m^2", 2ε, 5n)
    @test beam.ε isa SVector{3,ComplexF64} && beam.n isa SVector{3,Float64}
    @test beam.ε ≈ ε && beam.n ≈ n   # normalised
    @test iszero(beam.linewidth)
    @test beam.frequency == RelativeFrequency("S_1/2" => "P_1/2", 2π * 10.0u"MHz")
    @test LaserBeam(beam.frequency, 3.0u"mW/cm^2", ε, n; linewidth=2π * 0.1u"MHz").linewidth ==
          2π * 0.1u"MHz"

    @test_throws ArgumentError LaserBeam(beam.frequency, 3.0u"W", ε, n)        # not an intensity
    @test_throws ArgumentError LaserBeam(beam.frequency, 3.0u"W/m^2", ε, ε)      # ε ∥ n
    @test_throws ArgumentError LaserBeam(beam.frequency, 3.0u"W/m^2", 0ε, n)     # zero polarisation
    @test_throws ArgumentError LaserBeam(
        beam.frequency,
        3.0u"W/m^2",
        ε,
        n;
        linewidth=-1.0u"µs^-1",
    )
    @test_throws ArgumentError LaserBeam(
        beam.frequency,
        3.0u"W/m^2",
        ε,
        n;
        linewidth=0.1u"m",
    )
    @test_throws ArgumentError LaserBeam(beam.frequency, 3.0u"W/m^2", ε, [1.0, 2.0])

    basis = StateBasis(["S_1/2", "P_1/2"])
    scheme = LaserScheme(basis, [beam]; static_field=0.5u"mT")
    @test scheme.frame_reference == convert(NoHyperfineNumberSpec, "P_1/2")   # upper of first beam
    @test LaserScheme(basis, [beam]; static_field=0.5u"mT", frame_reference="S_1/2").frame_reference ==
          convert(NoHyperfineNumberSpec, "S_1/2")
    @test LaserScheme(basis, LaserBeam[]; static_field=0.5u"mT").frame_reference ==
          convert(NoHyperfineNumberSpec, "S_1/2")
    @test_throws ArgumentError LaserScheme(basis, [beam]; static_field=0.5u"m")
    @test_throws ArgumentError LaserScheme(
        basis,
        [beam];
        static_field=0.5u"mT",
        frame_reference="D_3/2",
    )
    other = LaserBeam("D_3/2" => "P_1/2", 0.0u"µs^-1", 3.0u"W/m^2", ε, n)
    @test_throws ArgumentError LaserScheme(basis, [beam, other]; static_field=0.5u"mT")
    # Hyperfine reference levels on a fine-structure basis (and vice versa).
    hf_beam = LaserBeam("S_1/2 F=4" => "P_1/2 F=4", 0.0u"µs^-1", 3.0u"W/m^2", ε, n)
    @test_throws ArgumentError LaserScheme(basis, [hf_beam]; static_field=0.5u"mT")
    @test_throws ArgumentError LaserScheme(
        StateBasis(ca43, "S_1/2", "P_1/2"),
        [beam];
        static_field=0.5u"mT",
    )

    mode = MotionalMode(2π * 1.2u"MHz", [1.0, 0.0, 1.0])
    @test mode.direction ≈ [1, 0, 1] / sqrt(2)
    @test_throws ArgumentError MotionalMode(1.2u"MHz" * 0, [0.0, 0.0, 1.0])
    @test_throws ArgumentError MotionalMode(1.2u"m", [0.0, 0.0, 1.0])
    @test_throws ArgumentError MotionalMode(2π * 1.2u"MHz", [0.0, 0.0, 0.0])
end

@testitem "Rotating frame over the beam graph" tags=[:unit, :fast] setup=[
    OpticalBlochSetup,
] begin
    basis = StateBasis(["S_1/2", "P_1/2", "D_3/2"])
    n, ε = beam_vectors(π / 2, π / 4)
    Δ1 = 2π * 20.0u"MHz"
    Δ2 = 2π * -35.0u"MHz"
    I0 = 1.0u"W/m^2"
    b422 = LaserBeam("S_1/2" => "P_1/2", Δ1, I0, ε, n)
    b1092 = LaserBeam("D_3/2" => "P_1/2", Δ2, I0, ε, n)
    S, P, D = (convert(NoHyperfineNumberSpec, l) for l in ("S_1/2", "P_1/2", "D_3/2"))

    # Reference P (default): the lower levels sit at minus their beams' detunings.
    frame =
        rotating_frame(sr88, LaserScheme(basis, [b422, b1092]; static_field=0.5u"mT"))
    @test frame.reference == P
    @test frame.frame_offsets[P] == 0.0u"µs^-1"
    @test frame.frame_offsets[S] ≈ -Δ1
    @test frame.frame_offsets[D] ≈ -Δ2
    @test frame.tree_beams == [1, 2]
    @test isempty(frame.beats)
    @test frame.traversals[P] == [0, 0] &&
          frame.traversals[S] == [-1, 0] &&
          frame.traversals[D] == [0, -1]

    # Reference S: same frame up to a common shift; the path to D runs up the
    # 422 nm beam and down the 1092 nm one.
    frame_s = rotating_frame(
        sr88,
        LaserScheme(
            basis,
            [b422, b1092];
            static_field=0.5u"mT",
            frame_reference="S_1/2",
        ),
    )
    @test frame_s.frame_offsets[S] == 0.0u"µs^-1"
    @test frame_s.frame_offsets[P] ≈ Δ1
    @test frame_s.frame_offsets[D] ≈ Δ1 - Δ2
    @test frame_s.traversals[D] == [1, -1]
    for l in (S, P, D)
        @test frame_s.frame_offsets[l] - frame.frame_offsets[l] ≈ Δ1
    end

    # A second beam on an already connected pair beats at the detuning
    # difference; an identical-frequency one does not (and adds coherently).
    Δ3 = 2π * 26.0u"MHz"
    b422b = LaserBeam("S_1/2" => "P_1/2", Δ3, I0, ε, n)
    frame_b = rotating_frame(
        sr88,
        LaserScheme(basis, [b422, b1092, b422b]; static_field=0.5u"mT"),
    )
    @test frame_b.tree_beams == [1, 2]
    @test length(frame_b.beats) == 1 && frame_b.beats[1][1] == 3
    @test frame_b.beats[1][2] ≈ Δ3 - Δ1
    model_b = lindblad_model(
        sr88,
        LaserScheme(basis, [b422, b1092, b422b]; static_field=0.5u"mT"),
    )
    @test length(model_b.harmonics) == 1
    @test model_b.harmonics[1][1] ≈ Δ3 - Δ1
    model_1 =
        lindblad_model(sr88, LaserScheme(basis, [b422, b1092]; static_field=0.5u"mT"))
    model_2 = lindblad_model(
        sr88,
        LaserScheme(basis, [b422, b1092, b422]; static_field=0.5u"mT"),
    )
    @test isempty(model_2.harmonics)
    s_range, p_range = staterange(basis, "S_1/2"), staterange(basis, "P_1/2")
    @test model_2.hamiltonian[p_range, s_range] ≈
          2 .* model_1.hamiltonian[p_range, s_range]
    @test model_2.hamiltonian[s_range, s_range] ≈ model_1.hamiltonian[s_range, s_range]

    # A loop of beams: the fourth beam closes S–P₁/₂–D₃/₂–P₃/₂–S and beats at the
    # frequency mismatch around the loop.
    basis4 = StateBasis(["S_1/2", "P_1/2", "D_3/2", "P_3/2"])
    Δ4, Δ5 = 2π * 12.0u"MHz", 2π * -7.0u"MHz"
    b408 = LaserBeam("S_1/2" => "P_3/2", Δ4, I0, ε, n)
    b1004 = LaserBeam("D_3/2" => "P_3/2", Δ5, I0, ε, n)
    frame_4 = rotating_frame(
        sr88,
        LaserScheme(basis4, [b422, b1092, b408, b1004]; static_field=0.5u"mT"),
    )
    @test frame_4.tree_beams == [1, 2, 3]
    @test length(frame_4.beats) == 1 && frame_4.beats[1][1] == 4
    @test frame_4.beats[1][2] ≈ Δ5 + Δ1 - Δ4 - Δ2

    # A level no beam reaches keeps its centroid as frame frequency.
    basis5 = StateBasis(["S_1/2", "P_1/2", "D_3/2", "D_5/2"])
    frame_5 =
        rotating_frame(sr88, LaserScheme(basis5, [b422, b1092]; static_field=0.5u"mT"))
    @test frame_5.frame_offsets[convert(NoHyperfineNumberSpec, "D_5/2")] == 0.0u"µs^-1"

    # Hyperfine reference levels: the frame offset across a beam includes the
    # zero-field hyperfine shifts of its reference F levels.
    hb = StateBasis(ca43, "S_1/2", "P_1/2")
    beam = LaserBeam("S_1/2 F=4" => "P_1/2 F=3", Δ1, I0, ε, n)
    frame_hf = rotating_frame(ca43, LaserScheme(hb, [beam]; static_field=0.5u"mT"))
    @test frame_hf.frame_offsets[S] ≈ -(
        Δ1 + Levels.hyperfine_shift(ca43, "P_1/2 F=3") -
        Levels.hyperfine_shift(ca43, "S_1/2 F=4")
    )
end

@testitem "Lindblad model of the Λ system" tags=[:unit, :fast] setup=[OpticalBlochSetup] begin
    Δ397, Δ866 = 2π * 300.0u"MHz", 2π * 290.0u"MHz"
    scheme = lambda_scheme(Δ397, Δ866)
    model = lindblad_model(CA40_NONREL_G, scheme)
    basis = LAMBDA_BASIS
    u = uconvert(u"µs^-1", Levels.BOHR_MAGNETON * B_LAMBDA / u"ħ")
    H = model.hamiltonian
    @test ishermitian(strip_h(H))

    # Frame diagonal with the P level as reference: Zeeman shifts g m u, plus the
    # detuning of the beam connecting to the reference for the lower levels.
    g = Dict("S_1/2" => 2.0, "D_3/2" => 4 / 5, "P_1/2" => 2 / 3)
    detuning = Dict("S_1/2" => Δ397, "D_3/2" => Δ866, "P_1/2" => 0.0u"µs^-1")
    for (i, state) in enumerate(basis)
        key =
            only(k for k in keys(g) if convert(NoHyperfineNumberSpec, k) == state.level)
        @test real(H[i, i]) ≈ g[key] * state.m * u + detuning[key] rtol = 1e-12
    end
    # With the S level as reference the whole diagonal shifts by −Δ397.
    model_s = lindblad_model(
        CA40_NONREL_G,
        lambda_scheme(Δ397, Δ866; frame_reference="S_1/2"),
    )
    @test real.(diag(model_s.hamiltonian)) ≈ real.(diag(H)) .- Δ397
    @test model_s.hamiltonian - Diagonal(model_s.hamiltonian) ≈ H - Diagonal(H)

    # Couplings: the (C + C†)/2 of each beam's coupling_matrix.
    C = sum(
        coupling_matrix(
            CA40_NONREL_G,
            basis,
            OpticalBloch.beam_levels(beam)[1] => OpticalBloch.beam_levels(beam)[2],
            beam.intensity,
            beam.ε,
            beam.n,
        ) for beam in scheme.beams
    )
    @test H - Diagonal(H) ≈ (C + C') / 2
    # The 866 nm beam along ẑ has no π component; the 397 nm beam drives all three.
    p, s, d = staterange(basis, "P_1/2"),
    staterange(basis, "S_1/2"),
    staterange(basis, "D_3/2")
    @test all(iszero, [H[k, i] for i in d, k in p if basis[k].m == basis[i].m])
    @test all(!iszero, H[p, s])

    # Decay operators: one per emitted component, ∑ L†L = Γ P_upper, Clebsch–
    # Gordan branching fractions of the components.
    ops, labels = decay_operators(model)
    @test all(l isa DecayLabel && isnothing(l.component) for l in labels)
    A_ps = einstein_a(CA40_NONREL_G, "S_1/2", "P_1/2")
    A_pd = einstein_a(CA40_NONREL_G, "D_3/2", "P_1/2")
    A_ds = einstein_a(CA40_NONREL_G, "S_1/2", "D_3/2")
    total = sum(L' * L for L in ops)
    expected = zeros(ComplexF64, 8, 8) .* u"µs^-1"
    for i in p
        expected[i, i] = A_ps + A_pd
    end
    for i in d
        expected[i, i] = A_ds
    end
    @test strip_h(total) ≈ strip_h(expected) atol = 1e-12
    fraction(L, i, k) = ustrip(u"µs^-1", abs2(L[i, k]))
    P_lo(m) = stateindex(basis, "P_1/2", m)
    S_lo(m) = stateindex(basis, "S_1/2", m)
    D_lo(m) = stateindex(basis, "D_3/2", m)
    idx(lo, q) = findfirst(
        l ->
            l.lower == convert(NoHyperfineNumberSpec, lo) &&
            l.upper == convert(NoHyperfineNumberSpec, "P_1/2") &&
            l.q == q,
        labels,
    )
    @test fraction(ops[idx("S_1/2", 0)], S_lo(-1//2), P_lo(-1//2)) /
          ustrip(u"µs^-1", A_ps) ≈ 1 / 3
    @test fraction(ops[idx("S_1/2", 1)], S_lo(-1//2), P_lo(1//2)) /
          ustrip(u"µs^-1", A_ps) ≈ 2 / 3
    @test fraction(ops[idx("D_3/2", 1)], D_lo(-3//2), P_lo(-1//2)) /
          ustrip(u"µs^-1", A_pd) ≈ 1 / 2
    @test fraction(ops[idx("D_3/2", 1)], D_lo(-1//2), P_lo(1//2)) /
          ustrip(u"µs^-1", A_pd) ≈ 1 / 6
    @test fraction(ops[idx("D_3/2", 0)], D_lo(-1//2), P_lo(-1//2)) /
          ustrip(u"µs^-1", A_pd) ≈ 1 / 3
    # The relative sign of the two components of the π operator is the CG sign.
    Lπ = ops[idx("S_1/2", 0)]
    @test real(Lπ[S_lo(-1//2), P_lo(-1//2)] / Lπ[S_lo(1//2), P_lo(1//2)]) ≈ -1
    @test length(dephasing_operators(model)[1]) == 0

    # Component-resolved decay: one operator per state pair, the same ∑ L†L.
    resolved = lindblad_model(CA40_NONREL_G, scheme; decay=:resolved)
    rops, rlabels = decay_operators(resolved)
    @test length(rops) == 4 + 6 + 8   # P→S (2 × 2), P→D (2 × 3), D→S (E2: all 8 pairs)
    @test all(l.component isa Pair for l in rlabels)
    @test all(count(!iszero, L) == 1 for L in rops)
    @test strip_h(sum(L' * L for L in rops)) ≈ strip_h(total) atol = 1e-12

    # Open decay channels: error naming the missing level, or dropped on request.
    open_basis = StateBasis(["S_1/2", "P_1/2"])
    n, ε = beam_vectors(π / 2, π / 4)
    open_scheme = LaserScheme(
        open_basis,
        [LaserBeam("S_1/2" => "P_1/2", Δ397, 1.0u"W/m^2", ε, n)];
        static_field=B_LAMBDA,
    )
    err = try
        lindblad_model(CA40_NONREL_G, open_scheme)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("D_3/2", err.msg) ||
          occursin("(2, 3//2)", err.msg)
    dropped = lindblad_model(CA40_NONREL_G, open_scheme; open_channels=:drop)
    @test length(decay_operators(dropped)[1]) == 3
    tot_dropped = sum(L' * L for L in decay_operators(dropped)[1])
    @test all(
        ustrip(u"µs^-1", real(tot_dropped[i, i])) ≈ ustrip(u"µs^-1", A_ps) for i in 3:4
    )
    # A partial manifold drops the decay channels into its missing states.
    partial = StateBasis([StateSpec("S_1/2", -1//2), StateSpec("P_1/2", 1//2)])
    n_z, ε_σ = beam_vectors(0.0, π / 4, π / 2)
    partial_model = lindblad_model(
        sr88,
        LaserScheme(
            partial,
            [LaserBeam("S_1/2" => "P_1/2", 0.0u"µs^-1", 1.0u"W/m^2", ε_σ, n_z)];
            static_field=0.0u"mT",
        );
        open_channels=:drop,
    )
    @test length(decay_operators(partial_model)[1]) == 1
    @test ustrip(
        u"µs^-1",
        real(sum(L' * L for L in partial_model.jump_operators)[2, 2]),
    ) ≈ 2 / 3 * ustrip(u"µs^-1", einstein_a(sr88, "S_1/2", "P_1/2"))
    @test_throws ArgumentError lindblad_model(CA40_NONREL_G, scheme; decay=:other)
    @test_throws ArgumentError lindblad_model(
        CA40_NONREL_G,
        scheme;
        open_channels=:ignore,
    )

    # Level projectors and the population bookkeeping helper.
    @test diag(level_projector(model, "P_1/2")) == [i in p for i in 1:8]
    @test sum(level_projector(model, l) for l in ("S_1/2", "D_3/2", "P_1/2")) == I(8)
end

@testitem "Laser dephasing operators" tags=[:unit, :fast] setup=[OpticalBlochSetup] begin
    γ = 2π * 0.1u"MHz"
    Δ397, Δ866 = 2π * 300.0u"MHz", 2π * 290.0u"MHz"
    model_p = lindblad_model(CA40_NONREL_G, lambda_scheme(Δ397, Δ866; linewidth=γ))
    model_s = lindblad_model(
        CA40_NONREL_G,
        lambda_scheme(Δ397, Δ866; linewidth=γ, frame_reference="S_1/2"),
    )
    ops_p, labels_p = dephasing_operators(model_p)
    ops_s, _ = dephasing_operators(model_s)
    @test [l.beam for l in labels_p] == [1, 2]
    basis = LAMBDA_BASIS
    # Reference P: each beam dephases the lower level it connects to P.
    @test ustrip.(u"µs^-1", ops_p[1]' * ops_p[1]) ≈
          ustrip(u"µs^-1", γ) .* level_projector(model_p, "S_1/2")
    @test ustrip.(u"µs^-1", ops_p[2]' * ops_p[2]) ≈
          ustrip(u"µs^-1", γ) .* level_projector(model_p, "D_3/2")
    # Reference S: the 397 nm beam acts on P and D, the 866 nm one on D alone —
    # yet every level-pair coherence decays at the same rate as before,
    # ∑_l γ_l (n_{l,i} − n_{l,j})² / 2.
    @test ustrip.(u"µs^-1", ops_s[1]' * ops_s[1]) ≈
          ustrip(u"µs^-1", γ) .*
          (level_projector(model_s, "P_1/2") + level_projector(model_s, "D_3/2"))
    rate(ops, i, j) = sum(abs2(ustrip(u"µs^-0.5", L[i, i] - L[j, j])) for L in ops) / 2
    for i in 1:8, j in 1:8
        @test rate(ops_s, i, j) ≈ rate(ops_p, i, j) atol = 1e-12
    end
    s, d, p = (first(staterange(basis, l)) for l in ("S_1/2", "D_3/2", "P_1/2"))
    @test rate(ops_p, s, p) ≈ ustrip(u"µs^-1", γ) / 2
    @test rate(ops_p, s, d) ≈ ustrip(u"µs^-1", γ)          # Raman coherence: both linewidths
    @test rate(ops_p, s, s + 1) == 0                       # within a level: nothing

    # A beat-note beam cannot carry a linewidth.
    n, ε = beam_vectors(deg2rad(-140.0), 0.0)
    extra =
        LaserBeam("S_1/2" => "P_1/2", 2π * 310.0u"MHz", 1.0u"W/m^2", ε, n; linewidth=γ)
    scheme = lambda_scheme(Δ397, Δ866)
    @test_throws ArgumentError lindblad_model(
        CA40_NONREL_G,
        LaserScheme(LAMBDA_BASIS, [scheme.beams; extra]; static_field=B_LAMBDA),
    )
end

@testitem "Hyperfine Lindblad models" tags=[:unit] setup=[OpticalBlochSetup] begin
    B = 0.5u"mT"
    n, ε = beam_vectors(π / 2, π / 4)
    I0 = 2.0u"W/m^2"
    Δ = 2π * 5.0u"MHz"

    # Zero hyperfine constants: the hyperfine model is the fine-structure one
    # tensored with the nuclear spin — the Hamiltonian spectrum doubles, each
    # copy shifted by the nuclear Zeeman energy g_I m_I μ_B B / ħ — with the
    # same decay-rate sum rule.
    hf_basis = StateBasis(SR88_TOY_HF, "S_1/2", "P_1/2", "D_3/2")
    fs_basis = StateBasis(["S_1/2", "P_1/2", "D_3/2"])
    hf_scheme = LaserScheme(
        hf_basis,
        [
            LaserBeam("S_1/2 F=1" => "P_1/2 F=1", Δ, I0, ε, n),
            LaserBeam("D_3/2 F=2" => "P_1/2 F=1", -Δ, I0, ε, n),
        ];
        static_field=B,
    )
    fs_scheme = LaserScheme(
        fs_basis,
        [
            LaserBeam("S_1/2" => "P_1/2", Δ, I0, ε, n),
            LaserBeam("D_3/2" => "P_1/2", -Δ, I0, ε, n),
        ];
        static_field=B,
    )
    hf = lindblad_model(SR88_TOY_HF, hf_scheme)
    fs = lindblad_model(sr88, fs_scheme)
    @test length(hf_basis) == 2 * length(fs_basis)
    @test ishermitian(strip_h(hf.hamiltonian))
    e_hf = sort(eigvals(Hermitian(strip_h(hf.hamiltonian))))
    e_fs = eigvals(Hermitian(strip_h(fs.hamiltonian)))
    u_I = SR88_TOY_HF.nuclear_g * ustrip(u"µs^-1", Levels.BOHR_MAGNETON * B / u"ħ")
    @test e_hf ≈ sort([e + m_I * u_I for e in e_fs for m_I in (-1 / 2, 1 / 2)]) atol = 1e-8
    total_hf = sum(L' * L for L in hf.jump_operators)
    for (i, state) in enumerate(hf_basis)
        fs_level = fine_structure(state.level)
        expected = 1 / something(lifetime(sr88, fs_level), Inf * u"s")
        @test ustrip(u"µs^-1", real(total_hf[i, i])) ≈ ustrip(u"µs^-1", expected) atol =
            1e-12
    end
    @test all(
        abs(ustrip(u"µs^-1", total_hf[i, k])) < 1e-12 for
        i in 1:length(hf_basis), k in 1:length(hf_basis) if i != k
    )
    # Each fine-structure decay component appears with its rate split over the
    # nuclear-spin copies.
    @test length(decay_operators(hf)[1]) == length(decay_operators(fs)[1])

    # ⁴³Ca⁺ S₁/₂ ⊕ P₁/₂ at 0.5 mT with the D₃/₂ branch dropped: Hermitian,
    # trace-preserving on the kept channel, and undefined at zero field.
    basis = StateBasis(ca43, "S_1/2", "P_1/2")
    beam = LaserBeam("S_1/2 F=4" => "P_1/2 F=4", Δ, I0, ε, n)
    scheme = LaserScheme(basis, [beam]; static_field=B)
    @test_throws ArgumentError lindblad_model(ca43, scheme)
    model = lindblad_model(ca43, scheme; open_channels=:drop)
    @test ishermitian(strip_h(model.hamiltonian))
    total = sum(L' * L for L in model.jump_operators)
    A = ustrip(u"µs^-1", einstein_a(ca43, "S_1/2", "P_1/2"))
    for i in staterange(basis, "P_1/2")
        @test ustrip(u"µs^-1", real(total[i, i])) ≈ A rtol = 1e-10
    end
    for i in staterange(basis, "S_1/2")
        @test abs(ustrip(u"µs^-1", total[i, i])) < 1e-12
    end
    @test_throws ArgumentError lindblad_model(
        ca43,
        LaserScheme(basis, [beam]; static_field=0.0u"mT");
        open_channels=:drop,
    )
end

@testitem "Motional coupling" tags=[:unit, :fast] setup=[OpticalBlochSetup] begin
    Δ397, Δ866 = 2π * 300.0u"MHz", 2π * 290.0u"MHz"
    scheme = lambda_scheme(Δ397, Δ866)
    model = lindblad_model(CA40_NONREL_G, scheme)
    ω_m = 2π * 1.1u"MHz"
    mode_z = MotionalMode(ω_m, [0, 0, 1.0])
    mc = motional_coupling(CA40_NONREL_G, scheme, model, mode_z)

    # Lamb–Dicke factors: the projected value is k·e_m x₀ (signed), the
    # unprojected one |k| x₀.
    x0 = sqrt(u"ħ" / (2 * CA40_NONREL_G.mass * ω_m))
    k397 =
        2π / (2π * u"c" / Levels.transition_frequency(CA40_NONREL_G, "S_1/2", "P_1/2"))
    η397 = uconvert(NoUnits, k397 * x0)
    @test mc.projected_lamb_dicke[1] ≈ η397 * cos(deg2rad(140.0)) rtol = 1e-6   # thesis-like −0.130
    @test mc.projected_lamb_dicke[1] ≈ -0.130 atol = 5e-3
    @test mc.projected_lamb_dicke[2] ≈ 0.0778 atol = 5e-4
    @test lamb_dicke(CA40_NONREL_G, mode_z, scheme.beams[1]; projected=false) ≈ η397 rtol =
        1e-6
    @test lamb_dicke(CA40_NONREL_G, mode_z, scheme.beams[1]) ≈
          mc.projected_lamb_dicke[1]
    # Sign flips with the beam direction, magnitude scales with cos of the
    # beam–mode angle, and a transverse mode sees nothing.
    flipped = LaserBeam(
        scheme.beams[1].frequency,
        scheme.beams[1].intensity,
        scheme.beams[1].ε,
        -scheme.beams[1].n,
    )
    @test lamb_dicke(CA40_NONREL_G, mode_z, flipped) ≈ -mc.projected_lamb_dicke[1]
    @test lamb_dicke(CA40_NONREL_G, MotionalMode(ω_m, [0, 1.0, 0]), scheme.beams[1]) ≈ 0 atol =
        1e-15
    @test lamb_dicke(CA40_NONREL_G, MotionalMode(ω_m, [1.0, 0, 0]), scheme.beams[1]) ≈
          η397 * sin(deg2rad(-140.0)) rtol = 1e-6
    # Counter-propagating beams: the two-photon factor is the difference.
    @test mc.projected_lamb_dicke[1] - mc.projected_lamb_dicke[2] ≈
          -(abs(mc.projected_lamb_dicke[1]) + mc.projected_lamb_dicke[2])

    # Sideband Hamiltonian: ∑ η_l i (C_l − C_l†)/2, Hermitian, purely off-diagonal.
    C = [
        coupling_matrix(
            CA40_NONREL_G,
            LAMBDA_BASIS,
            OpticalBloch.beam_levels(beam)[1] => OpticalBloch.beam_levels(beam)[2],
            beam.intensity,
            beam.ε,
            beam.n,
        ) for beam in scheme.beams
    ]
    expected =
        sum(η * im .* (Cl .- Cl') ./ 2 for (η, Cl) in zip(mc.projected_lamb_dicke, C))
    @test strip_h(mc.sideband_hamiltonian) ≈ strip_h(expected)
    @test ishermitian(strip_h(mc.sideband_hamiltonian))
    @test all(iszero, diag(mc.sideband_hamiltonian))
    @test isempty(mc.sideband_harmonics)

    # Recoil: one operator per decay operator, √α η₀ L with the unprojected η₀
    # of the transition and the emission-pattern moment of its component.
    ops, labels = decay_operators(model)
    @test mc.recoil_indices == findall(l -> l isa DecayLabel, model.jump_labels)
    @test length(mc.recoil_operators) == length(ops)
    for (k, (L, label)) in enumerate(zip(ops, labels))
        k_photon =
            Levels.transition_frequency(CA40_NONREL_G, label.lower, label.upper) / u"c"
        @test mc.recoil_lamb_dicke[k] ≈ uconvert(NoUnits, k_photon * x0) rtol = 1e-9
        rank = Levels.multipole_rank(label.lower, label.upper)
        @test mc.recoil_moments[k] == recoil_moment(rank, label.q, mode_z.direction)
        @test mc.recoil_operators[k] ≈
              (sqrt(mc.recoil_moments[k]) * mc.recoil_lamb_dicke[k]) .* L
    end
    # Emission-pattern moments: 1/5 (π) and 2/5 (σ) along ẑ, 2/5 and 3/10
    # transverse, 1/3 on average in any direction; E2 patterns likewise
    # average to 1/3.
    @test recoil_moment(1, 0, [0, 0, 1.0]) ≈ 1 / 5
    @test recoil_moment(1, 1, [0, 0, 1.0]) ≈ 2 / 5
    @test recoil_moment(1, -1, [0, 0, 1.0]) ≈ 2 / 5
    @test recoil_moment(1, 0, [1.0, 0, 0]) ≈ 2 / 5
    @test recoil_moment(1, 1, [0, 1.0, 0]) ≈ 3 / 10
    for direction in ([0, 0, 1.0], [1.0, 0, 0], [0.3, -0.4, 0.5])
        @test sum(recoil_moment(1, q, direction) for q in -1:1) / 3 ≈ 1 / 3
        @test sum(recoil_moment(2, q, direction) for q in -2:2) / 5 ≈ 1 / 3
    end
    @test_throws ArgumentError recoil_moment(1, 2, [0, 0, 1.0])
    @test_throws ArgumentError recoil_moment(3, 0, [0, 0, 1.0])
    # A constant moment replaces the pattern average.
    mc_const =
        motional_coupling(CA40_NONREL_G, scheme, model, mode_z; recoil_moment=2 / 5)
    @test all(==(2 / 5), mc_const.recoil_moments)
    @test_throws ArgumentError motional_coupling(
        CA40_NONREL_G,
        scheme,
        model,
        mode_z;
        recoil_moment=:isotropic,
    )

    # Several modes at once.
    modes = [mode_z, MotionalMode(2π * 2.0u"MHz", [1.0, 0, 1.0])]
    mcs = motional_coupling(CA40_NONREL_G, scheme, model, modes)
    @test length(mcs) == 2
    @test mcs[1].sideband_hamiltonian == mc.sideband_hamiltonian
    @test mcs[2].projected_lamb_dicke ≈
          [lamb_dicke(CA40_NONREL_G, modes[2], b) for b in scheme.beams]
end
