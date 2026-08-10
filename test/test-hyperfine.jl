@testitem "Hyperfine levels and species basis" tags=[:unit, :fast] begin
    using Unitful

    @test hyperfine_levels(ca43, "S_1/2") ==
          [HyperfineNumberSpec(0, 1//2, f) for f in 3:4]
    @test hyperfine_levels(ca43, "D_5/2") ==
          [HyperfineNumberSpec(2, 5//2, f) for f in 1:6]

    # Manifold expansion: ascending F, each ascending m_F.
    basis = StateBasis(ca43, "S_1/2")
    @test length(basis) == 16
    @test basis[1] == StateSpec("S_1/2 F=3", -3)
    @test basis[7] == StateSpec("S_1/2 F=3", 3)
    @test basis[8] == StateSpec("S_1/2 F=4", -4)
    @test basis[end] == StateSpec("S_1/2 F=4", 4)

    both = StateBasis(ca43, "S_1/2", "D_5/2")
    @test length(both) == 16 + 48
    @test staterange(both, "S_1/2") == 1:16
    @test staterange(both, "D_5/2") == 17:64
    @test staterange(both, "D_5/2 F=6") == (64-12):64

    # Single-F entries stay as given; invalid F values are rejected.
    @test length(StateBasis(ca43, "S_1/2 F=4")) == 9
    @test_throws ArgumentError StateBasis(ca43, "S_1/2 F=5")
    @test_throws ArgumentError StateBasis(ca43, "S_1/2 F=7/2")

    # Measured g_J overrides apply; other levels fall back to the LS formula.
    @test lande_g(ca43, "S_1/2") == 2.00225664
    @test lande_g(ca43, "D_5/2") == 1.2003340
    @test lande_g(ca43, "P_1/2") ≈ 2//3 rtol = 2e-3

    # Low-field g_F from the Landé projection formula, including the nuclear
    # term: g_F = g_J⋅[F(F+1)+J(J+1)−I(I+1)]/(2F(F+1)) + g_I⋅[…].
    g_f4 = lande_g(ca43, "S_1/2 F=4")
    g_f3 = lande_g(ca43, "S_1/2 F=3")
    g_i = ca43.nuclear_g
    @test g_f4 ≈ 2.00225664 * (1 / 8) + g_i * (7 / 8) rtol = 1e-12
    @test g_f3 ≈ 2.00225664 * (-1 / 8) + g_i * (9 / 8) rtol = 1e-12
    @test_throws ArgumentError lande_g(ca43, "S_1/2 F=5")
    @test_throws ArgumentError lande_g(sr88, "S_1/2 F=4")
end

@testitem "Zero-field hyperfine structure" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful
    using Levels: manifold_hamiltonian

    # The measured ⁴³Ca⁺ ground-state splitting [Arbes1994] is recovered from
    # the Casimir shifts (A(I + 1/2), F = 3 above F = 4 as A < 0).
    splitting = hyperfine_shift(ca43, "S_1/2 F=3") - hyperfine_shift(ca43, "S_1/2 F=4")
    @test splitting ≈ 2π * 3225.60828640u"MHz" rtol = 1e-12

    # The (2F+1)-weighted mean of the Casimir shifts vanishes: both the A and
    # the B term are traceless over the manifold, so the centroid is preserved.
    for level in ("S_1/2", "D_5/2", "D_3/2", "P_3/2")
        levels = hyperfine_levels(ca43, level)
        mean = sum(Int(2l.f + 1) * hyperfine_shift(ca43, l) for l in levels)
        @test abs(mean) < 2π * 1e-3u"Hz"
    end

    # The zero-field manifold Hamiltonian is diagonal with the Casimir shifts —
    # D_5/2 exercises both the A and the B term.
    basis = StateBasis(ca43, "D_5/2")
    h = manifold_hamiltonian(ca43, "D_5/2", 0.0u"mT")
    expected = Diagonal([hyperfine_shift(ca43, s.level) for s in basis])
    @test all(abs.(h .- expected) .< 2π * 1e-6u"Hz")

    # Casimir shift needs known constants and a genuine F level.
    @test_throws ArgumentError hyperfine_shift(ca43, "F_5/2 F=3")
    @test_throws ArgumentError hyperfine_shift(ca43, "D_5/2")

    # Hyperfine-resolved transition frequency = centroid + Casimir offsets.
    @test Levels.transition_frequency(ca43, "S_1/2 F=4", "D_5/2 F=6") ≈
          Levels.transition_frequency(ca43, "S_1/2", "D_5/2") -
          hyperfine_shift(ca43, "S_1/2 F=4") + hyperfine_shift(ca43, "D_5/2 F=6")
end

@testitem "S1/2 Breit–Rabi" tags=[:unit, :fast] begin
    using Unitful

    # Closed-form Breit–Rabi energies for the J = 1/2 ground manifold (in the
    # convention H_Z = μ_B B (g_J m_J + g_I m_I); cf. e.g. Woodgate,
    # "Elementary Atomic Structure", §9.4), with the stretched states written
    # out exactly (they are product states, so no square root).
    i_nuc = 7 / 2
    a_hz = ustrip(u"Hz", ca43.hyperfine[NoHyperfineNumberSpec(0, 1//2)].a / u"h")
    g_j = lande_g(ca43, "S_1/2")
    g_i = ca43.nuclear_g
    μb_h = ustrip(u"Hz/T", Levels.BOHR_MAGNETON / u"h")

    function breit_rabi(b, f, m_f)
        Δe = a_hz * (i_nuc + 1 / 2)
        x = (g_j - g_i) * μb_h * b / Δe
        e0 = -Δe / (2 * (2 * i_nuc + 1)) + g_i * μb_h * m_f * b
        if f == i_nuc + 1 / 2 && abs(m_f) == i_nuc + 1 / 2
            return e0 + (Δe / 2) * (1 + sign(m_f) * x)
        end
        s = f == i_nuc + 1 / 2 ? 1 : -1
        e0 + s * (Δe / 2) * sqrt(1 + 4 * m_f * x / (2 * i_nuc + 1) + x^2)
    end

    for b_mt in (0.01, 0.1, 0.4965, 1.0)
        m = hyperfine_manifold(ca43, "S_1/2", b_mt * u"mT")
        for state in m.basis
            e_pkg = ustrip(u"s^-1", state_energy(m, state)) / 2π
            e_ref = breit_rabi(b_mt * 1e-3, Float64(state.level.f), Float64(state.m))
            @test e_pkg ≈ e_ref rtol = 1e-12 atol = 1e-3
        end
    end
end

@testitem "Adiabatic labels and eigenstructure" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful
    using Levels: manifold_hamiltonian

    m = hyperfine_manifold(ca43, "D_5/2", 0.496u"mT")

    # m_F is exact: the eigenvector labelled (F, m_F) lives entirely in the
    # m_F block of the coupled basis.
    fz = [Float64(s.m) for s in m.basis]
    for (i, state) in enumerate(m.basis)
        v = m.states[:, i]
        @test sum(abs2.(v) .* fz) ≈ Float64(state.m) atol = 1e-12
    end

    # Eigenvectors are orthonormal, and diagonalise the manifold Hamiltonian
    # with the stored (label-aligned) energies.
    @test m.states' * m.states ≈ I atol = 1e-12
    h = ustrip.(u"µs^-1", manifold_hamiltonian(ca43, "D_5/2", 0.496u"mT"))
    @test m.states' * h * m.states ≈ diagm(ustrip.(u"µs^-1", m.energies)) atol = 1e-10

    # Labels are undefined at zero field.
    @test_throws ArgumentError hyperfine_manifold(ca43, "D_5/2", 0.0u"mT")

    # At low field, the exact shift approaches g_F m_F μ_B B / ħ. The residual
    # second-order term is ∼ g_J μ_B B / ΔE_hfs relative to the linear shift
    # (measured to scale as such): ~2e-6 for S_1/2 at 1e-4 mT, and ~1000×
    # larger for D_5/2, whose F intervals are MHz-scale — hence its much
    # smaller test field.
    for (state, b) in (
        (StateSpec("S_1/2 F=4", 3), 1e-4u"mT"),
        (StateSpec("S_1/2 F=3", -2), 1e-4u"mT"),
        (StateSpec("D_5/2 F=2", 1), 3e-7u"mT"),
    )
        expected = uconvert(
            u"µs^-1",
            lande_g(ca43, state.level) * state.m * Levels.BOHR_MAGNETON * b / u"ħ",
        )
        @test zeeman_shift(ca43, state, b) ≈ expected rtol = 1e-5
    end
    @test iszero(zeeman_shift(ca43, StateSpec("S_1/2 F=4", 4), 0.0u"mT"))
    @test_throws ArgumentError zeeman_shift(ca43, StateSpec("S_1/2", 1//2), 0.1u"mT")
end

@testitem "Hyperfine Zeeman Hamiltonian" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful
    using Levels: jx_matrix, jy_matrix, jz_matrix, coupling_transform

    # Independent product-basis construction (nuclear factor slow, electronic
    # fast), conjugated with the Clebsch–Gordan unitary, must reproduce the
    # coupled-basis matrix elements.
    b = [0.2, -0.1, 0.3]u"mT"
    basis = StateBasis(ca43, "S_1/2")
    h = zeeman_hamiltonian(ca43, basis, b)
    @test h ≈ h'

    i_nuc = 7//2
    eye_i = Matrix{Float64}(I, 8, 8)
    eye_j = Matrix{Float64}(I, 2, 2)
    product =
        (Levels.BOHR_MAGNETON / u"ħ") * (
            lande_g(ca43, "S_1/2") * (
                b[1] .* kron(eye_i, jx_matrix(1//2)) .+
                b[2] .* kron(eye_i, jy_matrix(1//2)) .+
                b[3] .* kron(eye_i, jz_matrix(1//2))
            ) .+
            ca43.nuclear_g * (
                b[1] .* kron(jx_matrix(i_nuc), eye_j) .+
                b[2] .* kron(jy_matrix(i_nuc), eye_j) .+
                b[3] .* kron(jz_matrix(i_nuc), eye_j)
            )
        )
    u = coupling_transform(ca43, "S_1/2")
    @test maximum(abs.(h .- u' * product * u)) < 2π * 1e-8u"Hz"

    # A transverse field mixes neighbouring F levels within the manifold…
    @test abs(h[stateindex(basis, "S_1/2 F=3", 2), stateindex(basis, "S_1/2 F=4", 3)]) >
          2π * 1u"Hz"
    # …while a longitudinal field is diagonal in m_F.
    hz = zeeman_hamiltonian(ca43, basis, [0.0, 0.0, 0.5]u"mT")
    for i in 1:length(basis), k in 1:length(basis)
        if basis[i].m != basis[k].m
            @test abs(hz[i, k]) < 2π * 1e-9u"Hz"
        end
    end

    # Different fine-structure manifolds do not couple, and scrambled/partial
    # bases pick out the corresponding elements of the full matrix.
    scrambled = StateBasis([
        StateSpec("D_5/2 F=4", 3),
        StateSpec("S_1/2 F=4", 3),
        StateSpec("S_1/2 F=3", 2),
    ])
    hs = zeeman_hamiltonian(ca43, scrambled, b)
    @test iszero(hs[1, 2])
    @test iszero(hs[1, 3])
    @test hs[3, 2] ==
          h[stateindex(basis, "S_1/2 F=3", 2), stateindex(basis, "S_1/2 F=4", 3)]

    # A hyperfine basis must not silently pass through the no-hyperfine code
    # path (which would read m_F as m_J).
    @test_throws MethodError zeeman_hamiltonian(sr88, basis, b)
end

@testitem "Hellmann–Feynman sensitivity" tags=[:unit, :fast] begin
    using FiniteDifferences
    using Unitful

    s = StateSpec("S_1/2 F=4", 4)
    d = StateSpec("D_5/2 F=4", 3)

    # The Hellmann–Feynman expectation value matches numerical differentiation
    # of the exact eigen-energies. (Differentiating the manifold energies
    # rather than the full optical transition_frequency keeps the finite
    # differences clear of the ~4e14 Hz offset, whose float eps of ~0.06 Hz
    # would dominate the derivative otherwise; rtol set by the FD noise floor
    # on the remaining GHz-scale energies.)
    energy_hz(manifold, state, b) = ustrip(
        u"Hz",
        state_energy(hyperfine_manifold(ca43, manifold, b * u"mT"), state) / 2π,
    )
    for b_mt in (0.3, 0.496)
        χ_hf = ustrip(u"Hz/mT", zeeman_sensitivity(ca43, s, d, b_mt * u"mT") / 2π)
        f(b) = energy_hz("D_5/2", d, b) - energy_hz("S_1/2", s, b)
        χ_fd = central_fdm(5, 1)(f, b_mt)
        @test χ_hf ≈ χ_fd rtol = 1e-6
    end

    # Pair form, and a same-manifold (microwave) transition.
    pair = StateSpec("S_1/2 F=4", 0) => StateSpec("S_1/2 F=3", 1)
    @test zeeman_sensitivity(ca43, pair, 0.3u"mT") ==
          zeeman_sensitivity(ca43, pair.first, pair.second, 0.3u"mT")
    f_mw(b) = energy_hz("S_1/2", pair.second, b) - energy_hz("S_1/2", pair.first, b)
    @test ustrip(u"Hz/mT", zeeman_sensitivity(ca43, pair, 14.6u"mT") / 2π) ≈
          central_fdm(5, 1)(f_mw, 14.6) rtol = 1e-6

    # Zero-field limit of the frequency reproduces the ground-state splitting.
    @test Levels.transition_frequency(
        ca43,
        StateSpec("S_1/2 F=4", 0),
        StateSpec("S_1/2 F=3", 0),
        1e-9u"mT",
    ) ≈ 2π * 3225.60828640u"MHz" rtol = 1e-9
end

@testitem "Field-insensitive points" tags=[:unit, :fast] begin
    using Logging
    using Unitful

    s = StateSpec("S_1/2 F=4", 4)
    d = StateSpec("D_5/2 F=4", 3)

    # 729 nm clock transition: both zeros quoted as 3.38 G and 4.96 G with
    # ±16 kHz/G² curvature in [Benhelm2007] (slightly older constants).
    upper = insensitive_field(ca43, s => d, (0.45u"mT", 0.55u"mT"))
    @test upper ≈ 0.496u"mT" rtol = 1e-2
    lower = insensitive_field(ca43, s => d, (0.30u"mT", 0.40u"mT"))
    @test lower ≈ 0.338u"mT" rtol = 1e-2

    # Microwave clock qubits: field-independent points at 146.0942 G and
    # 287.7827 G (atomic_physics test_ca43.py, after T. Harty's DPhil thesis).
    p1 = insensitive_field(
        ca43,
        StateSpec("S_1/2 F=4", 0) => StateSpec("S_1/2 F=3", 1),
        (14.0u"mT", 15.0u"mT"),
    )
    @test abs(p1 - 14.60942u"mT") < 1e-5u"mT"
    p2 = insensitive_field(
        ca43,
        StateSpec("S_1/2 F=4", 1) => StateSpec("S_1/2 F=3", 1),
        (28.0u"mT", 29.0u"mT"),
    )
    @test abs(p2 - 28.77827u"mT") < 1e-5u"mT"

    # A bracket without a sign change is rejected (silencing the solver's own
    # warning about the non-enclosing interval).
    @test_throws ArgumentError with_logger(NullLogger()) do
        insensitive_field(ca43, s => d, (0.6u"mT", 0.7u"mT"))
    end
end
