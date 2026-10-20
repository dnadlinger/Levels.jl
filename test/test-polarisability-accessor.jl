@testitem "State polarisability accessor" tags=[:unit, :fast] begin
    using Unitful

    _, ε = beam_vectors(0.4, 0.6, 0.8)
    w = Levels.polarisation_weights(ε)
    intensity = 12.0u"W/m^2"
    laser = 422.0u"nm"
    for (species, state) in (
        (sr88, StateSpec("D_5/2", 3//2)),
        (ca43, StateSpec("D_5/2 F=4", 3)),
        (ca43, StateSpec("S_1/2 F=4", 4)),
    )
        level, m = state.level, state.m
        j = level isa HyperfineNumberSpec ? level.f : level.j
        α0 = scalar_polarisability(species, level, laser)
        α1 = vector_polarisability(species, level, laser)
        α2 = tensor_polarisability(species, level, laser)
        tensor = (3m^2 - j * (j + 1)) / (j * (2j - 1))

        # Zero-field decomposition: π is parallel, the σ mean perpendicular,
        # and the σ difference the vector part (j → F for hyperfine levels).
        α = polarisability(species, state, laser)
        @test α[2] ≈ α0 + tensor * α2
        @test (α[1] + α[3]) / 2 ≈ α0 - tensor * α2 / 2
        @test α[3] - α[1] ≈ m / j * α1
        @test polarisability(species, state, laser, [0, 0, 1]) ≈ α[2]
        @test polarisability(species, state, laser, [1, 0, 0]) ≈ (α[1] + α[3]) / 2

        contracted = polarisability(species, state, laser, ε)
        @test contracted ≈ sum(w .* α)
        @test contracted ≈
              α0 + (w[3] - w[1]) * m / (2j) * α1 + (3w[2] - 1) / 2 * tensor * α2
        @test polarisability(species, state, laser, 3ε) ≈ contracted
        @test polarisability(species, state, 2π * u"c" / laser, ε) ≈ contracted
        @test light_shift(species, state, laser, intensity, ε; parts=:background) ≈
              -intensity * contracted / (2 * u"c" * u"ε0" * u"ħ")
    end

    # The 397 nm line is guarded at 0.03 nm (57 GHz), but not 1 nm (2 THz).
    s = StateSpec("S_1/2 F=4", 4)
    λ_397 =
        uconvert(u"nm", 2π * u"c" / Levels.transition_frequency(ca43, "S_1/2", "P_1/2"))
    @test_throws "'S_1/2' → 'P_1/2'" polarisability(ca43, s, λ_397 + 0.03u"nm")
    @test_throws ArgumentError polarisability(ca43, s, λ_397 + 0.03u"nm", ε)
    ω_397 = Levels.transition_frequency(ca43, "S_1/2", "P_1/2")
    @test_throws ArgumentError polarisability(ca43, s, -ω_397, ε)
    ω = 2π * u"c" / laser
    @test polarisability(ca43, s, -ω, ε) ≈ polarisability(ca43, s, ω, conj.(ε))
    @test_throws ArgumentError polarisability(
        ca43,
        s,
        RelativeFrequency("S_1/2 F=4" => "P_1/2 F=4", 2π * 1.0u"GHz"),
        ε,
    )
    @test length(polarisability(ca43, s, λ_397 + 1.0u"nm")) == 3

    state = StateSpec("D_5/2 F=4", 3)
    for B in (0.0u"mT", [0, 0, 1]u"mT", 1u"s", NaN * u"mT", (1 + im)u"mT")
        @test_throws ArgumentError polarisability(ca43, state, 0u"s^-1", ε; B)
    end
    @test_throws ArgumentError polarisability(ca43, state, 0u"s^-1", zeros(3))
    @test_throws ArgumentError polarisability(
        sr88,
        StateSpec("D_3/2", 1//2),
        422u"nm",
        ε,
    )

    # Reject kind mismatches rather than reinterpret m_J as m_F or vice versa.
    @test_throws ArgumentError polarisability(
        ca43,
        StateSpec("D_5/2", 3//2),
        0u"s^-1",
        ε,
    )
    @test_throws ArgumentError polarisability(
        sr88,
        StateSpec("S_1/2 F=2", 2),
        0u"s^-1",
        ε,
    )

    # Fine-structure eigenvectors are unchanged. A large field exposes the
    # otherwise tiny Zeeman correction to the denominators.
    fine = StateSpec("D_5/2", 3//2)
    @test polarisability(sr88, fine, 0u"s^-1", ε; B=0.5u"mT") ≈
          polarisability(sr88, fine, 0u"s^-1", ε) rtol = 1e-6
    B = 10.0u"T"
    w = Levels.polarisation_weights(ε)
    data = level_polarisability(sr88, fine.level)
    j, m = fine.level.j, fine.m
    # In the static limit the co- and counter-rotating terms share the
    # denominator, and together weight the Δm = q component by w_q + w_{-q}.
    channel(upper, d, q) = uconvert(
        u"C*m^2/V",
        d^2 / (2 * upper.j + 1) *
        Levels.dipole_cg(j, m, q, upper.j)^2 *
        (w[q+2] + w[2-q]) / (
            sr88.energies[upper] - sr88.energies[fine.level] +
            u"ħ" * (
                zeeman_shift(sr88, StateSpec(upper, m + q), B) -
                zeeman_shift(sr88, fine, B)
            )
        ),
    )
    expected =
        data.static_scalar +
        data.static_tensor * (3w[2] - 1) / 2 * (3m^2 - j * (j + 1)) / (j * (2j - 1)) +
        sum(
            channel(upper, d, q) for
            (upper, d) in data.reduced_dipoles, q in -1:1 if abs(m + q) <= upper.j
        )
    @test polarisability(sr88, fine, 0u"s^-1", ε; B) ≈ expected rtol = 1e-12
    @test !isapprox(
        polarisability(sr88, fine, 0u"s^-1", ε; B),
        polarisability(sr88, fine, 0u"s^-1", ε);
        rtol=1e-6,
    )
end

@testitem "Relative amplitude matrices" tags=[:unit, :fast] begin
    using Unitful

    # Check E1/E2 amplitudes for both species kinds and fields. The tolerance
    # allows rotation roundoff on entries that vanish by symmetry.
    cache = Dict{NoHyperfineNumberSpec,Levels.ManifoldSolution}()
    for (species, lo, hi, B) in (
        (sr88, "S_1/2", "P_3/2", nothing),
        (sr88, "S_1/2", "D_5/2", 0.5u"mT"),
        (ca43, "D_5/2", "P_3/2", nothing),
        (ca43, "S_1/2", "D_5/2", 0.5u"mT"),
    )
        empty!(cache)
        m_lo = Levels.manifold_solution(species, fine_structure(lo), B, cache)
        m_hi = Levels.manifold_solution(species, fine_structure(hi), B, cache)
        rank = Levels.multipole_rank(lo, hi)
        amplitudes = Levels.relative_amplitudes(species, m_lo, m_hi, rank)
        for (i, s) in enumerate(m_lo.basis), (k, t) in enumerate(m_hi.basis)
            expected =
                isnothing(B) ? transition_amplitude(species, s, t; rank) :
                transition_amplitude(species, s, t, B; rank)
            @test amplitudes[k, i] ≈ expected atol = 1e-12
        end
    end
end

@testitem "Field-resolved polarisability vs product-basis E1 sum" tags=[:unit] begin
    using LinearAlgebra
    using Unitful
    using WignerSymbols

    # Independent construction of dipole operators: electronic CG matrices
    # tensored with the nuclear identity, rather than F-basis CG × β factors.
    function product_sum(state, laser, B, ε)
        fs = fine_structure(state.level)
        source = hyperfine_manifold(ca43, fs, B)
        v = Levels.coupling_transform(ca43, fs) * source.states
        index = stateindex(source.basis, state)
        data = level_polarisability(ca43, fs)
        i_nuc, j = ca43.nuclear_spin, fs.j
        eye_i = Matrix{Float64}(I, Int(2i_nuc + 1), Int(2i_nuc + 1))
        w = Levels.polarisation_weights(ε)
        ħω = Levels.photon_energy(laser)
        α = 0.0u"C*m^2/V"
        for (other, d) in data.reduced_dipoles
            target = hyperfine_manifold(ca43, other, B)
            u = Levels.coupling_transform(ca43, other) * target.states
            for q in -1:1
                electronic = [
                    m_upper == m_lower + q ?
                    Float64(clebschgordan(j, m_lower, 1, q, other.j, m_upper)) :
                    0.0 for m_upper in (-other.j):other.j, m_lower in (-j):j
                ]
                amplitudes = u' * kron(electronic, eye_i) * v[:, index]
                for k in eachindex(amplitudes)
                    Δ =
                        ca43.energies[other] - ca43.energies[fs] +
                        u"ħ" * (target.energies[k] - source.energies[index])
                    α += uconvert(
                        u"C*m^2/V",
                        d^2 / (2 * other.j + 1) *
                        abs2(amplitudes[k]) *
                        (w[q+2] / (Δ - ħω) + w[2-q] / (Δ + ħω)),
                    )
                end
            end
        end
        # Only the unresolved electronic remainder uses the tensor operator.
        tensor =
            j > 1//2 ?
            [
                Float64((3m^2 - j * (j + 1)) / (j * (2j - 1))) for m in (-j):j for
                _ in (-i_nuc):i_nuc
            ] : zeros(size(v, 1))
        α +
        data.static_scalar +
        data.static_tensor * (3w[2] - 1) / 2 * sum(abs2.(v[:, index]) .* tensor)
    end

    _, elliptical = beam_vectors(0.4, 0.6, 0.8)
    for state in (StateSpec("S_1/2 F=4", 4), StateSpec("D_5/2 F=4", 3)),
        laser in (0.0u"s^-1", 422.0u"nm"),
        ε in ([0, 0, 1], [1, im, 0], elliptical)

        B = 0.4958333u"mT"
        α = polarisability(ca43, state, laser, ε; B)
        @test α ≈ product_sum(state, laser, B, ε) rtol = 1e-12
        @test polarisability(ca43, state, laser, ε; B=1e-8u"mT") ≈
              polarisability(ca43, state, laser, ε) rtol = 1e-7
        # Simultaneous field/projection reversal exchanges the circular parts.
        mirror = StateSpec(state.level, -state.m)
        @test polarisability(ca43, mirror, laser, conj.(ε); B=(-B)) ≈ α rtol = 1e-11
        intensity = 10.0u"W/m^2"
        c = LightShiftCoefficients(ca43, [state], laser; B)
        expected = -intensity * α / (2 * u"c" * u"ε0" * u"ħ")
        @test light_shift(c, state, intensity, ε; parts=:background) ≈ expected
        @test light_shift(ca43, state, laser, intensity, ε; B) ≈ expected
    end

    # Partial/scrambled bases address full-manifold eigenstates correctly.
    states = [
        StateSpec("D_5/2 F=4", 3),
        StateSpec("S_1/2 F=4", 4),
        StateSpec("D_5/2 F=3", -2),
    ]
    B, laser = 0.5u"mT", 422.0u"nm"
    c = LightShiftCoefficients(ca43, states, laser; B)
    for state in states
        @test light_shift(c, state, 1u"W/m^2", elliptical; parts=:background) ≈
              -1u"W/m^2" * polarisability(ca43, state, laser, elliptical; B) /
              (2 * u"c" * u"ε0" * u"ħ")
    end
    d = first(states)
    @test !isapprox(
        polarisability(ca43, d, 0u"s^-1", [0, 0, 1]; B),
        polarisability(ca43, d, 0u"s^-1", [0, 0, 1]);
        rtol=1e-3,
    )
end

@testitem "Field polarisability and near-resonant E1 split" tags=[:unit] begin
    using Unitful

    state = StateSpec("S_1/2 F=4", 4)
    B = 0.5u"mT"
    laser = RelativeFrequency("S_1/2 F=4" => "P_1/2 F=4", 2π * 50.0u"GHz")
    c = LightShiftCoefficients(ca43, StateBasis(ca43, "S_1/2", "P_1/2"), laser; B)
    ω = Levels.photon_energy(ca43, laser) / u"ħ"
    _, ε = beam_vectors(0.4, 0.6, 0.8)
    intensity = 10.0u"W/m^2"
    # The split must recover the plain sum on the same eigenstates. Use the
    # internal form because the public accessor guards this 50 GHz detuning.
    w = Levels.polarisation_weights(ε)
    plain = sum(w .* Levels.state_polarisabilities(ca43, state, u"ħ" * ω, B))
    full = -intensity * plain / (2 * u"c" * u"ε0" * u"ħ")
    @test_throws ArgumentError polarisability(ca43, state, ω, ε; B)
    @test light_shift(c, state, intensity, ε) ≈ full rtol = 1e-8
    @test light_shift(c, state, intensity, ε) ≈
          light_shift(c, state, intensity, ε; parts=:background) +
          light_shift(c, state, intensity, ε; parts=:resonant)
end

@testitem "Field polarisability remainder and missing intermediate constants" tags=[
    :unit,
] begin
    using LinearAlgebra
    using Unitful

    # A remainder-only species isolates the off-diagonal-in-F operator.
    fs = fine_structure("D_5/2")
    au = Levels.POLARISABILITY_AU
    remainder = HyperfineOneElectronSpecies(;
        mass=ca43.mass,
        nuclear_spin=ca43.nuclear_spin,
        nuclear_g=ca43.nuclear_g,
        energies=ca43.energies,
        hyperfine=ca43.hyperfine,
        einstein_as=ca43.einstein_as,
        lande_g_overrides=ca43.lande_g_overrides,
        polarisabilities=Dict(
            fs =>
                LevelPolarisability(Dict(); static_scalar=3au, static_tensor=-7au),
        ),
    )
    B, state = 0.5u"mT", StateSpec("D_5/2 F=4", 3)
    manifold = hyperfine_manifold(remainder, fs, B)
    index = stateindex(manifold.basis, state)
    U = Levels.coupling_transform(remainder, fs)
    j, i_nuc = fs.j, remainder.nuclear_spin
    operator = Diagonal([
        3 - 7 * (3m^2 - j * (j + 1)) / (j * (2j - 1)) for m in (-j):j for
        _ in (-i_nuc):i_nuc
    ])
    H0 = ustrip.(u"µs^-1", Levels.manifold_hamiltonian(remainder, fs, B))
    # Symmetric finite perturbation of the full Hamiltonian is independent of
    # taking an expectation in the field eigenvector inside the implementation.
    perturbation = 1e-4
    target_energy = ustrip(u"µs^-1", manifold.energies[index])
    energy_index = argmin(abs.(eigvals(Symmetric(H0)) .- target_energy))
    plus = eigvals(Symmetric(H0 + perturbation * (U' * operator * U)))[energy_index]
    minus = eigvals(Symmetric(H0 - perturbation * (U' * operator * U)))[energy_index]
    @test polarisability(remainder, state, 0u"s^-1", [0, 0, 1]; B) / au ≈
          (plus - minus) / (2perturbation) rtol = 1e-6

    # The established approximation for an intermediate without constants is
    # invariant under changing the basis used for that degenerate manifold.
    missing = HyperfineOneElectronSpecies(;
        mass=ca43.mass,
        nuclear_spin=ca43.nuclear_spin,
        nuclear_g=ca43.nuclear_g,
        energies=ca43.energies,
        einstein_as=ca43.einstein_as,
        hyperfine=Dict(fs => ca43.hyperfine[fs]),
        lande_g_overrides=ca43.lande_g_overrides,
        polarisabilities=ca43.polarisabilities,
    )
    v = Levels.coupling_transform(missing, fs) * manifold.states[:, index]
    source_shift = u"ħ" * manifold.energies[index]
    data = level_polarisability(missing, fs)
    other, d = only(data.reduced_dipoles)
    Δ = missing.energies[other] - missing.energies[fs] - source_shift
    # π electric field: electronic CG strengths, summed over the nuclear index.
    expected =
        sum(
            abs2(v[k]) * Levels.dipole_cg(j, m, 0, other.j)^2 for
            (k, m) in enumerate([m for m in (-j):j for _ in (-i_nuc):i_nuc])
        ) * d^2 / (2 * other.j + 1) * 2 / Δ
    tensor = sum(
        abs2(v[k]) * (3m^2 - j * (j + 1)) / (j * (2j - 1)) for
        (k, m) in enumerate([m for m in (-j):j for _ in (-i_nuc):i_nuc])
    )
    expected += data.static_scalar + data.static_tensor * tensor
    @test polarisability(missing, state, 0u"s^-1", [0, 0, 1]; B) ≈ expected rtol = 1e-12
end
