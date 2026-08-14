@testitem "Clebsch–Gordan factors" tags=[:unit, :fast] begin
    # E1 vs E2 from the level parities.
    @test Levels.multipole_rank("S_1/2", "P_3/2") == 1
    @test Levels.multipole_rank("D_3/2", "P_1/2") == 1
    @test Levels.multipole_rank("S_1/2", "D_5/2") == 2

    # Stretch transitions carry unit Clebsch–Gordan factors.
    @test Levels.clebsch_gordan(StateSpec("S_1/2", 1//2), StateSpec("P_3/2", 3//2)) ≈ 1
    @test Levels.clebsch_gordan(StateSpec("S_1/2", 1//2), StateSpec("D_5/2", 5//2)) ≈ 1

    # Non-stretch entries against table values: ⟨½ ½; 1 0 | ³⁄₂ ½⟩ = √(2/3),
    # ⟨½ -½; 2 1 | ⁵⁄₂ ½⟩ = √(2/5).
    @test Levels.clebsch_gordan(StateSpec("S_1/2", 1//2), StateSpec("P_3/2", 1//2)) ≈
          sqrt(2 / 3)
    @test Levels.clebsch_gordan(StateSpec("S_1/2", -1//2), StateSpec("D_5/2", 1//2)) ≈
          sqrt(2 / 5)

    # Δm beyond the multipole rank carries no amplitude.
    @test iszero(
        Levels.clebsch_gordan(StateSpec("S_1/2", -1//2), StateSpec("P_3/2", 3//2)),
    )
end

@testitem "E1 Rabi frequency vs saturation-intensity relation" tags=[:unit, :fast] begin
    using Unitful

    # A stretch transition driven by σ⁺ light is an effective two-level system
    # with the full line strength, so the textbook relation I/I_sat = 2 Ω²/Γ²
    # with I_sat = π h c Γ/(3 λ³) (e.g. Foot, "Atomic Physics" (2005), §7.6.1;
    # Metcalf & van der Straten, "Laser Cooling and Trapping" (1999), §2.4)
    # fixes the absolute scale of Ω.
    a = einstein_a(sr88, "S_1/2", "P_3/2")
    λ = 2π * u"c" / Levels.transition_frequency(sr88, "S_1/2", "P_3/2")
    isat = π * u"h" * u"c" * a / (3 * λ^3)

    intensity = 250.0u"W/m^2"
    n, ε = beam_vectors(0.0, π / 4, π / 2) # σ⁺ along the quantisation axis
    Ω = rabi_frequency(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("P_3/2", 3//2),
        intensity,
        ε,
        n,
    )
    @test Ω ≈ a * sqrt(intensity / (2 * isat)) rtol = 1e-9

    # Unknown transitions are rejected.
    @test_throws ArgumentError rabi_frequency(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("S_1/2", -1//2),
        intensity,
        ε,
        n,
    )
end

@testitem "E2 Rabi frequency vs James (1998)" tags=[:unit, :fast] begin
    using Unitful
    using WignerSymbols: wigner3j

    # Absolute quadrupole Rabi frequency for the ⁸⁸Sr⁺ S₁/₂(m = ½) → D₅/₂(m' =
    # ⁵⁄₂) stretch component at beam angle φ = π/2 and (linear) polarisation
    # angle γ = π/4, against D. F. V. James, Appl. Phys. B 66, 181 (1998), Eqs.
    # (5.11)–(5.13) [arXiv:quant-ph/9702053 numbering]:
    #
    #   Ω = (e |E| / (ħ √(c α))) √(A/k³) σ⁽ᴱ²⁾,
    #   σ⁽ᴱ²⁾ = √(15 (2j'+1)/4) |Σ_q (j 2 j'; -m q m')₃ⱼ c⁽q⁾_il ε_i n_l|,
    #
    # where |E| is the amplitude of one running-wave component of James's
    # standing wave; an ion at the peak field of a single running wave of
    # intensity I sees |E| = √(2 I/(ε₀ c)), for which his 2Ω₀-splitting
    # Hamiltonian gives the conventional sin²(Ω t/2) Rabi flopping.
    lower = StateSpec("S_1/2", 1//2)
    upper = StateSpec("D_5/2", 5//2)
    a = einstein_a(sr88, "S_1/2", "D_5/2")
    ω = Levels.transition_frequency(sr88, "S_1/2", "D_5/2")
    intensity = 1e4u"W/m^2"
    n, ε = beam_vectors(π / 2, π / 4)
    @test all(isreal, ε) # so James's e^{+iωt} field convention needs no conjugation

    # Only q = m - m' = -2 contributes; c⁽⁻²⁾ from James's Eq. (A.12).
    c2 = [1 im 0; im -1 0; 0 0 0] / sqrt(6)
    σ =
        sqrt(15 * (2 * (5 / 2) + 1) / 4) *
        abs(wigner3j(1//2, 2, 5//2, -1//2, -2, 5//2) * (transpose(ε) * c2 * n))
    α = u"q"^2 / (4π * u"ε0" * u"ħ" * u"c")
    E0 = sqrt(2 * intensity / (u"ε0" * u"c"))
    k = ω / u"c"
    Ω_james = u"q" * E0 / (u"ħ" * sqrt(u"c" * α)) * sqrt(a / k^3) * σ
    Ω = rabi_frequency(sr88, lower, upper, intensity, ε, n)
    @test Ω ≈ Ω_james rtol = 1e-9

    # Equivalent compact form Ω = Λ g √(15 A λ³ I/(2π h c)) with the geometry
    # factor g⁽±²⁾ = (1/√6)|cos(γ) sin(2φ)/2 + i sin(γ) sin(φ)| = 1/√12 and
    # Clebsch–Gordan factor Λ = 1 [Roos, PhD thesis, Innsbruck (2000)].
    λ = 2π * u"c" / ω
    @test Ω ≈ sqrt(15 * a * λ^3 * intensity / (2π * u"h" * u"c")) / sqrt(12) rtol = 1e-9

    # Only the directions of ε and n matter; the field amplitude is fixed by
    # the intensity, so any scale (and an overall phase) must drop out.
    @test rabi_frequency(sr88, lower, upper, intensity, cis(0.3) * 5ε, 3n) ≈ Ω
end

@testitem "Rabi frequencies vs coupling-matrix pipeline" tags=[:unit, :fast] begin
    using Unitful

    # Normalising the relative coupling matrix to the computed carrier Rabi
    # frequency of one probed transition must reproduce rabi_frequency for all
    # components, including for elliptical polarisation (Δm = ±q asymmetry).
    basis = StateBasis(["S_1/2", "D_5/2"])
    n, ε = beam_vectors(0.3, 0.7, 0.4)
    C = quadrupole_couplings(basis, "S_1/2", "D_5/2", ε, n)

    intensity = 5e4u"W/m^2"
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    Ω0 = rabi_frequency(sr88, probe.first, probe.second, intensity, ε, n)
    L = rabi_normalised(C, basis, probe, Ω0)
    for (lower, upper) in state_pairs("S_1/2", "D_5/2"; Δm=-2:2)
        @test abs(L[stateindex(basis, upper), stateindex(basis, lower)]) ≈
              rabi_frequency(sr88, lower, upper, intensity, ε, n) atol = 1e-12 * abs(Ω0)
    end
end

@testitem "Rabi normalisation" tags=[:unit, :fast] begin
    using Unitful

    basis = StateBasis(["S_1/2", "D_5/2"])
    # k ⊥ B₀, γ = 45°: Δm = 0 suppressed, Δm = ±1, ±2 present.
    n, ε = beam_vectors(π / 2, π / 4)
    C = quadrupole_couplings(basis, "S_1/2", "D_5/2", ε, n)

    Ω0 = 2π * 100.0u"kHz"
    probe = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2)
    L = rabi_normalised(C, basis, probe, Ω0)
    @test abs(L[stateindex(basis, probe.second), stateindex(basis, probe.first)]) ≈ Ω0

    # Relative amplitudes are preserved.
    i = stateindex(basis, "D_5/2", -1//2)
    k = stateindex(basis, "S_1/2", 1//2)
    i2 = stateindex(basis, probe.second)
    k2 = stateindex(basis, probe.first)
    @test L[i, k] / L[i2, k2] ≈ C[i, k] / C[i2, k2]

    # Normalising to a geometry-suppressed (Δm = 0) transition is refused.
    @test_throws ArgumentError rabi_normalised(
        C,
        basis,
        StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 1//2),
        Ω0,
    )
end

@testitem "Hyperfine reduction factor" tags=[:unit, :fast] begin
    using Unitful
    using WignerSymbols
    using Levels: clebsch_gordan, coupling_transform, hyperfine_reduction, jz_matrix

    # Definitive convention check: conjugating the product-basis electronic
    # Wigner–Eckart amplitudes ⟨J' m_J + q|T^k_q|J m_J⟩ = ⟨J m_J; k q|J' m_J'⟩
    # (identity on the nuclear factor) with the Clebsch–Gordan unitaries must
    # reproduce ⟨F m; k q|F' m'⟩ β^(k)(F → F') for every component — pinning
    # both the 6-j expression and its phase to the coupling_transform basis
    # convention.
    i_nuc = ca43.nuclear_spin
    n_i = Int(2 * i_nuc + 1)
    for (lo_str, hi_str, k) in (
        ("S_1/2", "P_1/2", 1),
        ("S_1/2", "P_3/2", 1),
        ("S_1/2", "D_5/2", 2),
        ("D_3/2", "P_1/2", 1),
    )
        j_lo = NoHyperfineNumberSpec(lo_str).j
        j_hi = NoHyperfineNumberSpec(hi_str).j
        d_lo = Int(2j_lo + 1)
        d_hi = Int(2j_hi + 1)
        u_lo = coupling_transform(ca43, lo_str)
        u_hi = coupling_transform(ca43, hi_str)
        lo_states = collect(StateBasis(ca43, lo_str))
        hi_states = collect(StateBasis(ca43, hi_str))

        maxdev = 0.0
        for q in (-k):k
            m = zeros(n_i * d_hi, n_i * d_lo)
            for i_i in 1:n_i
                for (a, m_j) in enumerate((-j_lo):j_lo),
                    (b, m_jp) in enumerate((-j_hi):j_hi)

                    m_jp == m_j + q || continue
                    m[(b-1)*n_i+i_i, (a-1)*n_i+i_i] =
                        Float64(clebschgordan(j_lo, m_j, k, q, j_hi, m_jp))
                end
            end
            amplitudes = u_hi' * m * u_lo
            for (ci, ls) in enumerate(lo_states), (ck, us) in enumerate(hi_states)
                expected = if us.m == ls.m + q
                    Float64(clebschgordan(ls.level.f, ls.m, k, q, us.level.f, us.m)) * hyperfine_reduction(ca43, ls.level, us.level; rank=k)
                else
                    0.0
                end
                maxdev = max(maxdev, abs(amplitudes[ck, ci] - expected))
            end
        end
        @test maxdev < 1e-12
    end

    # Line-strength sum rule: Σ_F β² = 1 for every upper F' — each hyperfine
    # sublevel decays at the full fine-structure rate, β² being the branching
    # fractions (which is why einstein_as must never hold F-resolved entries).
    for (lo_str, hi_str, k) in (("S_1/2", "P_3/2", 1), ("S_1/2", "D_5/2", 2))
        for hi in hyperfine_levels(ca43, hi_str)
            total = sum(
                hyperfine_reduction(ca43, lo, hi; rank=k)^2 for
                lo in hyperfine_levels(ca43, lo_str)
            )
            @test total ≈ 1.0 rtol = 1e-12
        end
    end

    # Stretched-to-stretched E2 goes through a single path, and the F = F'
    # amplitude has a simple closed form.
    @test hyperfine_reduction(ca43, "S_1/2 F=4", "D_5/2 F=6") ≈ 1.0 rtol = 1e-12
    @test hyperfine_reduction(ca43, "S_1/2 F=4", "D_5/2 F=4") ≈ 3 / (2 * sqrt(5)) rtol =
        1e-12
    # Triangle-forbidden pairs vanish.
    @test iszero(hyperfine_reduction(ca43, "S_1/2 F=4", "D_5/2 F=1"))

    # I → 0 limit: β ≡ +1 exactly, so the hyperfine amplitude degenerates —
    # including its sign — to the plain fine-structure Clebsch–Gordan form.
    toy = HyperfineOneElectronSpecies(;
        mass=1.0u"u",
        nuclear_spin=0//1,
        nuclear_g=0.0,
        energies=Dict(
            convert(NoHyperfineNumberSpec, k) => v for
            (k, v) in ["S_1/2" => 0.0u"J", "D_5/2" => u"h" * 411_000_000.0u"MHz"]
        ),
        hyperfine=Dict(
            convert(NoHyperfineNumberSpec, k) => HyperfineConstants(; a=0.0u"J") for
            k in ["S_1/2", "D_5/2"]
        ),
        einstein_as=Dict(
            convert(Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec}, k) => v for
            (k, v) in [("S_1/2", "D_5/2") => 1.0u"s^-1"]
        ),
    )
    @test hyperfine_reduction(toy, "S_1/2 F=1/2", "D_5/2 F=5/2") ≈ 1.0 atol = 1e-14
    for m_l in (-1//2):(1//2), m_u in (-5//2):(5//2)
        @test transition_amplitude(
            toy,
            StateSpec("S_1/2 F=1/2", m_l),
            StateSpec("D_5/2 F=5/2", m_u),
        ) ≈ clebsch_gordan(StateSpec("S_1/2", m_l), StateSpec("D_5/2", m_u)) atol =
            1e-14
    end

    # Kind guards.
    @test_throws ArgumentError transition_amplitude(
        ca43,
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", 3//2),
    )
    @test_throws ArgumentError hyperfine_reduction(ca43, "S_1/2", "D_5/2")
end

@testitem "Hyperfine Rabi frequency" tags=[:unit, :fast] begin
    using Unitful

    s = StateSpec("S_1/2 F=4", 4)
    d = StateSpec("D_5/2 F=4", 3)
    n, ε = beam_vectors(π / 2, π / 4)
    intensity = 250.0u"W/m^2"

    # Independent reconstruction of the James form with the F-basis angular
    # factor and the hyperfine-resolved transition frequency.
    a = einstein_a(ca43, "S_1/2", "D_5/2")
    ω = Levels.transition_frequency(ca43, s.level, d.level)
    geometry = quadrupole_geometry(ε, n)[Int(d.m-s.m)+3]
    expected = uconvert(
        u"µs^-1",
        sqrt(20π * u"c"^2 * intensity * a / (u"ħ" * ω^3)) *
        abs(transition_amplitude(ca43, s, d) * geometry),
    )
    @test rabi_frequency(ca43, s, d, intensity, ε, n) ≈ expected rtol = 1e-12

    # The at-field form scales the plain one by the exact at-field amplitude
    # ratio, and for a fine-structure species it is the identical (already
    # exact) result.
    B = 0.5u"mT"
    @test rabi_frequency(ca43, s, d, intensity, ε, n, B) ≈
          rabi_frequency(ca43, s, d, intensity, ε, n) *
          abs(transition_amplitude(ca43, s, d, B) / transition_amplitude(ca43, s, d)) rtol =
        1e-12
    @test rabi_frequency(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", 5//2),
        intensity,
        ε,
        n,
        B,
    ) == rabi_frequency(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", 5//2),
        intensity,
        ε,
        n,
    )

    # An undrivable |Δm| > 2 component vanishes identically, and
    # fine-structure states are rejected for a hyperfine species.
    @test iszero(
        rabi_frequency(
            ca43,
            StateSpec("S_1/2 F=4", -4),
            StateSpec("D_5/2 F=4", -1),
            intensity,
            ε,
            n,
        ),
    )
    @test_throws ArgumentError rabi_frequency(
        ca43,
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", 3//2),
        intensity,
        ε,
        n,
    )
end

@testitem "At-field transition amplitudes" tags=[:unit, :fast] begin
    using Unitful
    using WignerSymbols: clebschgordan
    using Levels: coupling_transform

    B = 0.5u"mT"
    m_s = hyperfine_manifold(ca43, "S_1/2", B)
    m_d = hyperfine_manifold(ca43, "D_5/2", B)
    s_states = collect(m_s.basis)
    d_states = collect(m_d.basis)

    # First-principles pin: the electronic E2 tensor components built in the
    # |m_J, m_I⟩ product basis (identity on the nucleus), conjugated into the
    # coupled basis with the Clebsch–Gordan unitaries and then into the field
    # eigenbasis with the manifold eigenvectors, must reproduce every at-field
    # amplitude — independently of the 6-j route the implementation takes —
    # and vanish for Δm ≠ q (m_F stays exact).
    i_nuc = ca43.nuclear_spin
    n_i = Int(2 * i_nuc + 1)
    j_lo = NoHyperfineNumberSpec("S_1/2").j
    j_hi = NoHyperfineNumberSpec("D_5/2").j
    d_lo = Int(2j_lo + 1)
    d_hi = Int(2j_hi + 1)
    u_lo = coupling_transform(ca43, "S_1/2")
    u_hi = coupling_transform(ca43, "D_5/2")
    maxdev = let dev = 0.0
        for q in -2:2
            m = zeros(n_i * d_hi, n_i * d_lo)
            for i_i in 1:n_i
                for (a, m_j) in enumerate((-j_lo):j_lo),
                    (b, m_jp) in enumerate((-j_hi):j_hi)

                    m_jp == m_j + q || continue
                    m[(b-1)*n_i+i_i, (a-1)*n_i+i_i] =
                        Float64(clebschgordan(j_lo, m_j, 2, q, j_hi, m_jp))
                end
            end
            amplitudes = m_d.states' * u_hi' * m * u_lo * m_s.states
            for (ci, ls) in enumerate(s_states), (ck, us) in enumerate(d_states)
                expected = if us.m == ls.m + q
                    transition_amplitude(m_s, m_d, ls, us)
                else
                    0.0
                end
                dev = max(dev, abs(amplitudes[ck, ci] - expected))
            end
        end
        dev
    end
    @test maxdev < 1e-12

    # The one-shot species form matches the manifold-pair form, B = 0 falls
    # back to the zero-field amplitudes, and Δm beyond the rank still carries
    # no amplitude.
    s = StateSpec("S_1/2 F=4", 4)
    d = StateSpec("D_5/2 F=4", 3)
    @test transition_amplitude(ca43, s, d, B) ≈ transition_amplitude(m_s, m_d, s, d) rtol =
        1e-12
    @test transition_amplitude(ca43, s, d, 0.0u"mT") == transition_amplitude(ca43, s, d)
    @test iszero(
        transition_amplitude(
            ca43,
            StateSpec("S_1/2 F=4", -4),
            StateSpec("D_5/2 F=4", -1),
            B,
        ),
    )

    # Low-field continuity (including sign) with the zero-field amplitudes.
    for (lower, upper) in state_pairs("S_1/2 F=4", "D_5/2 F=5"; Δm=-2:2)
        @test transition_amplitude(ca43, lower, upper, 1e-4u"mT") ≈
              transition_amplitude(ca43, lower, upper) atol = 1e-3
    end

    # At 0.5 mT the F mixing moves the ⁴³Ca⁺ clock-component amplitude at the
    # percent level — the rotation is not a silent no-op.
    ratio = transition_amplitude(ca43, s, d, B) / transition_amplitude(ca43, s, d)
    @test 1e-3 < abs(ratio - 1) < 0.2

    # For a fine-structure species a static field along ẑ changes nothing, so
    # the B form is the identical (already exact) amplitude.
    @test transition_amplitude(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", 5//2),
        B,
    ) == transition_amplitude(sr88, StateSpec("S_1/2", 1//2), StateSpec("D_5/2", 5//2))

    # Within one manifold (rank = 1, M1): the same manifold twice; at a small
    # field this recovers the zero-field F-basis amplitude.
    lo = StateSpec("S_1/2 F=4", 0)
    hi = StateSpec("S_1/2 F=3", 1)
    m_s_small = hyperfine_manifold(ca43, "S_1/2", 1e-4u"mT")
    @test transition_amplitude(m_s_small, m_s_small, lo, hi; rank=1) ≈
          transition_amplitude(ca43, lo, hi; rank=1) atol = 1e-4

    # Manifold solutions from different fields are rejected.
    m_d2 = hyperfine_manifold(ca43, "D_5/2", 2 * B)
    @test_throws ArgumentError transition_amplitude(m_s, m_d2, s, d)
end
