@testitem "Static polarisabilities vs Jiang (2009)" tags=[:unit, :fast] begin
    using Unitful

    # The stored data is a reduction of Tables 1 and 3 of [Jiang2009] (see
    # species_data.jl), so the static limit must reproduce that paper's totals:
    # α₀(5s) = 91.30, α₀(4d₅/₂) = 62.0 and α₂(4d₅/₂) = −47.7 a.u.
    au = Levels.POLARIZABILITY_AU
    dc = 0.0u"s^-1"
    @test scalar_polarisability(sr88, "S_1/2", dc) / au ≈ 91.30 rtol = 1e-3
    @test scalar_polarisability(sr88, "D_5/2", dc) / au ≈ 62.0 rtol = 1e-3
    @test tensor_polarisability(sr88, "D_5/2", dc) / au ≈ -47.7 rtol = 1e-3

    # S₁/₂ has no oriented sublevels, hence no tensor polarisability, and the
    # vector polarisability vanishes in the static limit for every level.
    @test iszero(tensor_polarisability(sr88, "S_1/2", dc))
    @test iszero(vector_polarisability(sr88, "S_1/2", dc))
    @test iszero(vector_polarisability(sr88, "D_5/2", dc))

    # The differential static scalar polarisability of the clock transition has
    # been measured to −4.8314(20)×10⁻⁴⁰ J m²/V² = −29.303(12) a.u.:
    # T. Lindvall, K. J. Hanhijärvi, T. Fordell, and A. E. Wallin, Phys. Rev.
    # Lett. 135, 043402 (2025), doi:10.1103/52by-28mr.
    Δα =
        scalar_polarisability(sr88, "D_5/2", dc) -
        scalar_polarisability(sr88, "S_1/2", dc)
    @test uconvert(u"C*m^2/V", Δα) ≈ -4.8314e-40u"C*m^2/V" rtol = 3e-3
end

@testitem "674 nm clock light shift vs Lindvall (2025)" tags=[:unit, :fast] begin
    using Unitful

    # T. Lindvall et al., "⁸⁸Sr⁺ optical clock with 7.9 × 10⁻¹⁹ systematic
    # uncertainty…", Phys. Rev. Applied 24, 044082 (2025),
    # doi:10.1103/cztf-bfvp, Sec. III F 1: for linear polarisation at 38° to the
    # magnetic field, the clock laser shifts the (±½ → ±½, ±½ → ±³⁄₂, ±½ → ±⁵⁄₂)
    # Zeeman pairs by (0.77, 0.74, 0.67) mHz/(W/m²).
    λ = 674.025u"nm" # 444 779 044 095 485 Hz
    θ = deg2rad(38.0)
    ε = [complex(sin(θ)), 0.0im, complex(cos(θ))]
    intensity = 1.0u"W/m^2"

    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], λ)
    for (m, reference) in zip((1//2, 3//2, 5//2), (0.77, 0.74, 0.67))
        transition = StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", m)
        shift = light_shift(coefficients, transition, intensity, ε)
        @test shift ≈ 2π * reference * u"mHz" rtol = 2e-2

        # The Zeeman pairs are quoted as pair averages; with a linear
        # polarisation the shift is the same for both signs of m by symmetry.
        mirrored = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -m)
        @test light_shift(coefficients, mirrored, intensity, ε) ≈ shift

        # The one-shot form must agree with the precomputed one.
        @test light_shift(sr88, transition.first, transition.second, λ, intensity, ε) ≈
              shift
    end

    # Appendix B of the same paper gives the differential scalar polarisability
    # at the 1092-nm repumper wavelength as 329(5) a.u.
    λ_repump = uconvert(u"nm", 2π * u"c" / (2π * 274.59u"THz"))
    Δα =
        scalar_polarisability(sr88, "D_5/2", λ_repump) -
        scalar_polarisability(sr88, "S_1/2", λ_repump)
    @test Δα / Levels.POLARIZABILITY_AU ≈ 329 rtol = 2e-2
end

@testitem "Light shift vs scalar/vector/tensor decomposition" tags=[:unit, :fast] begin
    using Unitful

    # The per-Δm-channel sum that light_shift evaluates must be equivalent to the
    # textbook decomposition, for arbitrary elliptical polarisation:
    #   α = α₀ + 𝒜 (m/2j) α₁ + (3|ε₀|² − 1)/2 (3m² − j(j+1))/(j(2j−1)) α₂.
    λ = 674.025u"nm"
    _, ε = beam_vectors(0.4, 0.9, 1.1) # generic elliptical polarisation
    w = Levels.polarisation_weights(ε)
    circular = w[3] - w[1] # 𝒜 = |ε₋₁|² − |ε₊₁|²

    for level in ("S_1/2", "D_5/2")
        j = convert(NoHyperfineNumberSpec, level).j
        α0 = scalar_polarisability(sr88, level, λ)
        α1 = vector_polarisability(sr88, level, λ)
        α2 = tensor_polarisability(sr88, level, λ)
        for m in (-j):j
            explicit = sum(
                w .* Levels.state_polarisabilities(
                    sr88,
                    StateSpec(level, m),
                    Levels.photon_energy(λ),
                ),
            )
            tensor = j > 1//2 ? (3m^2 - j * (j + 1)) / (j * (2j - 1)) : 0//1
            @test explicit ≈
                  α0 + circular * (m / 2j) * α1 + (3w[2] - 1) / 2 * tensor * α2
        end
    end
end

@testitem "Light shift signs and geometry" tags=[:unit, :fast] begin
    using Unitful

    λ = 674.025u"nm"
    intensity = 1e4u"W/m^2"
    basis = StateBasis(["S_1/2", "D_5/2"])
    coefficients = LightShiftCoefficients(sr88, basis, λ)

    # 674 nm is far red-detuned from the 422/408-nm S₁/₂ → P transitions, so the
    # ground state is pulled down; it is blue-detuned from the 1033-nm
    # D₅/₂ → P₃/₂ transition, which dominates D₅/₂ and pushes it up. Both add up
    # to a positive shift of the clock transition.
    _, ε = beam_vectors(π / 2, π / 4)
    s = light_shift(coefficients, StateSpec("S_1/2", 1//2), intensity, ε)
    d = light_shift(coefficients, StateSpec("D_5/2", 5//2), intensity, ε)
    @test s < zero(s)
    @test d > zero(d)
    @test light_shift(
        coefficients,
        StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 5//2),
        intensity,
        ε,
    ) ≈ d - s

    # Only the polarisation enters an E1 light shift, not the beam direction —
    # unlike the E2 coupling, which is why the shift and the Rabi frequency
    # constrain different combinations of the laser parameters. Two beams sharing
    # a polarisation but travelling in different directions are therefore
    # indistinguishable here, while their couplings differ.
    ε_d = [1.0 + 0im, 0.0im, 1.0 + 0im] / sqrt(2)
    n_1, n_2 = [0.0, 1.0, 0.0], [1.0, 0.0, -1.0] / sqrt(2)
    lower, upper = StateSpec("S_1/2", 1//2), StateSpec("D_5/2", 3//2)
    @test rabi_frequency(sr88, lower, upper, intensity, ε_d, n_1) ≉
          rabi_frequency(sr88, lower, upper, intensity, ε_d, n_2)

    # An overall phase of ε is a gauge choice and cannot matter.
    @test light_shift(coefficients, upper, intensity, cis(0.7) * ε_d) ≈
          light_shift(coefficients, upper, intensity, ε_d)

    # The scale of ε is irrelevant, as the intensity is given separately.
    @test light_shift(coefficients, StateSpec("D_5/2", 3//2), intensity, 7ε) ≈
          light_shift(coefficients, StateSpec("D_5/2", 3//2), intensity, ε)

    # Shifts are linear in the intensity.
    @test light_shift(coefficients, StateSpec("D_5/2", 3//2), 3intensity, ε) ≈
          3 * light_shift(coefficients, StateSpec("D_5/2", 3//2), intensity, ε)

    # σ⁺ and σ⁻ light shift ±m oppositely relative to the polarisation-averaged
    # value; that difference is what the vector polarisability describes.
    _, σ_plus = beam_vectors(0.0, π / 4, π / 2)
    _, σ_minus = beam_vectors(0.0, π / 4, -π / 2)
    up = light_shift(coefficients, StateSpec("S_1/2", 1//2), intensity, σ_plus)
    down = light_shift(coefficients, StateSpec("S_1/2", -1//2), intensity, σ_plus)
    @test up ≈ light_shift(coefficients, StateSpec("S_1/2", -1//2), intensity, σ_minus)
    @test down ≈ light_shift(coefficients, StateSpec("S_1/2", 1//2), intensity, σ_minus)
    @test !(up ≈ down)
end

@testitem "Light shift API" tags=[:unit, :fast] begin
    using Unitful

    λ = 674.025u"nm"
    intensity = 1e3u"W/m^2"
    _, ε = beam_vectors(π / 2, 0.3, 0.2)

    # Wavelengths and angular frequencies are interchangeable.
    ω = 2π * u"c" / λ
    @test LightShiftCoefficients(sr88, ["S_1/2"], ω).polarisabilities ≈
          LightShiftCoefficients(sr88, ["S_1/2"], λ).polarisabilities

    # Batch evaluation matches the single-transition form and preserves order.
    transitions = state_pairs("S_1/2", "D_5/2"; Δm=[-2, -1, 1, 2])
    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], λ)
    shifts = light_shift(coefficients, transitions, intensity, ε)
    @test length(shifts) == length(transitions)
    @test all(
        shifts[i] ≈ light_shift(coefficients, transitions[i], intensity, ε) for
        i in eachindex(transitions)
    )

    # Levels without polarisability data are rejected rather than silently
    # treated as unshifted.
    @test isnothing(level_polarisability(sr88, "D_3/2"))
    @test_throws ArgumentError LightShiftCoefficients(sr88, ["D_3/2"], λ)
    @test_throws ArgumentError scalar_polarisability(sr88, "P_1/2", λ)
    @test_throws ArgumentError light_shift(
        sr88,
        StateSpec("D_3/2", 1//2),
        λ,
        intensity,
        ε,
    )

    # A zero polarisation vector carries no information about the geometry.
    @test_throws ArgumentError light_shift(
        coefficients,
        StateSpec("S_1/2", 1//2),
        intensity,
        zeros(ComplexF64, 3),
    )
end

@testitem "Reduced dipole matrix elements vs Einstein A coefficients" tags=[
    :unit,
    :fast,
] begin
    using Unitful

    # The polarisability matrix elements come from a different literature source
    # ([Jiang2009]) than the Einstein A coefficients of the same species, so the
    # two can be cross-checked against each other via
    # A = ω³ |⟨j'‖d‖j⟩|² / (3π ε₀ ħ c³ (2j'+1)).
    #
    # The tolerances are those of the stored A coefficients: the two P → S decays
    # are known to ≈1%, and agree to better than that, whereas D₅/₂ → P₃/₂ is
    # only 8.7(15) µs⁻¹ in [Sansonetti2012] — for that channel the matrix element
    # is by far the better determined of the two.
    for (lower, upper, tolerance) in
        (("S_1/2", "P_1/2", 0.01), ("S_1/2", "P_3/2", 0.015), ("D_5/2", "P_3/2", 0.18))
        d = level_polarisability(sr88, lower).reduced_dipoles[convert(
            NoHyperfineNumberSpec,
            upper,
        )]
        ω = Levels.transition_frequency(sr88, lower, upper)
        j_upper = convert(NoHyperfineNumberSpec, upper).j
        a = ω^3 * d^2 / (3π * u"ε0" * u"ħ" * u"c"^3 * (2 * j_upper + 1))
        @test uconvert(u"µs^-1", a) ≈ einstein_a(sr88, lower, upper) rtol = tolerance
    end
end
