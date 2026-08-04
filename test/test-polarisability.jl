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
    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], λ)

    # Wavelengths and angular frequencies are interchangeable.
    ω = 2π * u"c" / λ
    @test LightShiftCoefficients(sr88, ["S_1/2"], ω).polarisabilities ≈
          LightShiftCoefficients(sr88, ["S_1/2"], λ).polarisabilities

    # Levels without polarisability data are admitted — the near-resonant
    # quadrupole part does not need it — but their electric-dipole shift is
    # rejected at evaluation rather than silently treated as unshifted.
    @test isnothing(level_polarisability(sr88, "D_3/2"))
    mixed = LightShiftCoefficients(sr88, ["S_1/2", "D_3/2"], λ)
    @test light_shift(mixed, StateSpec("S_1/2", 1//2), intensity, ε) ≈
          light_shift(coefficients, StateSpec("S_1/2", 1//2), intensity, ε)
    @test_throws ArgumentError light_shift(
        mixed,
        StateSpec("D_3/2", 3//2),
        intensity,
        ε,
    )
    @test_throws ArgumentError light_shift(
        mixed,
        StateSpec("S_1/2", 1//2) => StateSpec("D_3/2", 3//2),
        intensity,
        ε,
    )
    @test_throws ArgumentError scalar_polarisability(sr88, "P_1/2", λ)
    @test_throws ArgumentError light_shift(
        sr88,
        StateSpec("D_3/2", 1//2),
        λ,
        intensity,
        ε,
    )

    # The beam direction and the field only make sense together.
    @test_throws ArgumentError light_shift(
        coefficients,
        StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", 3//2),
        intensity,
        ε;
        n=[0.0, 1.0, 0.0],
    )
    @test_throws ArgumentError light_shift(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("D_5/2", 3//2),
        λ,
        intensity,
        ε;
        B=0.5u"mT",
    )

    # A zero polarisation vector carries no information about the geometry.
    @test_throws ArgumentError light_shift(
        coefficients,
        StateSpec("S_1/2", 1//2),
        intensity,
        zeros(ComplexF64, 3),
    )
end

@testitem "Near-resonant quadrupole shift vs exact diagonalisation" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful

    λ = 674.025u"nm"
    B = 0.5u"mT"
    intensity = 500.0u"W/m^2"
    n, ε = beam_vectors(1.1, 0.7, 0.4) # generic elliptical polarisation
    basis = StateBasis(["S_1/2", "D_5/2"])
    coefficients = LightShiftCoefficients(sr88, basis, λ)
    couplings = quadrupole_couplings(basis, "S_1/2", "D_5/2", ε, n)
    zeeman = [zeeman_shift(sr88, state, B) for state in basis]

    # Rotating-frame Hamiltonian of the whole S₁/₂–D₅/₂ manifold with the laser
    # `detuning` away from the probed transition, with the line centre
    # subtracted so that the two probed states sit at zero for zero detuning.
    # Dropping the probed coupling itself — which only drives the Rabi flopping
    # — the observed resonance is where the two dressed states adiabatically
    # connected to the probed pair cross, so this returns their splitting.
    function splitting(probed, detuning)
        lower = stateindex(basis, probed.first)
        upper = stateindex(basis, probed.second)
        Ω = rabi_frequency(sr88, probed.first, probed.second, intensity, ε, n)
        scaled = rabi_normalised(couplings, basis, probed, Ω)
        scaled[upper, lower] = zero(eltype(scaled))
        h = ustrip.(u"µs^-1", (scaled .+ scaled') ./ 2)
        for (i, state) in enumerate(basis)
            reference = state.level == probed.first.level ? lower : upper
            h[i, i] = ustrip(u"µs^-1", zeeman[i] - zeeman[reference])
            state.level == probed.second.level && (h[i, i] -= detuning)
        end
        values, vectors = eigen(Hermitian(h))
        dressed(k) = values[argmax(abs2.(vectors[k, :]))]
        dressed(upper) - dressed(lower)
    end

    for probed in state_pairs("S_1/2", "D_5/2"; Δm=[-2, -1, 0, 1, 2])
        # The splitting falls by exactly one per unit of detuning to the order
        # that matters, so a secant through two points brackets the crossing in
        # one step. Sampling either side of it rather than iterating onto it
        # matters: removing the probed coupling leaves l and u weakly coupled
        # through the other components, so the crossing stays narrowly avoided,
        # and right at its centre the two states are equal mixtures and the
        # identification by largest overlap above stops being meaningful.
        offset = 1e-2
        below, above = splitting(probed, -offset), splitting(probed, offset)
        exact = (-offset + 2offset * below / (below - above)) * u"µs^-1"

        # What second-order perturbation theory drops is the next order in the
        # coupling: fourth-order energies, the two-photon couplings between the
        # probed pair and within each manifold, all of them O(Ω⁴/Δ³). Sweeping
        # the intensity confirms the residual is that and nothing else — it
        # tracks (Ω/Δ)² with a coefficient of order one, here (4 × 10⁻⁴)² ≈
        # 1.6 × 10⁻⁷. Note the reference is what limits how tightly this can be
        # asserted: a shift of ~10⁻⁶ µs⁻¹ is being read off eigenvalues of a
        # matrix with ~10² µs⁻¹ entries, so at intensities much below this one
        # the comparison hits the double-precision floor long before the
        # perturbation theory does.
        @test quadrupole_light_shift(coefficients, probed, intensity, ε; n, B) ≈ exact rtol =
            1e-5
    end
end

@testitem "Near-resonant quadrupole shift vs explicit level sum" tags=[:unit, :fast] begin
    using Unitful

    λ = 674.025u"nm"
    B = 0.32u"mT"
    intensity = 12.0u"W/m^2"
    n, ε = beam_vectors(0.9, 0.6, 0.35)
    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], λ)

    # Second-order perturbation theory written out over the Zeeman components:
    # every other component sharing the probed upper state pushes it by
    # −Ω²/(4Δ), and every other one sharing the probed lower state pushes that
    # by +Ω²/(4Δ), so both enter the transition frequency with the same sign.
    # With the laser on resonance the detunings are Zeeman splittings alone.
    function reference(lower, upper)
        shift = 0.0u"µs^-1"
        for m in (-1//2):(1//2)
            m == lower.m && continue
            other = StateSpec("S_1/2", m)
            Ω = rabi_frequency(sr88, other, upper, intensity, ε, n)
            shift -=
                Ω^2 /
                (4 * (zeeman_shift(sr88, other, B) - zeeman_shift(sr88, lower, B)))
        end
        for m in (-5//2):(5//2)
            m == upper.m && continue
            other = StateSpec("D_5/2", m)
            Ω = rabi_frequency(sr88, lower, other, intensity, ε, n)
            shift -=
                Ω^2 /
                (4 * (zeeman_shift(sr88, upper, B) - zeeman_shift(sr88, other, B)))
        end
        shift
    end

    for (lower, upper) in state_pairs("S_1/2", "D_5/2"; Δm=[-2, -1, 0, 1, 2])
        shift = quadrupole_light_shift(coefficients, lower => upper, intensity, ε; n, B)
        @test shift ≈ reference(lower, upper)

        # The precomputed and one-shot forms must agree, …
        @test quadrupole_light_shift(sr88, lower, upper, intensity, ε; n, B) ≈ shift

        # … and adding the near-resonant part to light_shift must be exactly
        # that: the dipole shift plus the quadrupole one.
        @test light_shift(coefficients, lower => upper, intensity, ε; n, B) ≈
              light_shift(coefficients, lower => upper, intensity, ε) + shift
        @test light_shift(sr88, lower, upper, λ, intensity, ε; n, B) ≈
              light_shift(coefficients, lower => upper, intensity, ε; n, B)
    end
end

@testitem "Near-resonant quadrupole shift symmetries" tags=[:unit, :fast] begin
    using Unitful

    λ = 674.025u"nm"
    B = 0.24u"mT"
    intensity = 8.0u"W/m^2"
    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], λ)
    pairs =
        [StateSpec("S_1/2", 1//2) => StateSpec("D_5/2", m) for m in (1//2, 3//2, 5//2)]
    mirrored(t) = StateSpec("S_1/2", -t.first.m) => StateSpec("D_5/2", -t.second.m)

    # [Lindvall2025], Sec. III F 2: for a perfectly linear polarisation the
    # shifts of the two components of a Zeeman pair are equal and opposite, so
    # they cancel in the pair average that the clock is steered to.
    n_lin, ε_lin = beam_vectors(deg2rad(71.0), deg2rad(33.0))
    for t in pairs
        shift = quadrupole_light_shift(coefficients, t, intensity, ε_lin; n=n_lin, B)
        @test !isapprox(shift, zero(shift), atol=1e-12u"µs^-1")
        @test quadrupole_light_shift(
            coefficients,
            mirrored(t),
            intensity,
            ε_lin;
            n=n_lin,
            B,
        ) ≈ -shift
    end

    # An elliptical polarisation breaks that symmetry and leaves a net shift.
    n_ell, ε_ell = beam_vectors(deg2rad(71.0), deg2rad(33.0), 0.3)
    for t in pairs
        average =
            quadrupole_light_shift(coefficients, t, intensity, ε_ell; n=n_ell, B) +
            quadrupole_light_shift(
                coefficients,
                mirrored(t),
                intensity,
                ε_ell;
                n=n_ell,
                B,
            )
        @test !isapprox(average, zero(average), atol=1e-12u"µs^-1")
    end

    # A σ⁺ beam along the quantisation axis drives the Δm = +1 channel alone,
    # which splits the manifold into independent two-level systems: for the
    # transitions it can drive there is no other component left to couple to.
    n_σ, ε_σ = beam_vectors(0.0, π / 4, π / 2)
    for m in (-1//2, 1//2)
        driven = StateSpec("S_1/2", m) => StateSpec("D_5/2", m + 1)
        @test !iszero(rabi_frequency(sr88, driven..., intensity, ε_σ, n_σ))
        shift = quadrupole_light_shift(coefficients, driven, intensity, ε_σ; n=n_σ, B)
        @test isapprox(shift, zero(shift), atol=1e-12u"µs^-1")
    end

    # The shift is linear in the intensity and inverse in the field, which is
    # signed: reversing the field mirrors the Zeeman structure.
    t = first(pairs)
    shift = quadrupole_light_shift(coefficients, t, intensity, ε_ell; n=n_ell, B)
    @test quadrupole_light_shift(coefficients, t, 3intensity, ε_ell; n=n_ell, B) ≈
          3shift
    @test quadrupole_light_shift(coefficients, t, intensity, ε_ell; n=n_ell, B=B / 4) ≈
          4shift
    @test quadrupole_light_shift(coefficients, t, intensity, ε_ell; n=n_ell, B=(-B)) ≈
          -shift

    # Only the directions of ε and n matter, and an overall phase of ε is a
    # gauge choice — but unlike the dipole shift, the beam direction does
    # enter: the same polarisation sent along two different beam directions
    # gives different shifts.
    @test quadrupole_light_shift(
        coefficients,
        t,
        intensity,
        cis(0.7) * 5ε_ell;
        n=3n_ell,
        B,
    ) ≈ shift
    ε_y = [0.0im, 1.0 + 0.0im, 0.0im]
    along_z = quadrupole_light_shift(coefficients, t, intensity, ε_y; n=[0, 0, 1.0], B)
    along_x = quadrupole_light_shift(coefficients, t, intensity, ε_y; n=[1.0, 0, 0], B)
    @test !isapprox(along_z, along_x)
end

@testitem "674 nm quadrupole light shift vs Lindvall (2025)" tags=[:unit, :fast] begin
    using Unitful

    # [Lindvall2025] never publishes the clock-beam direction: Sec. III F 1 only
    # states the polarisation to be linear at 38° to the magnetic field, and
    # then refers to "our angles" through the intensities that reproduce a Rabi
    # frequency of 2π × 5 Hz on the (±½ → ±½, ±½ → ±³⁄₂, ±½ → ±⁵⁄₂) Zeeman
    # pairs, (2.6, 1.7, 1.9) mW/m². Those four numbers between them leave only
    # φ_k ≈ 71°, γ_pol ≈ 33° in the `beam_vectors` convention — plus the mirror
    # image φ_k ≈ 109°, which none of them can distinguish, and which changes
    # the shift below by under 10%. Note that γ_pol is a rotation out of the
    # (k, ẑ) plane rather than an angle to the field; it is the resulting ε
    # that sits at the quoted 38° to the quantisation axis.
    n, ε = beam_vectors(deg2rad(71.0), deg2rad(33.0))
    @test acosd(abs(ε[3])) ≈ 38 rtol = 2e-2

    lower = StateSpec("S_1/2", 1//2)
    reference = (2.6u"mW/m^2", 1.7u"mW/m^2", 1.9u"mW/m^2")
    intensities = map((1//2, 3//2, 5//2)) do m
        Ω = rabi_frequency(sr88, lower, StateSpec("D_5/2", m), 1.0u"W/m^2", ε, n)
        1.0u"W/m^2" * (2π * 5.0u"Hz" / Ω)^2
    end
    for (intensity, quoted) in zip(intensities, reference)
        @test intensity ≈ quoted rtol = 2e-2
    end

    # At their 4.8-µT field and a pulse area of 1.1π over a 100-ms probe, the
    # individual Zeeman components then shift by up to ≈170 µHz. Sec. III F 2
    # quotes ≈500 µHz, obtained after accounting for resolved sidebands and
    # thermal dephasing, which raise the intensity needed for a given pulse area
    # (they measure a factor of two directly); this only checks the scale.
    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], 674.025u"nm")
    shifts = map(zip((1//2, 3//2, 5//2), intensities)) do (m, intensity)
        upper = StateSpec("D_5/2", m)
        area = 1.1π / rabi_frequency(sr88, lower, upper, intensity, ε, n)
        scaled = intensity * (area / 100u"ms")^2
        quadrupole_light_shift(coefficients, lower => upper, scaled, ε; n, B=4.8u"µT")
    end
    @test 50u"µHz" < maximum(abs, shifts) / 2π < 1u"mHz"
end

@testitem "Near-resonant quadrupole shift API" tags=[:unit, :fast] begin
    using Unitful

    λ = 674.025u"nm"
    B = 0.24u"mT"
    intensity = 8.0u"W/m^2"
    n, ε = beam_vectors(0.9, 0.6, 0.35)
    lower, upper = StateSpec("S_1/2", 1//2), StateSpec("D_5/2", 3//2)
    coefficients = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], λ)

    # The Zeeman components are degenerate without a field, where second-order
    # perturbation theory has nothing to expand in.
    @test_throws ArgumentError quadrupole_light_shift(
        coefficients,
        lower => upper,
        intensity,
        ε;
        n,
        B=0.0u"mT",
    )

    # B is the signed scalar component along the quantisation axis, not the
    # Cartesian field vector zeeman_hamiltonian() accepts.
    @test_throws ArgumentError quadrupole_light_shift(
        coefficients,
        lower => upper,
        intensity,
        ε;
        n,
        B=[0.0, 0.0, 0.5] .* u"mT",
    )

    # Transitions that are not electric-quadrupole ones with known data are
    # rejected rather than silently treated as unshifted — including a pair
    # given the wrong way round, …
    @test_throws ArgumentError quadrupole_light_shift(
        coefficients,
        upper => lower,
        intensity,
        ε;
        n,
        B,
    )
    @test_throws ArgumentError quadrupole_light_shift(
        sr88,
        upper,
        lower,
        intensity,
        ε;
        n,
        B,
    )
    @test_throws ArgumentError quadrupole_light_shift(
        sr88,
        StateSpec("S_1/2", 1//2),
        StateSpec("P_3/2", 3//2),
        intensity,
        ε;
        n,
        B,
    )

    # … and so are Δm = ±3 components, which no beam geometry can drive, so
    # there is no resonance whose shift could be observed.
    undrivable = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", 5//2)
    @test iszero(Levels.clebsch_gordan(undrivable.first, undrivable.second))
    @test_throws ArgumentError quadrupole_light_shift(
        coefficients,
        undrivable,
        intensity,
        ε;
        n,
        B,
    )
    @test_throws ArgumentError quadrupole_light_shift(
        sr88,
        undrivable.first,
        undrivable.second,
        intensity,
        ε;
        n,
        B,
    )
    @test !haskey(
        coefficients.quadrupole_shifts,
        (
            stateindex(coefficients.basis, undrivable.first),
            stateindex(coefficients.basis, undrivable.second),
        ),
    )

    # A vanishing beam direction carries no information about the geometry, and
    # a polarisation that is not transverse to the beam direction does not
    # describe a physical beam.
    @test_throws ArgumentError quadrupole_light_shift(
        coefficients,
        lower => upper,
        intensity,
        ε;
        n=zeros(3),
        B,
    )
    @test_throws ArgumentError quadrupole_light_shift(
        coefficients,
        lower => upper,
        intensity,
        ε;
        n=[0.0, 1.0, 0.0],
        B,
    )

    # The near-resonant coefficients do not depend on the laser frequency, which
    # the model instead fixes to resonance with the probed transition — but for
    # that very reason, evaluating them for a laser inconsistent with that
    # premise is refused.
    off_resonant = LightShiftCoefficients(sr88, ["S_1/2", "D_5/2"], 1092.0u"nm")
    @test off_resonant.quadrupole_shifts == coefficients.quadrupole_shifts
    @test !isempty(coefficients.quadrupole_shifts)
    @test_throws ArgumentError quadrupole_light_shift(
        off_resonant,
        lower => upper,
        intensity,
        ε;
        n,
        B,
    )
    @test_throws ArgumentError light_shift(
        off_resonant,
        lower => upper,
        intensity,
        ε;
        n,
        B,
    )

    # Bases without an E2 pair simply carry no near-resonant data.
    @test isempty(LightShiftCoefficients(sr88, ["S_1/2"], λ).quadrupole_shifts)

    # The quadrupole part needs no E1 polarisability data: the 687-nm
    # S₁/₂ → D₃/₂ line works precomputed just as well as one-shot.
    ω_687 = Levels.transition_frequency(sr88, "S_1/2", "D_3/2")
    d32 = LightShiftCoefficients(sr88, ["S_1/2", "D_3/2"], ω_687)
    t = StateSpec("S_1/2", 1//2) => StateSpec("D_3/2", 3//2)
    @test quadrupole_light_shift(d32, t, intensity, ε; n, B) ≈
          quadrupole_light_shift(sr88, t.first, t.second, intensity, ε; n, B)
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

@testitem "F-basis light shifts" tags=[:unit] begin
    using Unitful
    using Levels: photon_energy, state_polarisabilities

    # Toy hyperfine twins of ⁸⁸Sr⁺ — identical fine-structure data, fictitious
    # nuclear spin — so the F-basis machinery can be validated against the
    # exactly-known J-basis results. Only the levels whose states are queried
    # need hyperfine constants; intermediate levels without them keep their F
    # levels degenerate at the centroid.
    toy(i_nuc, consts) = HyperfineOneElectronSpecies(;
        mass=sr88.mass,
        nuclear_spin=i_nuc,
        nuclear_g=2e-4,
        energies=sr88.energies,
        hyperfine=Dict(convert(NoHyperfineNumberSpec, k) => v for (k, v) in consts),
        einstein_as=sr88.einstein_as,
        polarisabilities=sr88.polarisabilities,
    )
    zero_consts = [
        "S_1/2" => HyperfineConstants(; a=0.0u"J"),
        "D_5/2" => HyperfineConstants(; a=0.0u"J"),
    ]
    ħω = photon_energy(674.0u"nm")

    # I → 0: the single-F states reproduce the fine-structure channel
    # polarisabilities exactly.
    toy0 = toy(0//1, zero_consts)
    for (level, hf_level, j) in
        (("S_1/2", "S_1/2 F=1/2", 1//2), ("D_5/2", "D_5/2 F=5/2", 5//2))
        for m in (-j):j
            @test state_polarisabilities(toy0, StateSpec(hf_level, m), ħω) ≈
                  state_polarisabilities(sr88, StateSpec(level, m), ħω) rtol = 1e-12
        end
    end

    # With degenerate intermediate levels (zeroed constants), the F-basis
    # result is exactly the 6-j re-projection of the J-basis polarisability:
    # the (2F+1)-weighted, channel-averaged manifold mean recovers α₀, and
    # J = 1/2 manifolds have exactly zero tensor component.
    toy32 = toy(3//2, zero_consts)
    for level in ("S_1/2", "D_5/2")
        states = collect(StateBasis(toy32, level))
        mean =
            sum(sum(state_polarisabilities(toy32, s, ħω)) / 3 for s in states) /
            length(states)
        @test mean ≈ scalar_polarisability(sr88, level, 674.0u"nm") rtol = 1e-12
    end
    for f in 1:2
        @test abs(tensor_polarisability(toy32, "S_1/2 F=$f", 674.0u"nm")) <
              1e-15 * abs(scalar_polarisability(sr88, "S_1/2", 674.0u"nm"))
    end

    # The α₀/α₁/α₂ decomposition (with j → F) against the per-channel sum, for
    # an arbitrary elliptical polarisation.
    _, ε = beam_vectors(0.4, 0.6, 0.8)
    w = Levels.polarisation_weights(ε)
    laser = 674.0u"nm"
    for (level, m) in (("D_5/2 F=4", 3//1), ("D_5/2 F=2", -1//1), ("S_1/2 F=2", 2//1))
        spec = HyperfineNumberSpec(level)
        f = spec.f
        α = state_polarisabilities(toy32, StateSpec(spec, m), ħω)
        α0 = scalar_polarisability(toy32, level, laser)
        α1 = vector_polarisability(toy32, level, laser)
        α2 = tensor_polarisability(toy32, level, laser)
        circ = w[3] - w[1]
        tensor_angle = (3 * w[2] - 1) / 2 * (3m^2 - f * (f + 1)) / (f * (2f - 1))
        @test sum(w .* α) ≈ α0 + circ * (m / (2f)) * α1 + tensor_angle * α2 rtol = 1e-12
    end

    # Hyperfine-resolved intermediate detunings give an S_1/2 F level a small
    # but genuinely non-zero tensor polarisability (∝ A_P/Δ) that pure J-basis
    # re-projection cannot produce.
    split_consts = [
        "S_1/2" => HyperfineConstants(; a=0.0u"J"),
        "D_5/2" => HyperfineConstants(; a=0.0u"J"),
        "P_1/2" => HyperfineConstants(; a=u"h" * -100.0u"MHz"),
        "P_3/2" => HyperfineConstants(; a=u"h" * -20.0u"MHz", b=u"h" * -5.0u"MHz"),
    ]
    toy_split = toy(3//2, split_consts)
    α2_split = tensor_polarisability(toy_split, "S_1/2 F=2", laser)
    α0_ref = scalar_polarisability(sr88, "S_1/2", laser)
    @test abs(α2_split) > 1e-9 * abs(α0_ref)
    @test abs(α2_split) < 1e-4 * abs(α0_ref)

    # Hyperfine levels on a species without hyperfine structure fail loudly.
    @test_throws ArgumentError scalar_polarisability(sr88, "S_1/2 F=2", laser)
end

@testitem "Hyperfine near-resonant quadrupole shift API" tags=[:unit] begin
    using Unitful
    using Levels: quadrupole_shift_coefficients

    s = StateSpec("S_1/2 F=4", 4)
    d = StateSpec("D_5/2 F=4", 3)
    laser = uconvert(
        u"nm",
        2π * u"c" * u"ħ" / (u"ħ" * Levels.transition_frequency(ca43, s.level, d.level)),
    )
    B = 0.3u"mT"
    intensity = 10.0u"W/m^2"
    n, ε = beam_vectors(π / 2, π / 4, 0.3)

    # One-shot form against the precomputed coefficients (E1 data is absent for
    # ca43, so the E1 rows are NaN, but the E2 machinery must work regardless).
    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    c = LightShiftCoefficients(ca43, basis, laser; B)
    @test quadrupole_light_shift(c, s => d, intensity, ε; n, B) ≈
          quadrupole_light_shift(ca43, s, d, intensity, ε; n, B) rtol = 1e-12

    # Unlike the fine-structure case ([Lindvall2025] Sec. III F 2), the ±m
    # Zeeman-pair average does *not* cancel for linear polarisation: the
    # spectators in other F levels sit at hyperfine-interval detunings that
    # are even under m → −m (only the intra-F Zeeman detunings are odd), so
    # their contribution survives — and here dominates — the pair average.
    # The rigorous mirror symmetry is instead E(F, m, −B) = E(F, −m, B): the
    # −m pair at +B equals the +m pair at −B.
    s_m = StateSpec("S_1/2 F=4", -4)
    d_m = StateSpec("D_5/2 F=4", -3)
    n_lin, ε_lin = beam_vectors(π / 2, π / 4)
    plus = quadrupole_light_shift(ca43, s, d, intensity, ε_lin; n=n_lin, B)
    minus = quadrupole_light_shift(ca43, s_m, d_m, intensity, ε_lin; n=n_lin, B)
    @test !isapprox(plus + minus, zero(plus); atol=abs(plus))
    @test minus ≈ quadrupole_light_shift(ca43, s, d, intensity, ε_lin; n=n_lin, B=(-B)) rtol =
        1e-9

    # Field-consistency and error paths.
    @test_throws ArgumentError quadrupole_light_shift(
        c,
        s => d,
        intensity,
        ε;
        n,
        B=0.4u"mT",
    )
    @test_throws ArgumentError quadrupole_shift_coefficients(ca43, s, d)
    @test_throws ArgumentError quadrupole_shift_coefficients(ca43, s, d, 0.0u"mT")
    @test_throws ArgumentError quadrupole_shift_coefficients(
        ca43,
        s,
        d,
        [0.0, 0.0, 0.3]u"mT",
    )
    @test_throws ArgumentError quadrupole_light_shift(
        ca43,
        StateSpec("S_1/2 F=4", -4),
        StateSpec("D_5/2 F=6", -1),
        intensity,
        ε;
        n,
        B,
    )
    @test_throws ArgumentError quadrupole_light_shift(ca43, d, s, intensity, ε; n, B)

    # Without B, the constructor precomputes the E1 part only.
    @test isempty(LightShiftCoefficients(ca43, basis, laser).quadrupole_shifts)
end

@testitem "Hyperfine near-resonant quadrupole shift vs exact model" tags=[
    :integration,
    :slow,
] setup=[] begin
    using Unitful
    using Levels.PeriodicDynamics

    # The perturbative κ contraction against the exact resonance position of
    # the laser-probed two-manifold model (monodromy propagation, no ac
    # drives): the shift of the carrier resonance from the off-resonant
    # couplings to the spectator Zeeman components is exactly what the
    # near-resonant E2 model describes. The laser-coupling scale is chosen
    # small enough (Ω₀ ≪ Zeeman splittings) that higher orders are negligible.
    s = StateSpec("S_1/2 F=4", 4)
    d = StateSpec("D_5/2 F=4", 3)
    B = 0.3u"mT"
    n, ε = beam_vectors(π / 2, π / 4, 0.3)
    # Large enough that the ~Ω₀·1e-7-scale resonance-finder tolerance of
    # exact_sideband stays well below the ∝ Ω₀² shift, small enough that
    # Ω₀ ≪ the ~MHz Zeeman splittings keeps higher orders negligible.
    Ω0 = 2π * 10.0u"kHz"

    # Intensity that gives the probed component the carrier Rabi frequency Ω0.
    ref_intensity = 1.0u"W/m^2"
    intensity = ref_intensity * (Ω0 / rabi_frequency(ca43, s, d, ref_intensity, ε, n))^2

    expected = quadrupole_light_shift(ca43, s, d, intensity, ε; n, B)

    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    coupling = rabi_normalised(
        quadrupole_couplings(ca43, basis, "S_1/2", "D_5/2", ε, n),
        basis,
        s => d,
        Ω0,
    )
    dt = DrivenTransition(
        ca43,
        basis,
        s => d;
        drive_frequency=2π * 30.0u"MHz",
        static_field=B,
        coupling,
    )
    exact = exact_sideband(dt; sideband=0, ngrid=16)
    @test exact.δ_res ≈ expected rtol = 2e-2
end
