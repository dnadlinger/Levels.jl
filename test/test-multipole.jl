@testitem "E2 geometry factors" tags=[:unit, :fast] begin
    # Moduli from Roos's thesis / James (1998), up to a common normalisation
    # constant: |Γ_0|, |Γ_±1|, |Γ_±2| for linear polarisation angle γ and beam
    # angle φ to the quantisation axis.
    roos(φ, γ) = (
        abs(cos(γ) * sin(2φ)) / 2,
        abs(cos(γ) * cos(2φ) + im * sin(γ) * cos(φ)) / sqrt(6),
        abs(cos(γ) * sin(2φ) / 2 + im * sin(γ) * sin(φ)) / sqrt(6),
    )
    ratios = Float64[]
    for φ in (0.3, 1.0, π / 2, 2.2), γ in (0.0, 0.7, π / 4, π / 2)
        n, ε = beam_vectors(φ, γ)
        Γ = quadrupole_geometry(ε, n)
        r = roos(φ, γ)
        for (q, rq) in ((0, r[1]), (1, r[2]), (-1, r[2]), (2, r[3]), (-2, r[3]))
            rq > 1e-12 && push!(ratios, abs(Γ[q+3]) / rq)
        end
    end
    @test all(x -> isapprox(x, ratios[1]; rtol=1e-12), ratios)

    # Δm = 0 vanishes for k ⊥ B₀; Δm = ±1, ±2 all present at γ = 45°.
    n, ε = beam_vectors(π / 2, π / 4)
    Γ = quadrupole_geometry(ε, n)
    @test abs(Γ[0+3]) < 1e-14
    @test all(abs(Γ[q+3]) > 0.1 for q in (-2, -1, 1, 2))
end

@testitem "Quadrupole coupling matrices" tags=[:unit, :fast] begin
    using WignerSymbols: clebschgordan

    basis = StateBasis(["S_1/2", "D_5/2"])
    n, ε = beam_vectors(0.3, 0.7, 0.4)
    Γ = quadrupole_geometry(ε, n)
    C = quadrupole_couplings(basis, "S_1/2", "D_5/2", ε, n)

    # Entries only in the upper⟨row|lower⟩⟨col| block.
    s_range = staterange(basis, "S_1/2")
    d_range = staterange(basis, "D_5/2")
    @test all(iszero, C[s_range, :])
    @test all(iszero, C[:, d_range])

    # Stretch transitions carry unit Clebsch–Gordan factors.
    up = C[stateindex(basis, "D_5/2", 5//2), stateindex(basis, "S_1/2", 1//2)]
    @test up ≈ Γ[2+3]
    down = C[stateindex(basis, "D_5/2", -5//2), stateindex(basis, "S_1/2", -1//2)]
    @test down ≈ Γ[-2+3]

    # A non-stretch entry against the explicit Clebsch–Gordan factor.
    c = C[stateindex(basis, "D_5/2", 1//2), stateindex(basis, "S_1/2", -1//2)]
    @test c ≈ clebschgordan(1//2, -1//2, 2, 1, 5//2, 1//2) * Γ[1+3]

    # |Δm| > 2 entries vanish.
    @test iszero(C[stateindex(basis, "D_5/2", -5//2), stateindex(basis, "S_1/2", 1//2)])
end

@testitem "Circular polarisation channel assignment" tags=[:unit, :fast] begin
    # A σ⁺ beam along the quantisation axis (ε = (x̂ + iŷ)/√2, photon spin +ħ
    # along ẑ) must drive Δm = +1 by angular momentum conservation, for E1 and
    # E2 alike (a plane wave along ẑ carries no orbital angular momentum about
    # ẑ, so the extra rank-2 index is taken up by the ν = 0 gradient).
    n, ε = beam_vectors(0.0, π / 4, π / 2)
    @test ε ≈ [1, im, 0] / sqrt(2)

    d = dipole_geometry(ε)
    @test abs(d[1+2]) ≈ 1
    @test all(abs(d[q+2]) < 1e-14 for q in (-1, 0))

    Γ = quadrupole_geometry(ε, n)
    @test abs(Γ[1+3]) ≈ 1 / sqrt(2)
    @test all(abs(Γ[q+3]) < 1e-14 for q in (-2, -1, 0, 2))
end

@testitem "Rank-2 Clebsch–Gordan table" tags=[:unit, :fast] begin
    using WignerSymbols: clebschgordan

    # The tabulated ⟨1 μ; 1 ν | 2 (μ + ν)⟩ must match WignerSymbols exactly.
    for μ in -1:1, ν in -1:1
        @test Levels.RANK2_CG[μ+2, ν+2] ≈ clebschgordan(1, μ, 1, ν, 2, μ + ν)
    end
end

@testitem "Quadrupole geometry sum rule" tags=[:unit, :fast] begin
    # Σ_q |Γ_q|² = 1/2 for any transverse beam geometry: of the unit total
    # norm of ε ⊗ n, the rank-0 part vanishes (ε ⊥ n) and the rank-1 part
    # carries |ε × n|²/2 = 1/2.
    for (φ, γ, η) in
        ((0.0, 0.0, 0.0), (0.3, 0.7, 0.4), (1.1, π / 4, π / 2), (2.0, 1.2, -0.8))
        n, ε = beam_vectors(φ, γ, η)
        @test sum(abs2, quadrupole_geometry(ε, n)) ≈ 0.5
    end
end

@testitem "Hyperfine quadrupole couplings" tags=[:unit, :fast] begin
    using Unitful
    using Levels: fine_structure

    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    n, ε = beam_vectors(0.3, 0.7, 0.2)
    Γ = quadrupole_geometry(ε, n)
    C = quadrupole_couplings(ca43, basis, "S_1/2", "D_5/2", ε, n)

    # Element-wise agreement with the zero-field transition amplitudes in the
    # upper⟨row|lower⟩⟨col| block, zero elsewhere.
    s_manifold = NoHyperfineNumberSpec("S_1/2")
    d_manifold = NoHyperfineNumberSpec("D_5/2")
    function expected_coupling(ls, us)
        if fine_structure(ls.level) == s_manifold &&
           fine_structure(us.level) == d_manifold
            q = us.m - ls.m
            abs(q) <= 2 ? transition_amplitude(ca43, ls, us) * Γ[Int(q)+3] : 0.0im
        else
            0.0im
        end
    end
    @test maximum(
        abs(C[k, i] - expected_coupling(basis[i], basis[k])) for
        i in 1:length(basis), k in 1:length(basis)
    ) < 1e-13

    # The eigenbasis rotation of the zero-field matrix reproduces the at-field
    # transition amplitudes component by component.
    B = 0.5u"mT"
    v = eigenbasis_transform(ca43, basis, B)
    C_B = v' * C * v
    m_s = hyperfine_manifold(ca43, "S_1/2", B)
    m_d = hyperfine_manifold(ca43, "D_5/2", B)
    function expected_at_field(ls, us)
        if fine_structure(ls.level) == s_manifold &&
           fine_structure(us.level) == d_manifold
            q = us.m - ls.m
            abs(q) <= 2 ? transition_amplitude(m_s, m_d, ls, us) * Γ[Int(q)+3] : 0.0im
        else
            0.0im
        end
    end
    @test maximum(
        abs(C_B[k, i] - expected_at_field(basis[i], basis[k])) for
        i in 1:length(basis), k in 1:length(basis)
    ) < 1e-13

    # Restricting to a single F pair picks the corresponding block of the full
    # matrix and leaves everything else zero.
    C46 = quadrupole_couplings(ca43, basis, "S_1/2 F=4", "D_5/2 F=6", ε, n)
    r4 = staterange(basis, "S_1/2 F=4")
    r6 = staterange(basis, "D_5/2 F=6")
    @test C46[r6, r4] == C[r6, r4]
    C46[r6, r4] .= 0
    @test iszero(C46)

    # The no-hyperfine species-first form is a passthrough.
    nh_basis = StateBasis("S_1/2", "D_5/2")
    @test quadrupole_couplings(sr88, nh_basis, "S_1/2", "D_5/2", ε, n) ==
          quadrupole_couplings(nh_basis, "S_1/2", "D_5/2", ε, n)
end

# Vector spherical harmonics Y^{(+1)}_{K,p}(k̂) for K = 1, 2, transcribed from
# the explicit linear-polarisation-basis forms in Appendix A of [Campbell2026]
# (W. C. Campbell, "Angular Geometry of Atomic Multipole Transitions",
# arXiv:2510.07451): an independent formulation of the beam-geometry factors,
# in which the coupling of a plane wave to a 2^K-pole Δm channel is the
# bilinear projection ε · Y^{(+1)}_{K,-Δm}(k̂) of its polarisation onto the
# channel's normalised far-field emission pattern. Shared with the Rabi
# cross-check in test-rates.jl.
@testsnippet CampbellVSH begin
    θhat(ϑ, φ) = [cos(ϑ) * cos(φ), cos(ϑ) * sin(φ), -sin(ϑ)]
    φhat(ϑ, φ) = [-sin(φ), cos(φ), 0.0]
    khat(ϑ, φ) = [sin(ϑ) * cos(φ), sin(ϑ) * sin(φ), cos(ϑ)]

    function vsh(K, p, ϑ, φ)
        θv = complex.(θhat(ϑ, φ))
        φv = complex.(φhat(ϑ, φ))
        if K == 1 && p == -1
            cis(-φ) * sqrt(3 / (16π)) .* (cos(ϑ) .* θv .- im .* φv)
        elseif K == 1 && p == 0
            -sqrt(3 / (8π)) * sin(ϑ) .* θv
        elseif K == 1 && p == 1
            -cis(φ) * sqrt(3 / (16π)) .* (cos(ϑ) .* θv .+ im .* φv)
        elseif K == 2 && p == -2
            cis(-2φ) * sqrt(5 / (16π)) * sin(ϑ) .* (cos(ϑ) .* θv .- im .* φv)
        elseif K == 2 && p == -1
            cis(-φ) * sqrt(5 / (16π)) .* (cos(2ϑ) .* θv .- im * cos(ϑ) .* φv)
        elseif K == 2 && p == 0
            -sqrt(15 / (32π)) * sin(2ϑ) .* θv
        elseif K == 2 && p == 1
            -cis(φ) * sqrt(5 / (16π)) .* (cos(2ϑ) .* θv .+ im * cos(ϑ) .* φv)
        elseif K == 2 && p == 2
            cis(2φ) * sqrt(5 / (16π)) * sin(ϑ) .* (cos(ϑ) .* θv .+ im .* φv)
        else
            error("Y^(+1)_{$K,$p} not tabulated")
        end
    end

    # The paper's Eq. (8) dot product is bilinear (no conjugation) in the
    # e^{-iωt} rotating polarisation, matching the field convention of
    # beam_vectors.
    bilinear(a, b) = sum(a .* b)

    # Beam directions and transverse polarisation coefficients (ε = a θ̂ + b φ̂:
    # linear along either axis, circular, and generic elliptical) for the
    # cross-check grids.
    const CAMPBELL_ANGLES = ((0.4, 0.0), (1.0, 0.9), (π / 2, 2.0), (2.2, 4.5))
    const CAMPBELL_POLS = (
        (1.0 + 0.0im, 0.0 + 0.0im),
        (0.0 + 0.0im, 1.0 + 0.0im),
        (sqrt(0.5) + 0.0im, sqrt(0.5) * im),
        (0.6 + 0.3im, -0.4 + 0.62im),
    )
end

@testitem "Geometry factors vs vector spherical harmonics" tags=[:unit, :fast] setup=[
    CampbellVSH,
] begin
    # d_q = (-1)^q √(8π/3) ε·Y^{(+1)}_{1,-q}(k̂) and Γ_q = (-1)^q √(4π/5)
    # ε·Y^{(+1)}_{2,-q}(k̂), channel phases included: the Clebsch–Gordan
    # tensor-product construction coincides with [Campbell2026]'s Wigner-D
    # tables, and the constants tie the [James1998] normalisation used by
    # rabi_frequency to the paper's Eq. (8) (checked in test-rates.jl).
    for (ϑ, φ) in CAMPBELL_ANGLES, (a, b) in CAMPBELL_POLS
        ε = a .* θhat(ϑ, φ) .+ b .* φhat(ϑ, φ)
        d = dipole_geometry(ε)
        for q in -1:1
            @test d[q+2] ≈ (-1)^q * sqrt(8π / 3) * bilinear(ε, vsh(1, -q, ϑ, φ)) atol =
                1e-13
        end
        Γ = quadrupole_geometry(ε, khat(ϑ, φ))
        for q in -2:2
            @test Γ[q+3] ≈ (-1)^q * sqrt(4π / 5) * bilinear(ε, vsh(2, -q, ϑ, φ)) atol =
                1e-13
        end
    end

    # Summed over two orthogonal transverse polarisations, the channel weights
    # are the far-field emission power patterns W_{K,|q|} = |Y^{(+1)}_{K,q}|²
    # of the Δm = q components ([Campbell2026] Eq. (A13)): the same geometry
    # factors govern absorption from a beam and emission into its direction.
    for (ϑ, φ) in CAMPBELL_ANGLES
        pols = (complex.(θhat(ϑ, φ)), complex.(φhat(ϑ, φ)))
        for q in -1:1
            w = sum(abs2(dipole_geometry(ε)[q+2]) for ε in pols)
            @test w ≈ 8π / 3 * sum(abs2, vsh(1, q, ϑ, φ)) rtol = 1e-12 atol = 1e-14
        end
        for q in -2:2
            w = sum(abs2(quadrupole_geometry(ε, khat(ϑ, φ))[q+3]) for ε in pols)
            @test w ≈ 4π / 5 * sum(abs2, vsh(2, q, ϑ, φ)) rtol = 1e-12 atol = 1e-14
        end
    end
end

@testitem "Two-beam E2 interference" tags=[:unit, :fast] begin
    # For E2, a single plane wave cannot isolate one |Δm| = 2 component: k ⊥ ẑ
    # with in-plane (φ̂) linear polarisation maximises the per-power |Δm| = 2
    # coupling but splits it evenly between Δm = ±2 ([Campbell2026] §IV A 2).
    φ_pol(φ_k) = complex.([-sin(φ_k), cos(φ_k), 0.0])
    k_xy(φ_k) = [cos(φ_k), sin(φ_k), 0.0]
    Γ1 = quadrupole_geometry(φ_pol(0.0), k_xy(0.0))
    @test abs(Γ1[-2+3]) ≈ 0.5
    @test abs(Γ1[2+3]) ≈ 0.5
    @test all(abs(Γ1[q+3]) < 1e-15 for q in -1:1)

    # A second phase-coherent beam of the same geometry at azimuth Δφ_k, with
    # drive phase π - 2Δφ_k relative to the first, cancels Δm = -2 exactly
    # while Δm = +2 adds fully ([Campbell2026] §IV B) — the complex channel
    # amplitudes support such coherent multi-beam geometries by summation.
    Δφ_k = π / 4
    Γ2 = quadrupole_geometry(φ_pol(Δφ_k), k_xy(Δφ_k))
    Γ = Γ1 .+ cis(π - 2Δφ_k) .* Γ2
    @test abs(Γ[-2+3]) < 1e-15
    @test abs(Γ[2+3]) ≈ 1.0
    @test all(abs(Γ[q+3]) < 1e-15 for q in -1:1)
end
