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
    using Levels: clebsch_gordan, fine_structure

    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    n, ε = beam_vectors(0.3, 0.7, 0.2)
    Γ = quadrupole_geometry(ε, n)
    C = quadrupole_couplings(ca43, basis, "S_1/2", "D_5/2", ε, n)

    # Element-wise agreement with the species-first Clebsch–Gordan amplitudes
    # in the upper⟨row|lower⟩⟨col| block, zero elsewhere.
    s_manifold = NoHyperfineNumberSpec("S_1/2")
    d_manifold = NoHyperfineNumberSpec("D_5/2")
    function expected_coupling(ls, us)
        if fine_structure(ls.level) == s_manifold &&
           fine_structure(us.level) == d_manifold
            q = us.m - ls.m
            abs(q) <= 2 ? clebsch_gordan(ca43, ls, us) * Γ[Int(q)+3] : 0.0im
        else
            0.0im
        end
    end
    @test maximum(
        abs(C[k, i] - expected_coupling(basis[i], basis[k])) for
        i in 1:length(basis), k in 1:length(basis)
    ) < 1e-14

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
