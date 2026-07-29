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
