@testitem "Angular momentum matrices" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Levels: jx_matrix, jy_matrix, jz_matrix, jplus_matrix

    # Pauli matrices over two for j = 1/2; note the increasing-m basis order flips
    # the sign of σ_y and σ_z relative to the usual (spin-up, spin-down) form.
    @test jx_matrix(1//2) ≈ [0 1; 1 0] / 2
    @test jy_matrix(1//2) ≈ [0 im; -im 0] / 2
    @test jz_matrix(1//2) ≈ [-1 0; 0 1] / 2

    # ⟨m + 1|J₊|m⟩ = √(j(j+1) − m(m+1)), e.g. ⟨−3/2|J₊|−5/2⟩ = √5.
    @test jplus_matrix(5//2)[2, 1] ≈ sqrt(5)
    @test jplus_matrix(1//2) == [0 0; 1 0]

    for j in (1//2, 1, 3//2, 5//2)
        jx, jy, jz = jx_matrix(j), jy_matrix(j), jz_matrix(j)
        # Commutation relations [Jx, Jy] = iJz (and cyclic).
        @test jx * jy - jy * jx ≈ im * Matrix(jz)
        @test jy * jz - jz * jy ≈ im * jx
        @test jz * jx - jx * jz ≈ im * jy
        # Casimir invariant J² = j(j+1)·1.
        @test jx^2 + jy^2 + Matrix(jz)^2 ≈ float(j * (j + 1)) * I
    end
end

@testitem "Spherical vector components" tags=[:unit, :fast] begin
    using Levels: spherical_component

    @test spherical_component([0, 0, 1], 0) ≈ 1
    @test spherical_component([0, 0, 1], 1) ≈ 0 atol = 1e-15
    @test spherical_component([1, 0, 0], 1) ≈ -1 / sqrt(2)
    @test spherical_component([1, 0, 0], -1) ≈ 1 / sqrt(2)
    @test spherical_component([0, 1, 0], 1) ≈ -im / sqrt(2)
    @test spherical_component([0, 1, 0], -1) ≈ -im / sqrt(2)

    # Norm is preserved: Σ_q |v_q|² = |v|².
    v = [0.3, -1.2, 0.7]
    @test sum(abs2(spherical_component(v, q)) for q in -1:1) ≈ sum(abs2, v)

    # Unitful vectors keep their units.
    using Unitful
    @test spherical_component([1.0, 0.0, 0.0]u"mT", -1) ≈ (1 / sqrt(2))u"mT"

    @test_throws ArgumentError spherical_component([1, 0, 0], 2)
end

@testitem "Clebsch–Gordan convention" tags=[:unit, :fast] begin
    # Levels relies on WignerSymbols.clebschgordan ⟨j₁ m₁; j₂ m₂|J M⟩ (Condon–
    # Shortley); anchor the convention with known values.
    using WignerSymbols: clebschgordan

    # Stretch states couple with unit amplitude.
    @test clebschgordan(1//2, 1//2, 2, 2, 5//2, 5//2) ≈ 1.0
    @test clebschgordan(1//2, -1//2, 2, -2, 5//2, -5//2) ≈ 1.0

    # Known integer-spin values.
    @test clebschgordan(1, 1, 1, -1, 2, 0) ≈ 1 / sqrt(6)
    @test clebschgordan(1, 0, 1, 0, 2, 0) ≈ sqrt(2 / 3)
    @test clebschgordan(1, 1, 1, 0, 2, 1) ≈ 1 / sqrt(2)

    # Completeness: Σ_J |⟨1/2 m; 2 q|J m+q⟩|² = 1. (clebschgordan throws a
    # DomainError for |m + q| > J rather than returning zero, so skip those.)
    for m in (-1//2, 1//2), q in -2:2
        Js = [J for J in (3//2, 5//2) if abs(m + q) <= J]
        @test sum(clebschgordan(1//2, m, 2, q, J, m + q)^2 for J in Js) ≈ 1.0
    end
end
