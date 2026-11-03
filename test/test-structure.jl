# Tests for the Eigenstructure: construction and validation, the accessor
# interface, the basis identities of the "Conventions" page, the
# zero-field limit and eigenvector gauge, partial and reordered bases, the
# adiabatic labelling across hyperfine-constant signs and I/J orderings, and the
# energy references.

@testitem "Eigenstructure construction and accessors" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful

    B = 0.5u"mT"
    s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B)
    basis = StateBasis(ca43, "S_1/2", "D_5/2")
    @test species(s) === ca43
    @test field(s) == B
    @test collect(Levels.basis(s)) == collect(basis)
    @test collect(coupled_basis(s)) == collect(basis)
    @test length(offsets(s)) == 64
    @test unit(eltype(offsets(s))) == u"rad/µs"
    @test occursin("Eigenstructure(S_1/2, D_5/2; 64 states; B = 0.5 mT)", repr(s))

    # Identity 1: orthonormal columns; square and orthogonal for a complete
    # canonical basis.
    W = eigenvectors(s)
    @test size(W) == (64, 64)
    @test W' * W ≈ I atol = 1e-12
    @test W * W' ≈ I atol = 1e-12
    # Columns live in their own manifold.
    s_range, d_range = staterange(basis, "S_1/2"), staterange(basis, "D_5/2")
    @test all(iszero, W[s_range, d_range]) && all(iszero, W[d_range, s_range])

    # Identity 2: the offset Hamiltonian rotates onto the labelled offsets.
    H = coupled_hamiltonian(s)
    @test unit(eltype(H)) == u"rad/µs"
    rotated = ustrip.(u"rad/µs", to_eigenbasis(s, H))
    @test rotated ≈ diagm(ustrip.(u"rad/µs", offsets(s))) atol = 1e-9

    # Identity 7: the product-basis moment operator conjugated with the
    # Clebsch–Gordan unitary gives the coupled moments, and together with the
    # Casimir shifts the manifold block of the offset Hamiltonian.
    for level in ("S_1/2", "D_5/2")
        U = Levels.coupling_transform(ca43, level)
        M = Levels.moment_operators(ca43, level)
        @test U' * M.z * U ≈ Levels.coupled_moments(ca43, level).z
        r = staterange(coupled_basis(s), level)
        casimir = [Levels.hyperfine_shift(ca43, st.level) for st in collect(basis)[r]]
        @test H[r, r] ≈ Diagonal(casimir) .+ B .* (U' * M.z * U)
    end

    # The constructor forms agree.
    @test eigenvectors(Eigenstructure(ca43, basis; B)) == W
    @test eigenvectors(Eigenstructure(ca43, collect(basis); B)) == W
    @test eigenvectors(Eigenstructure(ca43, ["S_1/2", "D_5/2"]; B)) == W
    @test offsets(Eigenstructure(ca43, "S_1/2 F=4"; B)) ==
          offsets(s)[staterange(basis, "S_1/2 F=4")]

    # Validation: basis kind, species levels, hyperfine constants, the field.
    @test_throws ArgumentError Eigenstructure(ca43, StateBasis(["S_1/2", "D_5/2"]); B)
    @test_throws ArgumentError Eigenstructure(sr88, StateBasis(ca43, "S_1/2"); B)
    @test_throws ArgumentError Eigenstructure(sr88, "G_7/2"; B)
    @test_throws ArgumentError Eigenstructure(ca43, "S_1/2 F=5"; B)
    @test_throws UndefKeywordError Eigenstructure(sr88, "S_1/2")
    for bad in ([0.0, 0.0, 0.5]u"mT", 1u"s", NaN * u"mT", Inf * u"mT", (1 + im)u"mT")
        @test_throws ArgumentError Eigenstructure(sr88, "S_1/2"; B=bad)
    end
    no_constants = HyperfineOneElectronSpecies(;
        mass=ca43.mass,
        nuclear_spin=ca43.nuclear_spin,
        nuclear_g=ca43.nuclear_g,
        energies=ca43.energies,
        einstein_as=ca43.einstein_as,
        hyperfine=Dict(
            convert(NoHyperfineNumberSpec, "S_1/2") =>
                ca43.hyperfine[convert(NoHyperfineNumberSpec, "S_1/2")],
        ),
    )
    @test Eigenstructure(no_constants, "S_1/2"; B) isa Eigenstructure
    @test_throws ArgumentError Eigenstructure(no_constants, "S_1/2", "D_5/2"; B)
    @test_throws ArgumentError Eigenstructure(no_constants, "D_5/2"; B=0.0u"mT")

    # A manifold outside the basis solves afresh; one inside is the stored one.
    @test Levels.manifold_solution(s, "S_1/2 F=4") ===
          s.manifolds[convert(NoHyperfineNumberSpec, "S_1/2")]
    p = Levels.manifold_solution(s, "P_1/2")
    @test p isa Levels.ManifoldSolution && length(p.basis) == 16
    @test_throws ArgumentError Levels.manifold_solution(
        Eigenstructure(no_constants, "S_1/2"; B),
        "D_5/2",
    )

    # The static field is stored as a Float64 flux density in mT, whatever it
    # was given as.
    @test field(Eigenstructure(sr88, "S_1/2"; B=5u"Gauss")) === 0.5u"mT"
end

@testitem "Fine-structure species: identity eigenvectors, field-dependent energies" tags=[
    :unit,
    :fast,
] begin
    using LinearAlgebra
    using Unitful

    μ_B = Levels.BOHR_MAGNETON
    for species in (sr88, ca40), B in (0.5u"mT", -0.3u"mT", 0.0u"mT")
        s = Eigenstructure(species, "S_1/2", "D_5/2"; B)
        @test eigenvectors(s) == I
        for (i, state) in enumerate(s.basis)
            expected = uconvert(
                u"rad/µs",
                lande_g(species, state.level) * state.m * μ_B * B / u"ħ",
            )
            @test offsets(s)[i] ≈ expected
            @test zeeman_shift(s, state) ≈ expected
            @test state_offset(s, state) == offsets(s)[i]
        end
        @test coupled_hamiltonian(s) ≈ Diagonal(offsets(s))
        @test iszero(coupled_hamiltonian(s) - Diagonal(coupled_hamiltonian(s)))
        # Species and structure forms of the field-independent angular
        # quantities agree in matching orders.
        n, ε = beam_vectors(0.3, 0.7, 0.4)
        @test relative_couplings(s, "S_1/2" => "D_5/2", ε, n) ==
              relative_couplings(species, s.basis, "S_1/2" => "D_5/2", ε, n)
        @test coupling_matrix(s, "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n) ≈
              coupling_matrix(species, s.basis, "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n)
        @test zeeman_hamiltonian(s, [0.2, -0.1, 0.3]u"mT") ≈
              zeeman_hamiltonian(species, s.basis, [0.2, -0.1, 0.3]u"mT")
        for pair in state_pairs("S_1/2" => "D_5/2"; Δm=-2:2)
            @test transition_amplitude(s, pair) == transition_amplitude(species, pair)
            @test zeeman_sensitivity(s, pair) ≈ zeeman_sensitivity(species, pair)
        end
    end
    # A scrambled or partial fine-structure basis selects columns of the identity.
    scrambled = Eigenstructure(
        sr88,
        [StateSpec("D_5/2", -1//2), StateSpec("S_1/2", 1//2), StateSpec("D_5/2", 5//2)];
        B=0.5u"mT",
    )
    W = eigenvectors(scrambled)
    @test size(W) == (8, 3)
    @test W' * W ≈ I
    for (k, state) in enumerate(scrambled.basis)
        @test W[:, k] == [c == state for c in coupled_basis(scrambled)]
    end
    # The transition frequency and the polarisability denominators do depend
    # on the field.
    clock = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -1//2)
    @test transition_frequency(
        Eigenstructure(sr88, "S_1/2", "D_5/2"; B=0.5u"mT"),
        clock,
    ) != transition_frequency(
        Eigenstructure(sr88, "S_1/2", "D_5/2"; B=0.0u"mT"),
        clock,
    )
end

@testitem "Zero field, low-field continuity and eigenvector gauge" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful

    # Identity 4: at zero field the eigenvectors are the identity and the
    # offsets the Casimir shifts; the structure accepts B = 0 as an ordinary
    # field.
    s0 = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=0.0u"mT")
    @test eigenvectors(s0) == I
    @test offsets(s0) == [Levels.hyperfine_shift(ca43, st.level) for st in s0.basis]
    @test all(iszero(zeeman_shift(s0, st)) for st in s0.basis)
    sp, dp = StateSpec("S_1/2 F=4", 4), StateSpec("D_5/2 F=4", 3)
    @test transition_amplitude(s0, sp => dp) ≈ transition_amplitude(ca43, sp => dp)
    @test transition_frequency(s0, sp => dp) ≈
          Levels.transition_frequency(ca43, sp.level => dp.level)

    # On the resolved ⁴³Ca⁺ spectrum, W → I linearly as B → 0.
    deviation(b) =
        maximum(abs, eigenvectors(Eigenstructure(ca43, "S_1/2", "D_5/2"; B=b)) - I)
    @test deviation(1e-6u"mT") < 1e-5
    @test deviation(1e-3u"mT") / deviation(1e-6u"mT") ≈ 1e3 rtol = 0.05
    @test deviation(1e-9u"mT") / deviation(1e-6u"mT") ≈ 1e-3 rtol = 0.05
    # Low-field continuity of the energies with the Casimir values.
    for b in (1e-9u"mT", 1e-6u"mT")
        sb = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=b)
        linear = [
            lande_g(ca43, st.level) * st.m * Levels.BOHR_MAGNETON * b / u"ħ" for
            st in sb.basis
        ]
        @test offsets(sb) ≈ offsets(s0) .+ linear rtol = 1e-3 atol = 1e-6u"rad/µs"
        # m_F is exact: every column is supported in its own m block.
        W = eigenvectors(sb)
        for (k, st) in enumerate(sb.basis), (i, c) in enumerate(coupled_basis(sb))
            c.m == st.m || @test iszero(W[i, k])
        end
    end

    # Own-label-positive gauge: the diagonal stays positive over the whole
    # working-field range up to the first gauge crossing (1.29 mT in D_5/2).
    for b in exp10.(range(-3, log10(1.2); length=60)) .* u"mT"
        @test all(>(0), diag(eigenvectors(Eigenstructure(ca43, "S_1/2", "D_5/2"; B=b))))
    end
    # At the documented crossing the fallback keeps every column normalised and
    # the rotation identity intact.
    sx = Eigenstructure(ca43, "D_5/2"; B=1.29u"mT")
    @test eigenvectors(sx)' * eigenvectors(sx) ≈ I atol = 1e-12
    @test ustrip.(u"rad/µs", to_eigenbasis(sx, coupled_hamiltonian(sx))) ≈
          diagm(ustrip.(u"rad/µs", offsets(sx))) atol = 1e-9

    # A whole-column rephasing is a gauge change: consistently transformed
    # operators give the same observables.
    s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=0.5u"mT")
    n, ε = beam_vectors(π / 2, π / 4, 0.3)
    C = ustrip.(u"rad/µs", coupling_matrix(s, "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n))
    D = Diagonal(rand((-1.0, 1.0), 64))
    C_gauge = D' * C * D
    @test abs.(C_gauge) ≈ abs.(C)
    @test eigvals(Hermitian((C_gauge + C_gauge') / 2)) ≈
          eigvals(Hermitian((C + C') / 2))
    @test any(!≈(0), C_gauge .- C) || D == I

    # Zero-constants counterexample: degenerate F branches do not approach the
    # coupled basis; the structure still builds, with exact m, deterministic
    # (ascending-F) labels and the own-label component positive.
    toy = HyperfineOneElectronSpecies(;
        mass=sr88.mass,
        nuclear_spin=1//2,
        nuclear_g=-2e-4 * uconvert(NoUnits, Unitful.mp / Unitful.me),
        energies=sr88.energies,
        lande_g_overrides=sr88.lande_g_overrides,
        hyperfine=Dict(
            convert(NoHyperfineNumberSpec, "S_1/2") =>
                HyperfineConstants(; a=0.0u"J"),
        ),
        einstein_as=sr88.einstein_as,
    )
    for b in (1e-6u"mT", 1e-9u"mT", 1e-12u"mT")
        st = Eigenstructure(toy, "S_1/2"; B=b)
        W = eigenvectors(st)
        @test maximum(abs, W - I) ≈ sqrt(0.5) atol = 1e-6
        @test all(>(0), diag(W))
        @test W' * W ≈ I atol = 1e-12
        @test state_offset(st, "S_1/2 F=0", 0) < state_offset(st, "S_1/2 F=1", 0)
    end
    @test eigenvectors(Eigenstructure(toy, "S_1/2"; B=0.0u"mT")) == I
end

@testitem "Adiabatic labels across constants, spins and field signs" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful

    # Toy species covering positive and negative A, a quadrupole constant, I > J
    # and I < J, both field signs: the block solve must give exact m, a
    # low-field limit on the zero-field labels, and energies that are exactly
    # the eigenvalues of the manifold Hamiltonian.
    function toy(i_nuc, a, b)
        HyperfineOneElectronSpecies(;
            mass=ca43.mass,
            nuclear_spin=i_nuc,
            nuclear_g=ca43.nuclear_g,
            energies=ca43.energies,
            einstein_as=ca43.einstein_as,
            hyperfine=Dict(
                convert(NoHyperfineNumberSpec, "S_1/2") => HyperfineConstants(; a),
                convert(NoHyperfineNumberSpec, "D_5/2") => HyperfineConstants(; a, b),
            ),
        )
    end
    for i_nuc in (3//2, 7//2),
        a in (u"h" * 100.0u"MHz", u"h" * -100.0u"MHz"),
        b in (0.0u"J", u"h" * 5.0u"MHz"),
        B in (0.3u"mT", -0.3u"mT")

        sp = toy(i_nuc, a, b)
        s = Eigenstructure(sp, "S_1/2", "D_5/2"; B)
        W = eigenvectors(s)
        @test W' * W ≈ I atol = 1e-12
        for (k, st) in enumerate(s.basis), (i, c) in enumerate(coupled_basis(s))
            c.m == st.m || @test iszero(W[i, k])
        end
        @test all(>(0), diag(W))
        H = ustrip.(u"rad/µs", coupled_hamiltonian(s))
        for level in ("S_1/2", "D_5/2")
            r = staterange(coupled_basis(s), level)
            @test sort(eigvals(Symmetric(H[r, r]))) ≈
                  sort(ustrip.(u"rad/µs", offsets(s)[r]))
        end
        # Labels follow the zero-field order into the field: at a tiny field
        # the eigenvectors are the coupled states.
        small = Eigenstructure(sp, "S_1/2", "D_5/2"; B=sign(B) * 1e-7u"mT")
        @test maximum(abs, eigenvectors(small) - I) < 1e-3
        # Sign of the field: the mirror state at the mirrored field.
        mirror = Eigenstructure(sp, "S_1/2", "D_5/2"; B=(-B))
        for st in s.basis
            @test state_offset(mirror, StateSpec(st.level, -st.m)) ≈ state_offset(s, st) rtol =
                1e-12 atol = 1e-9u"rad/µs"
        end
    end

    # Small-field robustness on every shipped species: structures at 1 pT, 1 nT
    # and 1 µT build, with exact m and offsets continuous with zero field.
    for sp in (sr88, ca40, ca43)
        levels_with_data = [
            l for l in levels(sp) if
            sp isa NoHyperfineOneElectronSpecies || haskey(sp.hyperfine, l)
        ]
        s0 = Eigenstructure(sp, levels_with_data; B=0.0u"mT")
        for b in (1e-9u"mT", 1e-6u"mT", 1e-3u"mT")
            s = Eigenstructure(sp, levels_with_data; B=b)
            W = eigenvectors(s)
            @test W' * W ≈ I atol = 1e-12
            for (k, st) in enumerate(s.basis), (i, c) in enumerate(coupled_basis(s))
                c.m == st.m || @test iszero(W[i, k])
            end
            scale = uconvert(u"rad/µs", 10 * Levels.BOHR_MAGNETON * b / u"ħ")
            @test maximum(abs, offsets(s) - offsets(s0)) < scale
        end
    end
end

@testitem "Partial and reordered bases" tags=[:unit, :fast] begin
    using LinearAlgebra
    using Unitful

    B = 0.5u"mT"
    full = Eigenstructure(ca43, "S_1/2", "D_5/2"; B)
    states = [
        StateSpec("D_5/2 F=4", 3),
        StateSpec("S_1/2 F=4", 4),
        StateSpec("D_5/2 F=3", -2),
    ]
    partial = Eigenstructure(ca43, states; B)
    @test collect(Levels.basis(partial)) == states
    # The coupled basis is complete, manifolds in order of first appearance.
    @test length(coupled_basis(partial)) == 64
    @test staterange(coupled_basis(partial), "D_5/2") == 1:48
    @test staterange(coupled_basis(partial), "S_1/2") == 49:64
    @test !iscomplete(partial, "S_1/2") && !iscomplete(partial, "D_5/2 F=4")
    @test iscomplete(full, "S_1/2") && iscomplete(full, "D_5/2 F=4")
    @test !iscomplete(full, "P_1/2")

    # Columns are those of the complete solution (up to the coupled-basis order).
    Wp = eigenvectors(partial)
    @test size(Wp) == (64, 3)
    @test Wp' * Wp ≈ I atol = 1e-12
    Wf = eigenvectors(full)
    for (k, st) in enumerate(states), c in coupled_basis(partial)
        @test Wp[stateindex(coupled_basis(partial), c), k] ==
              Wf[stateindex(coupled_basis(full), c), stateindex(Levels.basis(full), st)]
    end
    @test offsets(partial) == [state_offset(full, st) for st in states]

    # Projection: to_eigenbasis over the subset is the submatrix of the complete
    # rotation; a matrix over the subset itself is refused.
    b = [0.2, -0.1, 0.3]u"mT"
    sub = to_eigenbasis(partial, zeeman_hamiltonian(ca43, coupled_basis(partial), b))
    whole = to_eigenbasis(full, zeeman_hamiltonian(ca43, coupled_basis(full), b))
    for (i, si) in enumerate(states), (k, sk) in enumerate(states)
        @test sub[i, k] ≈ whole[
            stateindex(Levels.basis(full), si),
            stateindex(Levels.basis(full), sk),
        ] atol = 1e-12u"rad/µs"
    end
    @test zeeman_hamiltonian(partial, b) ≈ sub
    @test_throws ArgumentError to_eigenbasis(partial, zeros(3, 3))
    @test_throws ArgumentError to_eigenbasis(full, zeros(16, 16))

    # Matrix factories address the retained eigenstates with the exact
    # complete-manifold elements.
    n, ε = beam_vectors(π / 2, π / 4, 0.3)
    Cp = coupling_matrix(partial, "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n)
    Cf = coupling_matrix(full, "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n)
    for (i, si) in enumerate(states), (k, sk) in enumerate(states)
        @test Cp[k, i] ==
              Cf[stateindex(Levels.basis(full), sk), stateindex(Levels.basis(full), si)]
    end
    Rp = relative_couplings(partial, "S_1/2" => "D_5/2", ε, n)
    @test Rp[1, 2] ≈
          transition_amplitude(full, states[2] => states[1]) *
          quadrupole_geometry(ε, n)[Int(states[1].m-states[2].m)+3]
    @test iszero(Rp[2, 1]) && iszero(Rp[3, 1])
    # Scalar queries accept any state of the coupled basis.
    outside = StateSpec("S_1/2 F=3", 0)
    @test state_offset(partial, outside) == state_offset(full, outside)
    @test state_energy(partial, outside) == state_energy(full, outside)
    @test zeeman_shift(partial, outside) == zeeman_shift(full, outside)
    @test transition_amplitude(partial, outside => StateSpec("D_5/2 F=3", 1)) ==
          transition_amplitude(full, outside => StateSpec("D_5/2 F=3", 1))
    @test_throws ArgumentError state_offset(partial, StateSpec("P_1/2 F=4", 0))
    @test_throws ArgumentError state_offset(partial, StateSpec("S_1/2", 1//2))
    @test_throws ArgumentError state_offset(partial, StateSpec("S_1/2 F=5", 0))

    # Projection does not preserve products (rectangular W).
    X = Levels.coupled_moments(ca43, "D_5/2").x
    Xe = eigenvectors(partial)[1:48, [1, 3]]' * X * eigenvectors(partial)[1:48, [1, 3]]
    full_sq =
        eigenvectors(partial)[1:48, [1, 3]]' *
        (X * X) *
        eigenvectors(partial)[1:48, [1, 3]]
    @test !isapprox(Xe * Xe, full_sq; rtol=1e-6)
end

@testitem "Energy references and intervals" tags=[:unit, :fast] begin
    using Unitful

    B = 0.5u"mT"
    s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B)
    sp, dp = StateSpec("S_1/2 F=4", 4), StateSpec("D_5/2 F=4", 3)

    # Total energies include the centroid in the species' reference; offsets do
    # not.
    @test state_energy(s, sp) ≈ level_energy(ca43, "S_1/2") + u"ħ" * state_offset(s, sp)
    @test state_energy(s, dp) ≈ level_energy(ca43, "D_5/2") + u"ħ" * state_offset(s, dp)
    @test state_energy(s, dp.level, dp.m) == state_energy(s, dp)
    @test state_offset(s, dp.level, dp.m) == state_offset(s, dp)
    @test unit(state_energy(s, sp)) == u"J"
    @test abs(state_offset(s, dp)) < angular(100u"MHz")

    # The interval API separates centroid and offsets.
    @test transition_frequency(s, sp => dp) ≈
          Levels.transition_frequency(ca43, "S_1/2 F=4" => "D_5/2 F=4") +
          zeeman_shift(s, dp) - zeeman_shift(s, sp)
    @test transition_frequency(s, "S_1/2" => "D_5/2") ==
          Levels.transition_frequency(ca43, "S_1/2" => "D_5/2")
    # A same-manifold microwave interval needs no optical rate entry, and a
    # structure over the one manifold suffices.
    s_only = Eigenstructure(ca43, "S_1/2"; B)
    mw = StateSpec("S_1/2 F=4", 0) => StateSpec("S_1/2 F=3", 1)
    @test transition_frequency(s_only, mw) ==
          state_offset(s_only, mw.second) - state_offset(s_only, mw.first)
    @test transition_frequency(s_only, mw) ≈ transition_frequency(s, mw)
    @test transition_frequency(
        Eigenstructure(ca43, "S_1/2"; B=1e-9u"mT"),
        StateSpec("S_1/2 F=4", 0) => StateSpec("S_1/2 F=3", 0),
    ) ≈ angular(3225.60828640u"MHz") rtol = 1e-9
    @test_throws ArgumentError transition_frequency(s_only, sp => dp)

    # Zeeman shift and sensitivity are exact structure quantities.
    @test zeeman_shift(s, sp) ≈
          state_offset(s, sp) - Levels.hyperfine_shift(ca43, sp.level)
    @test zeeman_sensitivity(s, sp => dp) ≈ zeeman_sensitivity(s, [sp => dp][1])
    @test unit(zeeman_sensitivity(s, sp => dp)) == u"rad/µs/mT"
end
