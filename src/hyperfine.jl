# Hyperfine-structure machinery: the coupled |F, m_F⟩ basis is the canonical
# matrix/index convention (hyperfine interaction diagonal, matrices indexable by
# StateBasis{HyperfineNumberSpec}); the |m_J, m_I⟩ product basis appears only as
# an internal construction device, glued by the Clebsch–Gordan unitary of
# coupling_transform. The at-field eigen-solution lives in structure.jl.

using LinearAlgebra
using Unitful
using WignerSymbols
using SciMLBase: IntervalNonlinearProblem, solve, successful_retcode
using SimpleNonlinearSolve: ITP

"""
    hyperfine_levels(species, fs_level)

Returns the hyperfine ``F`` levels of the given fine-structure level, in order
of increasing ``F`` from ``|I - J|`` to ``I + J``.
"""
function hyperfine_levels(species::HyperfineOneElectronSpecies, fs_level)
    spec = fine_structure(fs_level)
    i = species.nuclear_spin
    [HyperfineNumberSpec(spec.l, spec.j, f) for f in abs(i-spec.j):(i+spec.j)]
end

"""
Returns the hyperfine levels a level argument stands for: a single (validated)
``F`` level, or all ``F`` levels of a fine-structure manifold.
"""
function hyperfine_level_list(species::HyperfineOneElectronSpecies, entry)
    spec = parse_level(entry)
    if spec isa HyperfineNumberSpec
        [validate_hyperfine(species, spec)]
    else
        hyperfine_levels(species, spec)
    end
end

"""
    StateBasis(species::HyperfineOneElectronSpecies, levels...)

Creates a hyperfine state basis, expanding any fine-structure level into all its
``F`` levels (in order of increasing ``F``, each contributing its ``2F + 1``
sublevels in order of increasing ``m_F``).

This is the canonical coupled basis the manifold solutions of an
[`Eigenstructure`](@ref) are expressed in.
"""
function StateBasis(species::HyperfineOneElectronSpecies, levels...)
    StateBasis([l for entry in levels for l in hyperfine_level_list(species, entry)])
end

"""
    hyperfine_shift(species, level)

Returns the zero-field hyperfine shift of the given ``F`` level relative to the
centroid of its fine-structure level, in angular frequency units.

The shift is ``E/ħ``, where the Casimir energy for
``K = F(F+1) - I(I+1) - J(J+1)`` is

```math
E = \\frac{A K}{2}
  + B \\, \\frac{\\frac{3}{2} K (K + 1) - 2 I (I + 1) J (J + 1)}
               {4 I (2I - 1) J (2J - 1)}
```

(the quadrupole term only for ``I, J > 1/2``), with the constants from the
species' [`HyperfineConstants`](@ref).
"""
function hyperfine_shift(species::HyperfineOneElectronSpecies, level)
    spec = parse_level(level)
    if !(spec isa HyperfineNumberSpec)
        throw(ArgumentError("Level '$level' does not specify a hyperfine F level"))
    end
    validate_hyperfine(species, spec)
    fs = fine_structure(spec)
    consts = get(species.hyperfine, fs) do
        throw(ArgumentError("No hyperfine constants known for level '$fs'"))
    end

    i = species.nuclear_spin
    j = spec.j
    f = spec.f
    k = f * (f + 1) - i * (i + 1) - j * (j + 1)
    e = consts.a * k / 2
    if !iszero(consts.b)
        if i <= 1//2 || j <= 1//2
            throw(
                ArgumentError(
                    "Quadrupole hyperfine constant requires I, J > 1/2 for '$fs'",
                ),
            )
        end
        e +=
            consts.b * (3//2 * k * (k + 1) - 2 * i * (i + 1) * j * (j + 1)) /
            (4 * i * (2i - 1) * j * (2j - 1))
    end
    uconvert(u"rad/µs", e / u"ħ")
end

# A fine-structure species has no hyperfine structure: the shift of any of its
# levels from its own centroid is zero, which lets frequency bookkeeping (cf.
# `Levels.OpticalBloch`) treat both species kinds alike.
function hyperfine_shift(species::NoHyperfineOneElectronSpecies, level)
    validate_level_kind(species, level)
    zero(1.0u"rad/µs")
end

"""
    moment_operators(species, fs_level) -> (; x, y, z)

Returns the magnetic operator
``\\vec M = (g_J μ_B \\vec J - g_I μ_N \\vec I)/ħ`` of one fine-structure
manifold, in angular-frequency units per flux density, so that the Zeeman
Hamiltonian in angular units is ``\\vec B ⋅ \\vec M``. Here ``J`` and ``I``
are dimensionless angular-momentum matrices.

The product basis is `kron(J-space, I-space)`: electronic projection varies
slowest, nuclear projection fastest, both increasing. Use
[`Levels.coupling_transform`](@ref) to express an operator in the coupled basis.
"""
function moment_operators(species::HyperfineOneElectronSpecies, fs_level)
    spec = fine_structure(parse_level(fs_level))
    i = species.nuclear_spin
    j = spec.j
    g_j = lande_g(species, spec)
    g_i = species.nuclear_g

    eye_i = Matrix{Float64}(I, Int(2i + 1), Int(2i + 1))
    eye_j = Matrix{Float64}(I, Int(2j + 1), Int(2j + 1))
    component(op) = uconvert.(
        u"rad/µs/mT",
        (g_j * BOHR_MAGNETON / u"ħ") .* kron(op(j), eye_i) .-
        (g_i * NUCLEAR_MAGNETON / u"ħ") .* kron(eye_j, op(i)),
    )
    (x=component(jx_matrix), y=component(jy_matrix), z=component(jz_matrix))
end

"""
    coupling_transform(species, fs_level)

Returns the real Clebsch–Gordan matrix
``U_{(m_J,m_I),(F,m_F)} = ⟨J m_J; I m_I | F m_F⟩`` for one manifold.
Rows use the product order of [`Levels.moment_operators`](@ref); columns use
`StateBasis(species, fs_level)`. A product-basis operator transforms as
``X_c = U^† X_p U``. Nonzero entries satisfy ``m_J + m_I = m_F``.

Coupling is electron first, with Condon–Shortley phases. Nuclear-first states
differ by ``(-1)^{I+J-F}``; see [`Levels.hyperfine_reduction`](@ref) for the
amplitude convention. The product storage order fixes indices independently
of this phase choice.
"""
function coupling_transform(species::HyperfineOneElectronSpecies, fs_level)
    spec = fine_structure(parse_level(fs_level))
    i = species.nuclear_spin
    j = spec.j
    d_i = Int(2i + 1)
    num_states = Int(2j + 1) * d_i

    U = zeros(num_states, num_states)
    col = 0
    for f in abs(i-j):(i+j), m_f in (-f):f
        col += 1
        for (j_row, m_j) in enumerate((-j):j), (i_row, m_i) in enumerate((-i):i)
            m_i + m_j == m_f || continue
            U[(j_row-1)*d_i+i_row, col] = Float64(clebschgordan(j, m_j, i, m_i, f, m_f))
        end
    end
    U
end

"""
Returns the Cartesian components of the magnetic-moment operator of the given
fine-structure manifold in the canonical coupled basis, i.e.
[`Levels.moment_operators`](@ref) conjugated with
[`Levels.coupling_transform`](@ref).
"""
function coupled_moments(species::HyperfineOneElectronSpecies, fs_level)
    ops = moment_operators(species, fs_level)
    u = coupling_transform(species, fs_level)
    (x=u' * ops.x * u, y=u' * ops.y * u, z=u' * ops.z * u)
end

# The hyperfine method of the coupled-basis Zeeman Hamiltonian; see the
# docstring in zeeman.jl.
function zeeman_hamiltonian(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    b,
)
    length(b) == 3 || throw(ArgumentError("The field must be a Cartesian 3-vector"))
    for level in basis.levels
        validate_hyperfine(species, level)
    end
    manifolds = unique!([fine_structure(s.level) for s in basis])
    blocks = Dict(
        fs => begin
            ops = coupled_moments(species, fs)
            b[1] .* ops.x .+ b[2] .* ops.y .+ b[3] .* ops.z
        end for fs in manifolds
    )

    num_states = length(basis)
    H = [
        begin
            bra = basis[i]
            ket = basis[k]
            fs = fine_structure(ket.level)
            if fine_structure(bra.level) != fs
                zero(first(first(values(blocks))))
            else
                blocks[fs][
                    manifold_index(species, bra.level, bra.m),
                    manifold_index(species, ket.level, ket.m),
                ]
            end
        end for i in 1:num_states, k in 1:num_states
    ]
    uconvert.(u"rad/µs", H)
end

"""
Returns the index of the state `(level, m)` in the canonical coupled basis of
its fine-structure manifold (``F`` ascending, each ``m_F`` ascending).
"""
function manifold_index(
    species::HyperfineOneElectronSpecies,
    level::HyperfineNumberSpec,
    m,
)
    i = species.nuclear_spin
    f_min = abs(i - level.j)
    offset = sum((Int(2f + 1) for f in f_min:(level.f-1)); init=0)
    offset + Int(m + level.f) + 1
end

"""
    insensitive_field(species, transition::Pair, bracket; abstol = 0u"mT")

Returns a static flux density along ẑ at which the hyperfine transition's
magnetic sensitivity vanishes. `bracket` is a required pair of flux densities
whose sensitivities have opposite signs; there is no default because a
transition may have several roots.

Each trial builds an [`Eigenstructure`](@ref) and evaluates the
Hellmann–Feynman derivative through [`zeeman_sensitivity`](@ref). `abstol` is
the absolute field tolerance of the bracketed root search; zero runs to
floating-point resolution. Physical accuracy also depends on the species data.
"""
function insensitive_field(
    species::HyperfineOneElectronSpecies,
    transition::Pair,
    bracket::NTuple{2,Unitful.BField};
    abstol=0u"mT",
)
    lower, upper = transition
    levels = unique!([
        fine_structure(parse_level(lower.level)),
        fine_structure(parse_level(upper.level)),
    ])
    basis = StateBasis(species, levels...)
    f(b, _) = ustrip(
        u"rad/µs/mT",
        zeeman_sensitivity(Eigenstructure(species, basis; B=b * u"mT"), lower => upper),
    )
    prob = IntervalNonlinearProblem{false}(f, ustrip.(u"mT", bracket))
    sol = solve(prob, ITP(); abstol=ustrip(u"mT", abstol))
    if !successful_retcode(sol)
        throw(
            ArgumentError(
                "No sensitivity zero found over the bracket $bracket " *
                "(solver returned $(sol.retcode))",
            ),
        )
    end
    sol.u * u"mT"
end

export hyperfine_levels, hyperfine_shift, insensitive_field
public coupling_transform, moment_operators, coupled_moments
