# Hyperfine-structure machinery: the coupled |F, m_F⟩ basis is the canonical
# matrix/index convention (hyperfine interaction diagonal, matrices indexable by
# StateBasis{HyperfineNumberSpec}); the |m_I, m_J⟩ product basis appears only as
# an internal construction device, glued by the Clebsch–Gordan unitary of
# coupling_transform.

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

This is the canonical basis the manifold machinery
([`hyperfine_manifold`](@ref)) and the periodic-driving models are expressed
in.
"""
function StateBasis(species::HyperfineOneElectronSpecies, levels...)
    StateBasis([l for entry in levels for l in hyperfine_level_list(species, entry)])
end

"""
    hyperfine_shift(species, level)

Returns the zero-field hyperfine shift of the given ``F`` level relative to the
centroid of its fine-structure level, in angular frequency units.

With ``K = F(F+1) - I(I+1) - J(J+1)``, this is the standard Casimir expression

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
    uconvert(u"µs^-1", e / u"ħ")
end

"""
    moment_operators(species, fs_level) -> (; x, y, z)

Returns the Cartesian components of the magnetic-moment operator
``\\vec{M} = (μ_B / ħ) (g_J \\vec{J} + g_I \\vec{I})`` of the given
fine-structure manifold — the Zeeman Hamiltonian is ``\\vec{B} ⋅ \\vec{M}`` —
in the ``|m_I, m_J⟩`` product basis, in angular frequency units per flux
density.

The product basis is ordered with the nuclear projection slow and the
electronic one fast (`kron(I-space, J-space)`), each in order of increasing
projection; [`Levels.coupling_transform`](@ref) maps it to the canonical
coupled basis.
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
        u"µs^-1/mT",
        (BOHR_MAGNETON / u"ħ") .*
        (g_j .* kron(eye_i, op(j)) .+ g_i .* kron(op(i), eye_j)),
    )
    (x=component(jx_matrix), y=component(jy_matrix), z=component(jz_matrix))
end

"""
    coupling_transform(species, fs_level)

Returns the Clebsch–Gordan unitary ``U`` relating the ``|m_I, m_J⟩`` product
basis of the given fine-structure manifold (cf.
[`Levels.moment_operators`](@ref)) to the canonical coupled ``|F, m_F⟩`` basis
(cf. [`StateBasis(species, levels...)`](@ref StateBasis)):
``U_{(m_I, m_J), (F, m_F)} = ⟨I m_I; J m_J | F m_F⟩``, so a product-basis
operator ``X`` transforms to the coupled basis as ``U^† X U``.

The matrix is real (Condon–Shortley phases) and block-diagonal in
``m_F = m_I + m_J``.
"""
function coupling_transform(species::HyperfineOneElectronSpecies, fs_level)
    spec = fine_structure(parse_level(fs_level))
    i = species.nuclear_spin
    j = spec.j
    d_j = Int(2j + 1)
    n = Int(2i + 1) * d_j

    U = zeros(n, n)
    col = 0
    for f in abs(i-j):(i+j), m_f in (-f):f
        col += 1
        for (i_row, m_i) in enumerate((-i):i), (j_row, m_j) in enumerate((-j):j)
            m_i + m_j == m_f || continue
            U[(i_row-1)*d_j+j_row, col] = Float64(clebschgordan(i, m_i, j, m_j, f, m_f))
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

"""
    zeeman_hamiltonian(species::HyperfineOneElectronSpecies, basis, B)

Returns the matrix of the Zeeman Hamiltonian
``\\vec{B} ⋅ (μ_B / ħ) (g_J \\vec{J} + g_I \\vec{I})`` over the given hyperfine
basis, in angular frequency units.

`B` is the static Cartesian magnetic-field 3-vector, with the quantisation axis
along z; it must be real for the result to be Hermitian. Unlike the
no-hyperfine case, the result is *not* diagonal in the ``F`` levels — the field
mixes ``F`` within each fine-structure manifold (while remaining block-diagonal
across manifolds, and diagonal in ``m_F`` for a field along z). Any basis
ordering or subset of sublevels is supported.
"""
function zeeman_hamiltonian(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    B,
)
    for level in basis.levels
        validate_hyperfine(species, level)
    end
    manifolds = unique!([fine_structure(s.level) for s in basis])
    blocks = Dict(
        fs => begin
            ops = coupled_moments(species, fs)
            B[1] .* ops.x .+ B[2] .* ops.y .+ B[3] .* ops.z
        end for fs in manifolds
    )

    n = length(basis)
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
        end for i in 1:n, k in 1:n
    ]
    uconvert.(u"µs^-1", H)
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
    manifold_hamiltonian(species, fs_level, B)

Returns the hyperfine + Zeeman Hamiltonian of one fine-structure manifold at
the static flux density `B` along the quantisation axis ẑ, over the canonical
coupled basis (`StateBasis(species, fs_level)`), in angular frequency units.

The energy zero is the hyperfine centroid. The hyperfine part is diagonal
(cf. [`hyperfine_shift`](@ref)); the Zeeman part mixes the ``F`` levels while
staying diagonal in ``m_F``.
"""
function manifold_hamiltonian(species::HyperfineOneElectronSpecies, fs_level, B)
    spec = fine_structure(parse_level(fs_level))
    basis = StateBasis(species, spec)
    shifts = [hyperfine_shift(species, s.level) for s in basis]
    moments = coupled_moments(species, spec)
    uconvert.(u"µs^-1", Diagonal(shifts) .+ B .* moments.z)
end

"""
Eigen-solution of one fine-structure manifold of a
[`HyperfineOneElectronSpecies`](@ref) at a static magnetic field along ẑ, with
the eigenstates carrying adiabatic ``(F, m_F)`` labels.

``m_F`` is exact (``[H, F_z] = 0``); ``F`` is a nominal label assigned
adiabatically — by energy order within each ``m_F`` block following the
zero-field ordering of the ``F`` levels, which is rigorous since levels of
equal ``m_F`` do not cross as a function of the field. Both `energies` and the
eigenvector columns of `states` are permuted to align with `basis`, so
[`stateindex`](@ref) addresses them directly and operators rotate into the
eigenbasis as `states' * X * states`.

Constructed via [`hyperfine_manifold`](@ref); requires ``B ≠ 0`` (at zero field
the labels within each degenerate ``F`` level are arbitrary).
"""
struct HyperfineManifold{S<:HyperfineOneElectronSpecies,B<:Quantity,E<:Quantity}
    "The species."
    species::S

    "The fine-structure manifold."
    level::NoHyperfineNumberSpec

    "The static flux density along the quantisation axis ẑ."
    field::B

    "The canonical coupled basis the solution is expressed in."
    basis::StateBasis{HyperfineNumberSpec}

    "Eigen-energies relative to the hyperfine centroid, aligned to `basis`."
    energies::Vector{E}

    "Eigenvector columns in the coupled basis, aligned to `basis`."
    states::Matrix{Float64}
end

"""
    hyperfine_manifold(species, fs_level, B) -> HyperfineManifold

Diagonalises [`manifold_hamiltonian`](@ref) at the static flux density `B ≠ 0`
along ẑ and assigns the adiabatic ``(F, m_F)`` labels (cf.
[`HyperfineManifold`](@ref)).
"""
function hyperfine_manifold(species::HyperfineOneElectronSpecies, fs_level, B)
    if iszero(B)
        throw(
            ArgumentError(
                "Adiabatic (F, m_F) labels are undefined at B = 0; " *
                "a (small) non-zero field is required",
            ),
        )
    end
    spec = fine_structure(parse_level(fs_level))
    basis = StateBasis(species, spec)
    n = length(basis)

    h = ustrip.(u"µs^-1", manifold_hamiltonian(species, spec, B))
    vals, vecs = eigen(Symmetric(real(h)))

    # Fix the arbitrary eigenvector signs for reproducibility: largest-magnitude
    # coupled-basis component positive.
    for k in 1:n
        if vecs[argmax(abs.(view(vecs, :, k))), k] < 0
            vecs[:, k] .*= -1
        end
    end

    # m_F is exact, and F_z is diagonal in the coupled basis.
    fz = [Float64(s.m) for s in basis]
    m_f_labels = [round(2 * sum(abs2.(view(vecs, :, k)) .* fz)) / 2 for k in 1:n]

    # Adiabatic F labels: within each m_F block (energy-ascending from eigen),
    # follow the zero-field energy ordering of the F levels.
    i = species.nuclear_spin
    f_max = i + spec.j
    zero_field_order = sort(
        collect(abs(i-spec.j):f_max);
        by=f -> hyperfine_shift(species, HyperfineNumberSpec(spec.l, spec.j, f)),
    )
    energies = zeros(n)
    states = zeros(n, n)
    for m_f in (-f_max):f_max
        ks = [k for k in 1:n if m_f_labels[k] == m_f]
        available = [f for f in zero_field_order if f >= abs(m_f)]
        @assert length(ks) == length(available)
        for (k, f) in zip(ks, available)
            idx = stateindex(basis, HyperfineNumberSpec(spec.l, spec.j, f), m_f)
            energies[idx] = vals[k]
            states[:, idx] = vecs[:, k]
        end
    end

    HyperfineManifold(species, spec, B, basis, energies .* u"µs^-1", states)
end

"""
    state_energy(m::HyperfineManifold, state)
    state_energy(m::HyperfineManifold, level, m_F)

Returns the eigen-energy (relative to the hyperfine centroid, in angular
frequency units) of the eigenstate with the given adiabatic ``(F, m_F)`` label.
"""
state_energy(m::HyperfineManifold, state::StateSpec) =
    m.energies[stateindex(m.basis, state)]
state_energy(m::HyperfineManifold, level, m_f) =
    m.energies[stateindex(m.basis, level, m_f)]

"""
    eigenbasis_transform(basis::StateBasis, manifolds::HyperfineManifold...)
    eigenbasis_transform(species::HyperfineOneElectronSpecies, basis::StateBasis, B)

Returns the orthogonal matrix ``V`` relating the canonical coupled ``|F, m_F⟩``
basis to the field eigenbasis over the given state basis: column ``k`` holds
the coupled-basis components of the eigenstate carrying the adiabatic
``(F, m_F)`` label `basis[k]` (cf. [`hyperfine_manifold`](@ref)), so a
coupled-basis operator ``X`` over `basis` rotates into the eigenbasis as
``V^† X V``. For example, the exact at-field counterpart of a zero-field
[`quadrupole_couplings`](@ref) matrix `C` is `V' * C * V`, with the basis
states then denoting the adiabatically-labelled eigenstates (individual
components are also available directly via [`transition_amplitude`](@ref)).

The basis may span any number of fine-structure manifolds, in any state order,
but must contain each spanned manifold **completely** — the eigen-solution
lives on the full manifold state space. ``V`` is block-diagonal in the
manifolds. The species form solves each manifold at the static flux density
`B` (along the quantisation axis ẑ); the manifold form takes pre-solved
[`HyperfineManifold`](@ref)s, which must share one species and field. The
column signs follow the eigenvector convention of
[`hyperfine_manifold`](@ref).
"""
function eigenbasis_transform(
    basis::StateBasis{HyperfineNumberSpec},
    manifolds::HyperfineManifold...,
)
    if isempty(manifolds)
        throw(ArgumentError("At least one manifold solution is required"))
    end
    reference = first(manifolds)
    for m in manifolds
        if m.species !== reference.species
            throw(ArgumentError("Manifold solutions must belong to one species"))
        end
        if m.field != reference.field
            throw(
                ArgumentError(
                    "Manifold solutions must share one static field, " *
                    "got $(m.field) and $(reference.field)",
                ),
            )
        end
    end
    by_level = Dict(m.level => m for m in manifolds)
    if length(by_level) != length(manifolds)
        throw(ArgumentError("Duplicate manifold solutions given"))
    end

    # For each basis state: its manifold and canonical index within it.
    slots = map(collect(basis)) do state
        fs = fine_structure(state.level)
        m = get(by_level, fs, nothing)
        if isnothing(m)
            throw(ArgumentError("No manifold solution given for level '$fs'"))
        end
        (; fs, m, index=stateindex(m.basis, state))
    end
    for (fs, m) in by_level
        present = count(slot -> slot.fs == fs, slots)
        iszero(present) && continue
        if present != length(m.basis)
            throw(
                ArgumentError(
                    "Basis must contain the complete '$fs' manifold — the " *
                    "eigenbasis rotation needs the full state space",
                ),
            )
        end
    end

    n = length(basis)
    v = zeros(n, n)
    for k in 1:n, i in 1:n
        slots[i].fs == slots[k].fs || continue
        v[i, k] = slots[k].m.states[slots[i].index, slots[k].index]
    end
    v
end

function eigenbasis_transform(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    B,
)
    manifolds = unique!([fine_structure(s.level) for s in basis])
    eigenbasis_transform(
        basis,
        (hyperfine_manifold(species, fs, B) for fs in manifolds)...,
    )
end

"""
    zeeman_shift(species::HyperfineOneElectronSpecies, state, B)

Returns the Zeeman shift ``E(B) - E(0)`` of the adiabatically-labelled
hyperfine eigenstate, in angular frequency units.

Unlike the no-hyperfine (first-order) case, this is the **exact** shift from
diagonalising the hyperfine + Zeeman Hamiltonian, and is nonlinear in `B` in
the Breit–Rabi regime.
"""
function zeeman_shift(
    species::HyperfineOneElectronSpecies,
    state::StateSpec{HyperfineNumberSpec},
    B,
)
    iszero(B) && return zero(1.0u"µs^-1")
    m = hyperfine_manifold(species, fine_structure(state.level), B)
    uconvert(u"µs^-1", state_energy(m, state) - hyperfine_shift(species, state.level))
end

function zeeman_shift(species::HyperfineOneElectronSpecies, state::StateSpec, B)
    throw(
        ArgumentError(
            "State '$state' must specify a hyperfine (F) level for a hyperfine species",
        ),
    )
end

"""
Returns the field derivative ``dE/dB`` of the given adiabatically-labelled
eigenstate at the manifold's field, via the Hellmann–Feynman theorem:
the expectation value of ``∂H/∂B = (μ_B/ħ)(g_J J_z + g_I I_z)`` in the
eigenstate (exact per ``m_F`` block, as both operators commute with ``F_z``).
"""
function state_moment(m::HyperfineManifold, state)
    v = view(m.states, :, stateindex(m.basis, state))
    moments = coupled_moments(m.species, m.level)
    uconvert(u"µs^-1/mT", dot(v, moments.z, v))
end

"""
    zeeman_sensitivity(species::HyperfineOneElectronSpecies, lower, upper, B)
    zeeman_sensitivity(species::HyperfineOneElectronSpecies, transition::Pair, B)

Returns the magnetic-field sensitivity ``χ = d(E_u - E_l)/dB`` of the
transition between the two adiabatically-labelled hyperfine eigenstates at the
static field `B`, in angular frequency units per flux density.

Evaluated exactly (Hellmann–Feynman, no numerical differentiation) from the
manifold eigenstates, so it is valid in the Breit–Rabi regime; the transition
may lie within one fine-structure manifold or connect two.
"""
function zeeman_sensitivity(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec{HyperfineNumberSpec},
    upper::StateSpec{HyperfineNumberSpec},
    B,
)
    fs_lower = fine_structure(lower.level)
    fs_upper = fine_structure(upper.level)
    m_lower = hyperfine_manifold(species, fs_lower, B)
    m_upper = fs_upper == fs_lower ? m_lower : hyperfine_manifold(species, fs_upper, B)
    uconvert(u"µs^-1/mT", state_moment(m_upper, upper) - state_moment(m_lower, lower))
end

zeeman_sensitivity(species::HyperfineOneElectronSpecies, transition::Pair, B) =
    zeeman_sensitivity(species, transition.first, transition.second, B)

"""
    transition_frequency(species::HyperfineOneElectronSpecies, lower::StateSpec,
                         upper::StateSpec, B)

Returns the exact frequency of the transition between two
adiabatically-labelled hyperfine eigenstates at the static field `B` along ẑ,
from the manifold eigen-energies on top of the centroid splitting (in angular
units).

Unlike the level form, the result is the signed difference upper − lower, as
within one fine-structure manifold either state may lie higher.
"""
function transition_frequency(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec{HyperfineNumberSpec},
    upper::StateSpec{HyperfineNumberSpec},
    B,
)
    fs_lower = fine_structure(lower.level)
    fs_upper = fine_structure(upper.level)
    m_lower = hyperfine_manifold(species, fs_lower, B)
    m_upper = fs_upper == fs_lower ? m_lower : hyperfine_manifold(species, fs_upper, B)
    centroid = (species.energies[fs_upper] - species.energies[fs_lower]) / u"ħ"
    uconvert(
        u"ps^-1",
        centroid + state_energy(m_upper, upper) - state_energy(m_lower, lower),
    )
end

"""
    insensitive_field(species, transition::Pair, bracket)

Returns the static field (along ẑ) at which the first-order magnetic-field
sensitivity of the given `lower => upper` hyperfine transition vanishes, found
by bracketed root finding ([ITP](https://doi.org/10.1145/3423597)) of
[`zeeman_sensitivity`](@ref) over `bracket` (a 2-tuple of flux densities).

There is deliberately no default bracket: sensitivity zeros are typically not
unique (e.g. the ⁴³Ca⁺ 729 nm clock transition has them at both 3.38 G and
4.96 G), so the bracket chooses which one is meant. An error is raised if the
sensitivity does not change sign over the bracket.
"""
function insensitive_field(
    species::HyperfineOneElectronSpecies,
    transition::Pair,
    bracket,
)
    lower, upper = transition
    f(b, _) = ustrip(u"µs^-1/mT", zeeman_sensitivity(species, lower, upper, b * u"mT"))
    prob = IntervalNonlinearProblem{false}(
        f,
        (ustrip(u"mT", bracket[1]), ustrip(u"mT", bracket[2])),
    )
    sol = solve(prob, ITP())
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

export HyperfineManifold,
    hyperfine_levels,
    hyperfine_shift,
    hyperfine_manifold,
    state_energy,
    eigenbasis_transform,
    insensitive_field
public coupling_transform, moment_operators, coupled_moments, manifold_hamiltonian
