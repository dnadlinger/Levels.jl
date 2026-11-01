# The solved level structure: every fine-structure manifold a state basis
# touches, diagonalised once at one static flux density along ẑ. The basis
# states then address the adiabatically labelled eigenstates, and every
# downstream model takes the structure in place of a species, a basis and a
# field.

using LinearAlgebra
using Unitful

"""
One fine-structure manifold solved at a static field along ẑ: its canonical
complete basis (``F`` ascending, ``m`` ascending within each level, or the
``|J, m_J⟩`` states of a fine-structure species), the eigen-energies relative
to the manifold centroid (angular frequencies, aligned to `basis`) and the
eigenvector columns in the coupled basis (aligned to `basis`). For a
fine-structure species the eigenvectors are the identity and the energies the
first-order Zeeman shifts.

Built by [`Levels.solve_manifold`](@ref); the solutions of all manifolds of an
[`Eigenstructure`](@ref) are available through [`Levels.manifold_solution`](@ref).
"""
struct ManifoldSolution{L<:LevelSpec}
    "The fine-structure level."
    level::NoHyperfineNumberSpec

    "The canonical complete basis of the manifold."
    basis::StateBasis{L}

    "Eigen-energies relative to the manifold centroid, aligned to `basis`."
    energies::Vector{typeof(1.0u"rad/µs")}

    "Eigenvector columns in the coupled basis, aligned to `basis`."
    states::Matrix{Float64}
end

identity_states(num_states) = Matrix{Float64}(I, num_states, num_states)

"""
Returns the static flux density as a `Float64` quantity in mT, raising an
error unless it is a real, finite, scalar flux density (zero included).
"""
function validate_field(B)
    if !(B isa Unitful.BField)
        throw(
            ArgumentError(
                "B must be the signed scalar magnetic flux density along the " *
                "quantisation axis ẑ (e.g. 0.5u\"mT\"), not e.g. the Cartesian " *
                "field vector zeeman_hamiltonian takes as a positional argument",
            ),
        )
    end
    if !isreal(B) || !isfinite(B)
        throw(ArgumentError("B must be a real, finite scalar flux density"))
    end
    uconvert(u"mT", 1.0 * real(B))
end

"""
Returns the hyperfine + Zeeman Hamiltonian block of one fine-structure
manifold at the flux density `B` along ẑ over its canonical coupled basis, in
angular frequency units, with the manifold centroid as the energy zero: the
Casimir shifts on the diagonal (zero for a fine-structure species) plus `B`
times the ẑ component of the magnetic moment.
"""
function manifold_hamiltonian_block(species::NoHyperfineOneElectronSpecies, fs, B)
    basis = StateBasis([fs])
    g = lande_g(species, fs)
    Matrix(
        Diagonal([
            uconvert(u"rad/µs", g * s.m * BOHR_MAGNETON * B / u"ħ") for s in basis
        ]),
    )
end

function manifold_hamiltonian_block(species::HyperfineOneElectronSpecies, fs, B)
    basis = StateBasis(species, fs)
    casimir = [hyperfine_shift(species, s.level) for s in basis]
    uconvert.(u"rad/µs", Diagonal(casimir) .+ B .* coupled_moments(species, fs).z)
end

"""
    solve_manifold(species, fs_level, B; require_constants = true) -> ManifoldSolution

Solves one complete fine-structure manifold at the static flux density `B`
along ẑ, returning centroid-relative angular energies and coupled-basis
eigenvectors.

For a fine-structure species the eigenvectors are the identity and the
energies are linear Zeeman shifts. For a hyperfine species the hyperfine +
Zeeman Hamiltonian is diagonalised separately in each conserved-``m_F`` block.
Within a block, ``F`` labels follow the energy order of the zero-field levels.
At `B = 0` the coupled states and their Casimir energies are returned directly.
Ties between zero-field levels are ordered by ascending ``F``; these labels
have no unique adiabatic continuation.

Each eigenvector's component on its own zero-field label is made positive.
If it is exactly zero, the largest-magnitude component is made positive.
This gives the identity limit for a resolved zero-field spectrum, but permits
sign jumps at finite field; see [Bases and representations](@ref bases).

Missing hyperfine constants raise an error. With `require_constants = false`,
that manifold is instead returned degenerate at its centroid, with identity
eigenvectors and no Zeeman shift. This approximation is used only for
far-detuned polarisability intermediates.
"""
function solve_manifold(
    species::NoHyperfineOneElectronSpecies,
    fs::NoHyperfineNumberSpec,
    B;
    require_constants=true,
)
    h = manifold_hamiltonian_block(species, fs, B)
    basis = StateBasis([fs])
    ManifoldSolution(fs, basis, collect(diag(h)), identity_states(length(basis)))
end

function solve_manifold(
    species::HyperfineOneElectronSpecies,
    fs::NoHyperfineNumberSpec,
    B;
    require_constants=true,
)
    basis = StateBasis(species, fs)
    num_states = length(basis)
    if !haskey(species.hyperfine, fs)
        if require_constants
            throw(
                ArgumentError(
                    "No hyperfine constants known for level '$fs', so its F " *
                    "levels cannot be located",
                ),
            )
        end
        return ManifoldSolution(
            fs,
            basis,
            zeros(num_states) .* u"rad/µs",
            identity_states(num_states),
        )
    end
    if iszero(B)
        casimir =
            [uconvert(u"rad/µs", hyperfine_shift(species, s.level)) for s in basis]
        return ManifoldSolution(fs, basis, casimir, identity_states(num_states))
    end

    h = ustrip.(u"rad/µs", manifold_hamiltonian_block(species, fs, B))
    i_nuc = species.nuclear_spin
    j = fs.j
    f_max = i_nuc + j
    # Zero-field energy order of the F levels (stable: ties stay in ascending F).
    zero_field_order = sort(
        collect(abs(i_nuc-j):f_max);
        by=f -> hyperfine_shift(species, HyperfineNumberSpec(fs.l, j, f)),
    )
    energies = zeros(num_states)
    states = zeros(num_states, num_states)
    for m_f in (-f_max):f_max
        available = [f for f in zero_field_order if f >= abs(m_f)]
        block =
            [stateindex(basis, HyperfineNumberSpec(fs.l, j, f), m_f) for f in available]
        vals, vecs = eigen(Symmetric(h[block, block]))
        for (k, col) in enumerate(block)
            v = vecs[:, k]
            own = v[k]
            reference = iszero(own) ? v[argmax(abs.(v))] : own
            reference < 0 && (v = -v)
            energies[col] = vals[k]
            states[block, col] = v
        end
    end
    ManifoldSolution(fs, basis, energies .* u"rad/µs", states)
end

"""
The solved hyperfine + Zeeman structure of a species at a static magnetic
flux density along ẑ.

Constructed with [`Eigenstructure(species, basis; B)`](@ref Eigenstructure(::OneElectronSpecies, ::StateBasis)).
Each fine-structure manifold touched by the basis is solved in full; the basis
selects and orders the labelled eigenstates used by operator matrices. For a
hyperfine species these are field eigenstates labelled by their zero-field
``F`` branch and conserved ``m_F``; for a fine-structure species they are
``|J, m_J⟩`` states with linear Zeeman energies.

[`offsets`](@ref) gives their centroid-relative angular energies,
[`eigenvectors`](@ref) their coupled-basis components, and
[`to_eigenbasis`](@ref) transforms coupled-basis operators. See
[Bases and representations](@ref bases) for labels, phases and partial bases.
Treat the structure and its stored arrays as read-only; accessors return them
directly.
"""
struct Eigenstructure{S<:OneElectronSpecies,L<:LevelSpec}
    "The species."
    species::S

    "The selected basis: any subset of the manifolds' states, in any order."
    basis::StateBasis{L}

    "The static flux density along the quantisation axis ẑ (signed; zero allowed)."
    B::typeof(1.0u"mT")

    "The solution of every fine-structure manifold the basis touches."
    manifolds::Dict{NoHyperfineNumberSpec,ManifoldSolution{L}}

    "The complete canonical zero-field basis of those manifolds, in order of first appearance."
    coupled_basis::StateBasis{L}

    "Eigen-energy of each basis state relative to its fine-structure centroid."
    offsets::Vector{typeof(1.0u"rad/µs")}

    "Eigenvector matrix `W`: rows over `coupled_basis`, columns over `basis`."
    eigenvectors::Matrix{Float64}
end

"""
    Eigenstructure(species, basis::StateBasis; B)
    Eigenstructure(species, levels...; B)
    Eigenstructure(species, states::AbstractVector{<:StateSpec}; B)

Solves the level structure of `species` at the static flux density `B` along
the quantisation axis ẑ (a signed scalar, e.g. `0.5u"mT"`; `0u"mT"` for the
zero-field structure — there is deliberately no default) for the states of
`basis`. The level form expands each level into all its states via
[`StateBasis(species, levels...)`](@ref StateBasis) (a fine-structure level
of a hyperfine species into all its ``F`` levels); the state form fixes an
explicit subset and order.

The basis kind must match the species (``F`` levels for a hyperfine species,
``J`` levels for a fine-structure one), every basis level must be a level of
the species, and — for a hyperfine species — every basis manifold must carry
hyperfine constants, since without them its ``F`` levels cannot be located.
Any subset of states is admitted: the eigenproblem is solved on the complete
manifolds regardless, so partial bases address exact eigenstates (cf.
[`iscomplete`](@ref) for the consumers that need whole manifolds).
"""
function Eigenstructure(species::OneElectronSpecies, basis::StateBasis; B)
    validate_basis_kind(species, basis)
    field = validate_field(B)
    for level in basis.levels
        fs = fine_structure(level)
        haskey(species.energies, fs) ||
            throw(ArgumentError("Level '$fs' is not a level of the species"))
        level isa HyperfineNumberSpec && validate_hyperfine(species, level)
    end
    fs_levels = unique!([fine_structure(level) for level in basis.levels])
    manifolds = Dict(fs => solve_manifold(species, fs, field) for fs in fs_levels)
    coupled = StateBasis(reduce(vcat, (manifolds[fs].basis.states for fs in fs_levels)))

    num_states = length(basis)
    offsets = Vector{typeof(1.0u"rad/µs")}(undef, num_states)
    W = zeros(length(coupled), num_states)
    for (k, state) in enumerate(basis)
        fs = fine_structure(state.level)
        m = manifolds[fs]
        index = stateindex(m.basis, state)
        offsets[k] = m.energies[index]
        W[staterange(coupled, fs), k] = m.states[:, index]
    end
    Eigenstructure(species, basis, field, manifolds, coupled, offsets, W)
end

Eigenstructure(species::OneElectronSpecies, states::AbstractVector{<:StateSpec}; B) =
    Eigenstructure(species, StateBasis(states); B)
Eigenstructure(species::OneElectronSpecies, levels::AbstractVector; B) =
    Eigenstructure(species, StateBasis(species, levels...); B)
Eigenstructure(species::OneElectronSpecies, levels...; B) =
    Eigenstructure(species, StateBasis(species, levels...); B)

"""
    species(s::Eigenstructure)

Returns the species of the structure.
"""
species(s::Eigenstructure) = s.species

"""
    Levels.basis(s::Eigenstructure)

Returns the selected [`StateBasis`](@ref) of the structure: the eigenstates
every operator matrix built on `s` is indexed by, in that order. Not
exported, as QuantumToolbox exports a `basis` of its own; `s.basis` reads the
same field.
"""
basis(s::Eigenstructure) = s.basis

"""
    field(s::Eigenstructure)

Returns the static flux density along ẑ the structure was solved at, in mT.
"""
field(s::Eigenstructure) = s.B

"""
    offsets(s::Eigenstructure)

Returns the eigen-energy of every basis state relative to the zero-field
centroid of its fine-structure level, as angular frequencies aligned to
`basis(s)` (cf. [`state_offset`](@ref)).
"""
offsets(s::Eigenstructure) = s.offsets

"""
    coupled_basis(s::Eigenstructure)

Returns the complete canonical zero-field basis of the structure's manifolds —
the index space of [`coupled_hamiltonian`](@ref), of the rows of
[`eigenvectors`](@ref) and of the input of [`to_eigenbasis`](@ref): every
manifold the basis touches, in order of first appearance, each with its
``F`` levels ascending and ``m`` ascending within each (the ``|J, m_J⟩``
states of a fine-structure species).
"""
coupled_basis(s::Eigenstructure) = s.coupled_basis

"""
    eigenvectors(s::Eigenstructure)

Returns the real eigenvector matrix ``W`` of the structure: column `k` holds
the coupled-basis components (rows indexed by [`coupled_basis`](@ref)) of the
eigenstate labelled `basis(s)[k]`. The columns are orthonormal and supported
within their own manifold; for a complete basis in canonical order ``W`` is
square and orthogonal, for a reordered or partial basis a column permutation
or subset of it. A coupled-basis operator rotates into the eigenbasis as
``W^† X W`` ([`to_eigenbasis`](@ref)).

For a fine-structure species ``W`` is a column selection of the identity. At
zero field it is the same for any species. The column signs follow the
own-label-positive convention of [`Levels.solve_manifold`](@ref): a whole-column
sign change is a gauge choice that leaves every consistently computed
observable unchanged.
"""
eigenvectors(s::Eigenstructure) = s.eigenvectors

"""
    to_eigenbasis(s::Eigenstructure, X)

Rotates the matrix `X`, given over the complete coupled basis
[`coupled_basis`](@ref)`(s)`, into the field eigenbasis of the structure:
``W^† X W`` with ``W`` = [`eigenvectors`](@ref)`(s)`, indexed by `basis(s)`.
For a partial basis the result holds the exact matrix elements between the
retained eigenstates; a matrix over `basis(s)` itself is refused, as the
eigenstates have components on coupled states outside the subset, so no exact
rotation exists from that input. Projection does not preserve products: for a
partial basis ``(W^† X W)(W^† Y W) ≠ W^† X Y W`` in general.
"""
function to_eigenbasis(s::Eigenstructure, X::AbstractMatrix)
    num_coupled = length(s.coupled_basis)
    if size(X) != (num_coupled, num_coupled)
        throw(
            ArgumentError(
                "The matrix must be given over the complete coupled basis of the " *
                "structure ($num_coupled states; see coupled_basis), not over its " *
                "selected basis: the eigenstates have components on coupled states " *
                "outside a partial basis, so there is no exact rotation from a " *
                "$(size(X, 1)) × $(size(X, 2)) input",
            ),
        )
    end
    W = s.eigenvectors
    W' * X * W
end

# The fine-structure manifolds of the structure, in order of first appearance.
manifold_levels(s::Eigenstructure) =
    unique!([fine_structure(level) for level in s.coupled_basis.levels])

"""
    coupled_hamiltonian(s::Eigenstructure)

Returns the static offset Hamiltonian ``H/ħ`` over [`coupled_basis`](@ref)`(s)`,
in angular-frequency units. Each fine-structure block contains the Casimir
hyperfine shifts plus `field(s)` times the magnetic operator
[`Levels.coupled_moments`](@ref)`(species(s), level).z`.

Each manifold's centroid is subtracted, giving
`to_eigenbasis(s, coupled_hamiltonian(s)) ≈ Diagonal(offsets(s))`.
[`state_energy`](@ref) adds the centroid; [`transition_frequency`](@ref)
forms intervals between labelled eigenstates.
"""
function coupled_hamiltonian(s::Eigenstructure)
    cb = s.coupled_basis
    num_states = length(cb)
    H = zeros(num_states, num_states) .* u"rad/µs"
    for fs in manifold_levels(s)
        r = staterange(cb, fs)
        H[r, r] = manifold_hamiltonian_block(s.species, fs, s.B)
    end
    H
end

"""
    iscomplete(s::Eigenstructure, level)

Returns whether `basis(s)` contains every state of the fine-structure manifold
enclosing `level`, which may name a single ``F`` level. Returns `false` for an
absent manifold. Consumers requiring complete manifolds include
`Levels.PeriodicDynamics.DrivenTransition` and
`Levels.PeriodicDynamics.quasienergy_shift`.
"""
function iscomplete(s::Eigenstructure, level)
    fs = fine_structure(parse_level(level))
    m = get(s.manifolds, fs, nothing)
    isnothing(m) && return false
    count(state -> fine_structure(state.level) == fs, s.basis.states) == length(m.basis)
end

"""
    Levels.manifold_solution(s::Eigenstructure, level)

Returns the [`Levels.ManifoldSolution`](@ref) of the fine-structure manifold
enclosing `level`: the stored one for a manifold of the structure, a fresh
solve at `field(s)` otherwise (which, for a hyperfine species, requires the
level's hyperfine constants).
"""
function manifold_solution(s::Eigenstructure, level)
    fs = fine_structure(parse_level(level))
    get(s.manifolds, fs) do
        solve_manifold(s.species, fs, s.B)
    end
end

kind_name(::Type{NoHyperfineNumberSpec}) = "fine-structure (J)"
kind_name(::Type{HyperfineNumberSpec}) = "hyperfine (F)"

"""
Resolves a state label to the manifold solution it belongs to and its index in
that manifold's canonical basis, accepting any state of `coupled_basis(s)`
(not only basis states). Raises for a state of the wrong kind, of a manifold
the structure does not contain, or with an invalid projection.
"""
function resolve_state(s::Eigenstructure{S,L}, state::StateSpec) where {S,L}
    spec = parse_level(state.level)
    if !(spec isa L)
        throw(
            ArgumentError(
                "State '$state' must specify a $(kind_name(L)) level like the " *
                "states of the structure",
            ),
        )
    end
    fs = fine_structure(spec)
    m = get(s.manifolds, fs, nothing)
    if isnothing(m)
        throw(
            ArgumentError(
                "Level '$fs' is not a manifold of the structure " *
                "($(join(string.(manifold_levels(s)), ", ")))",
            ),
        )
    end
    m, stateindex(m.basis, StateSpec(spec, state.m))
end

"""
    state_offset(s::Eigenstructure, state)
    state_offset(s::Eigenstructure, level, m)

Returns the eigen-energy of the labelled eigenstate relative to the zero-field
centroid of its fine-structure level, in angular frequency units — for a
hyperfine species the exact hyperfine + Zeeman eigen-energy of the
adiabatically labelled state, for a fine-structure species the first-order
Zeeman shift. Any state of [`coupled_basis`](@ref)`(s)` may be queried, not
only basis states. [`state_energy`](@ref) adds the centroid, giving an
absolute energy in joules.
"""
function state_offset(s::Eigenstructure, state::StateSpec)
    m, index = resolve_state(s, state)
    m.energies[index]
end
state_offset(s::Eigenstructure, level, m) = state_offset(s, StateSpec(level, m))

"""
    state_energy(s::Eigenstructure, state)
    state_energy(s::Eigenstructure, level, m)

Returns the energy of the labelled eigenstate in joules, relative to the
species' common energy reference: the centroid energy of its fine-structure
level ([`level_energy`](@ref)) plus ``ħ`` times [`state_offset`](@ref).
Subtracting total energies to form a small splitting can lose precision;
[`transition_frequency`](@ref) keeps centroid and offset terms separate.
"""
function state_energy(s::Eigenstructure, state::StateSpec)
    m, index = resolve_state(s, state)
    uconvert(u"J", level_energy(s.species, m.level) + u"ħ" * m.energies[index])
end
state_energy(s::Eigenstructure, level, m) = state_energy(s, StateSpec(level, m))

"""
    zeeman_shift(s::Eigenstructure, state)

Returns the Zeeman shift ``[E(B) - E(0)]/ħ`` of the labelled eigenstate at the
structure's field, in angular frequency units: [`state_offset`](@ref) minus the
zero-field hyperfine shift of the state's level. Exact for both species
kinds (for a hyperfine species nonlinear in the field in the Breit–Rabi
regime); zero at `B = 0`.
"""
function zeeman_shift(s::Eigenstructure, state::StateSpec)
    m, index = resolve_state(s, state)
    label = m.basis[index]
    uconvert(u"rad/µs", m.energies[index] - hyperfine_shift(s.species, label.level))
end

# The field derivative dE/dB of a labelled eigenstate at the structure's
# field: the Hellmann–Feynman expectation value of ∂H/∂B, exact per m block.
function state_moment(s::Eigenstructure{<:NoHyperfineOneElectronSpecies}, state)
    m, index = resolve_state(s, state)
    label = m.basis[index]
    uconvert(
        u"rad/µs/mT",
        lande_g(s.species, label.level) * label.m * BOHR_MAGNETON / u"ħ",
    )
end

function state_moment(s::Eigenstructure{<:HyperfineOneElectronSpecies}, state)
    m, index = resolve_state(s, state)
    v = view(m.states, :, index)
    uconvert(u"rad/µs/mT", dot(v, coupled_moments(s.species, m.level).z, v))
end

"""
    zeeman_sensitivity(s::Eigenstructure, lower::StateSpec => upper::StateSpec)

Returns the magnetic-field sensitivity ``χ = (1/ħ) d(E_u - E_l)/dB`` of the
transition between the two labelled eigenstates at the structure's field, in
angular frequency units per flux density. Exact (Hellmann–Feynman theorem on
the eigenstates, no numerical differentiation), so valid in the Breit–Rabi
regime of a hyperfine species; the transition may lie within one
fine-structure manifold or connect two. For a fine-structure species it is the
field-independent ``(g_u m_u - g_l m_l) μ_B/ħ``.
"""
function zeeman_sensitivity(
    s::Eigenstructure,
    transition::Pair{<:StateSpec,<:StateSpec},
)
    lower, upper = transition
    uconvert(u"rad/µs/mT", state_moment(s, upper) - state_moment(s, lower))
end

"""
    transition_frequency(s::Eigenstructure, lower::StateSpec => upper::StateSpec)

Returns the angular frequency of the transition between the two labelled
eigenstates at the structure's field: the centroid interval of their
fine-structure levels plus the difference of their [`state_offset`](@ref)s,
formed with the two terms separated so that same-manifold (microwave)
intervals and small differences keep their precision. A pair of levels
forwards to the zero-field species form.
"""
function transition_frequency(s::Eigenstructure, transition::Pair)
    lower, upper = transition
    if lower isa StateSpec && upper isa StateSpec
        m_lo, i = resolve_state(s, lower)
        m_hi, k = resolve_state(s, upper)
        centroid =
            (
                level_energy(s.species, m_hi.level) -
                level_energy(s.species, m_lo.level)
            ) / u"ħ"
        return uconvert(u"rad/µs", centroid + m_hi.energies[k] - m_lo.energies[i])
    end
    transition_frequency(s.species, transition)
end

"""
    zeeman_hamiltonian(s::Eigenstructure, b)

Returns the Zeeman Hamiltonian of a magnetic field with the Cartesian 3-vector
`b` (a perturbing field, e.g. the amplitude of an ac field — not the static
field, which is fixed in `s`) in the field eigenbasis of the structure,
indexed by `basis(s)`: the coupled-basis matrix of
[`zeeman_hamiltonian(species, basis, b)`](@ref) over
[`coupled_basis`](@ref)`(s)` rotated with [`to_eigenbasis`](@ref). Exact
matrix elements between the retained eigenstates for any subset or order.
"""
zeeman_hamiltonian(s::Eigenstructure, b) = uconvert.(
    u"rad/µs",
    to_eigenbasis(s, zeeman_hamiltonian(s.species, s.coupled_basis, b)),
)

Base.show(io::IO, s::Eigenstructure) = print(
    io,
    "Eigenstructure(",
    join(string.(manifold_levels(s)), ", "),
    "; ",
    length(s.basis),
    " states; B = ",
    s.B,
    ")",
)

export Eigenstructure,
    species,
    field,
    offsets,
    coupled_basis,
    eigenvectors,
    to_eigenbasis,
    coupled_hamiltonian,
    iscomplete,
    state_offset,
    state_energy
public basis, ManifoldSolution, manifold_solution, solve_manifold
