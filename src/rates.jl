using LinearAlgebra: I
using StaticArrays: SVector
using Unitful
using WignerSymbols

"""
    parity_multipole(lower => upper)

Returns the lowest electric multipole allowed by the parities of the two given
fine-structure levels, as determined from their orbital angular momenta: `:E1` if
the parities differ, `:E2` if they are equal.

This is a selection rule, not the multipolarity of a transition: a same-parity
pair may decay by M1 instead, which only the species data records (cf.
[`Levels.multipole`](@ref)). It serves where electric-dipole character is the
question, and as the default for rate entries not marked otherwise.
"""
function parity_multipole(transition::Pair)
    lower, upper = transition
    lo = convert(NoHyperfineNumberSpec, lower)
    hi = convert(NoHyperfineNumberSpec, upper)
    isodd(lo.l + hi.l) ? :E1 : :E2
end

"""
    multipole(species, lower => upper)

Returns the multipolarity of the transition between the two given levels — `:E1`,
`:E2` or `:M1` — as recorded with the pair's Einstein A coefficient (in either
orientation; a [`MagneticDipole`](@ref) entry is M1, any other the lowest electric
multipole allowed by parity). A pair without a tabulated rate gets that parity
default, [`Levels.parity_multipole`](@ref). Hyperfine levels resolve to their
fine-structure part.
"""
function multipole(species::OneElectronSpecies, transition::Pair)
    lo = fine_structure(transition.first)
    hi = fine_structure(transition.second)
    get(species.multipoles, lo => hi) do
        get(() -> parity_multipole(lo => hi), species.multipoles, hi => lo)
    end
end

"""
    multipole_rank(multipole::Symbol)
    multipole_rank(species, lower => upper)

Returns the tensor rank ``k`` of the given multipole — 1 for `:E1` and `:M1`, 2 for
`:E2` — or of the transition between two levels of a species, as per
[`Levels.multipole`](@ref).
"""
function multipole_rank(multipole::Symbol)
    multipole in (:E1, :M1) && return 1
    multipole == :E2 && return 2
    throw(ArgumentError("Unknown multipole :$multipole; expected :E1, :E2 or :M1"))
end

multipole_rank(species::OneElectronSpecies, transition::Pair) =
    multipole_rank(multipole(species, transition))

"""
    clebsch_gordan(lower::StateSpec => upper::StateSpec)

Returns the bare Clebsch–Gordan coefficient ``⟨j m; R Δm | j' m'⟩`` between two
fine-structure states, with the rank ``R`` that of the lowest electric multipole
allowed by parity ([`Levels.parity_multipole`](@ref)).

This is the purely geometric coupling factor; the general transition amplitude
relative to the reduced matrix element — species-aware, covering hyperfine
``F`` levels and their exact mixing at a static field — is
[`transition_amplitude`](@ref).
"""
function clebsch_gordan(transition::Pair{<:StateSpec,<:StateSpec})
    lower, upper = transition
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    rank = multipole_rank(parity_multipole(lo => hi))
    Δm = upper.m - lower.m
    abs(Δm) <= rank || return 0.0
    Float64(clebschgordan(lo.j, lower.m, rank, Δm, hi.j, upper.m))
end

# Triangle condition (|a - b| ≤ c ≤ a + b with integer perimeter) for angular
# momentum coupling.
is_triangle(a, b, c) = abs(a - b) <= c <= a + b && isinteger(a + b - c)

"""
    transition_amplitude(species, lower::StateSpec => upper::StateSpec; rank)
    transition_amplitude(s::Eigenstructure, lower::StateSpec => upper::StateSpec; rank)

Returns the dimensionless rank-`k` transition amplitude
``⟨u|T^k_q|l⟩ / (⟨J'‖T^k‖J⟩/\\sqrt{2J'+1})``, with ``q = m_u - m_l``.
`rank` defaults to that of the pair's tabulated multipole
([`Levels.multipole_rank`](@ref)); an explicit rank supports other operators,
e.g. `1` for M1 within a fine-structure level.

The species form uses zero-field states: the Clebsch–Gordan coefficient
``⟨J m; k q | J' m'⟩`` for fine structure, or
``⟨F m; k q | F' m'⟩ β^{(k)}(F → F')`` for hyperfine structure, with ``β``
from [`Levels.hyperfine_reduction`](@ref).

The structure form rotates the complete-manifold tensor component into the
field eigenbases before selecting the labelled states, including their ``F``
mixing. Any state of `coupled_basis(s)` may be named. Signs follow the
structure's eigenvector gauge; [Conventions](@ref conventions)
describes its zero-field limit and possible sign jumps in field sweeps.
"""
function transition_amplitude(
    species::OneElectronSpecies,
    transition::Pair{<:StateSpec,<:StateSpec};
    rank=nothing,
)
    lower = validate_state_kind(species, transition.first)
    upper = validate_state_kind(species, transition.second)
    zero_field_amplitude(species, lower, upper, rank)
end

function zero_field_amplitude(
    species::NoHyperfineOneElectronSpecies,
    lower,
    upper,
    rank,
)
    lo = lower.level
    hi = upper.level
    k = something(rank, multipole_rank(species, lo => hi))
    Δm = upper.m - lower.m
    (abs(Δm) <= k && is_triangle(lo.j, k, hi.j)) || return 0.0
    Float64(clebschgordan(lo.j, lower.m, k, Δm, hi.j, upper.m))
end

function zero_field_amplitude(species::HyperfineOneElectronSpecies, lower, upper, rank)
    validate_hyperfine(species, lower.level)
    validate_hyperfine(species, upper.level)
    k = something(rank, multipole_rank(species, lower.level => upper.level))
    Δm = upper.m - lower.m
    abs(Δm) <= k || return 0.0
    β = hyperfine_reduction(species.nuclear_spin, lower.level => upper.level; rank=k)
    iszero(β) && return 0.0
    Float64(clebschgordan(lower.level.f, lower.m, k, Δm, upper.level.f, upper.m)) * β
end

function transition_amplitude(
    s::Eigenstructure,
    transition::Pair{<:StateSpec,<:StateSpec};
    rank=nothing,
)
    m_lo, i = resolve_state(s, transition.first)
    m_hi, k = resolve_state(s, transition.second)
    relative_amplitudes(
        s.species,
        m_lo,
        m_hi,
        something(rank, multipole_rank(s.species, m_lo.level => m_hi.level)),
    )[
        k,
        i,
    ]
end

"""
Returns the multiplier of the rank-`rank` Clebsch–Gordan coefficient between
two levels: the triangle condition (0 or 1) for fine-structure levels, or the
[`Levels.hyperfine_reduction`](@ref) factor ``β`` for hyperfine levels.
"""
reduction_factor(species, lo::NoHyperfineNumberSpec, hi::NoHyperfineNumberSpec, rank) =
    is_triangle(lo.j, rank, hi.j) ? 1.0 : 0.0

reduction_factor(
    species::HyperfineOneElectronSpecies,
    lo::HyperfineNumberSpec,
    hi::HyperfineNumberSpec,
    rank,
) = hyperfine_reduction(species.nuclear_spin, lo => hi; rank)

"""
    relative_amplitudes(species, source::ManifoldSolution, target::ManifoldSolution, rank)

Returns the `[target, source]` matrix of rank-`rank` transition amplitudes
relative to the fine-structure reduced matrix element between the eigenstates
of two manifold solutions (over their canonical bases): the complete-manifold
zero-field matrix of [`transition_amplitude`](@ref) values, rotated into both
field eigenbases. This is the one primitive every at-field coupling matrix,
Rabi frequency, physical dipole component and decay operator selects from —
a restriction to particular ``F`` levels or states is always applied *after*
this rotation. Definite ``m`` makes each entry belong to a single component
``q = m_{target} - m_{source}``.
"""
function relative_amplitudes(
    species,
    source::ManifoldSolution,
    target::ManifoldSolution,
    rank,
)
    c = zeros(length(target.basis), length(source.basis))
    for lo in source.basis.levels, hi in target.basis.levels
        β = reduction_factor(species, lo, hi, rank)
        iszero(β) && continue
        j_lo, j_hi = momentum(lo), momentum(hi)
        for i in staterange(source.basis, lo), k in staterange(target.basis, hi)
            q = target.basis[k].m - source.basis[i].m
            abs(q) <= rank || continue
            cg =
                clebschgordan(j_lo, source.basis[i].m, rank, q, j_hi, target.basis[k].m)
            c[k, i] = Float64(cg) * β
        end
    end
    target.states' * c * source.states
end

# Returns a function (lower_state, upper_state) -> relative amplitude for
# states of the two given fine-structure manifolds of the structure, from one
# rotated matrix.
function amplitude_lookup(s::Eigenstructure, fs_lo, fs_hi, rank)
    m_lo = manifold_solution(s, fs_lo)
    m_hi = manifold_solution(s, fs_hi)
    a = relative_amplitudes(s.species, m_lo, m_hi, rank)
    (lo, hi) -> a[stateindex(m_hi.basis, hi), stateindex(m_lo.basis, lo)]
end

"""
    hyperfine_reduction(nuclear_spin, lower => upper; rank)
    hyperfine_reduction(species::HyperfineOneElectronSpecies, lower => upper; rank)

Returns the factor relating an electronic rank-`k` multipole amplitude between
hyperfine states to its fine-structure reduced element:

```math
β^{(k)}(F → F') = (-1)^{J' + I + F + k}
    \\sqrt{(2F + 1)(2J' + 1)}
    \\begin{Bmatrix} J & J' & k \\\\ F' & F & I \\end{Bmatrix},
```

so ``⟨F' m'|T^k_q|F m⟩ = ⟨F m; k q|F' m'⟩ β^{(k)}
⟨J'‖T^k‖J⟩/\\sqrt{2J'+1}``. The coupled states use electron-first
Condon–Shortley coefficients ``⟨J m_J; I m_I | F m_F⟩``, as in
[`Levels.coupling_transform`](@ref) and [Campbell2026; Appendix E](@cite).
Nuclear-first states differ by ``(-1)^{I+J-F}`` per state; amplitudes between
them acquire the corresponding initial- and final-state phases.

For fixed upper ``F'``, ``\\sum_F β^2 = 1``: ``β^2`` gives the branching
fractions, and every upper hyperfine sublevel has the fine-structure decay
rate. Species tables therefore store fine-structure Einstein A coefficients.
For ``I = 0``, ``β = +1`` on allowed transitions.

`rank` defaults to that of the lowest electric multipole allowed by parity
([`Levels.parity_multipole`](@ref)), or for the species form to that of the pair's
tabulated multipole ([`Levels.multipole_rank`](@ref)); an explicit rank supports
other operators such as M1 within a fine-structure level.
"""
function hyperfine_reduction(nuclear_spin, transition::Pair; rank=nothing)
    lower, upper = transition
    lo = parse_level(lower)
    hi = parse_level(upper)
    if !(lo isa HyperfineNumberSpec && hi isa HyperfineNumberSpec)
        throw(ArgumentError("Levels must specify hyperfine (F) levels"))
    end
    k = something(
        rank,
        multipole_rank(parity_multipole(fine_structure(lo) => fine_structure(hi))),
    )
    if !is_triangle(lo.f, k, hi.f) || !is_triangle(lo.j, k, hi.j)
        return 0.0
    end
    Float64(
        (-1)^Int(hi.j + nuclear_spin + lo.f + k) *
        sqrt((2 * lo.f + 1) * (2 * hi.j + 1)) *
        wigner6j(lo.j, hi.j, k, hi.f, lo.f, nuclear_spin),
    )
end

function hyperfine_reduction(
    species::HyperfineOneElectronSpecies,
    transition::Pair;
    rank=nothing,
)
    lower, upper = transition
    lo = parse_level(lower)
    hi = parse_level(upper)
    lo isa HyperfineNumberSpec && validate_hyperfine(species, lo)
    hi isa HyperfineNumberSpec && validate_hyperfine(species, hi)
    if isnothing(rank) && lo isa HyperfineNumberSpec && hi isa HyperfineNumberSpec
        rank = multipole_rank(species, lo => hi)
    end
    hyperfine_reduction(species.nuclear_spin, lo => hi; rank)
end

# The levels a level argument stands for in a basis of the species' kind: the
# single fine-structure level, or — for a hyperfine species — a single F level
# or all F levels of a fine-structure manifold.
level_list(species::NoHyperfineOneElectronSpecies, level) =
    [validate_level_kind(species, level)]
level_list(species::HyperfineOneElectronSpecies, level) =
    hyperfine_level_list(species, level)

# Resolves a `lower => upper` level pair to the lists of levels on each side
# and their fine-structure manifolds, multipole and multipole rank.
function level_pair(species, pair::Pair)
    lower_levels = level_list(species, pair.first)
    upper_levels = level_list(species, pair.second)
    fs_lo = fine_structure(first(lower_levels))
    fs_hi = fine_structure(first(upper_levels))
    m = multipole(species, fs_lo => fs_hi)
    (; lower_levels, upper_levels, fs_lo, fs_hi, multipole=m, rank=multipole_rank(m))
end

# The absolute scale of the James formula for one fine-structure pair at the
# given intensity: Ω = scale × amplitude × geometry, with the Einstein A
# coefficient and the **centroid** transition frequency of the pair, so that
# one physical normalisation serves every component of the pair. M1 shares the
# E1 prefactor: the 1/c of the wave's magnetic field cancels against the c⁻⁵ of
# the M1 rate.
function rabi_scale(species, fs_lo, fs_hi, intensity)
    a = einstein_a(species, fs_lo => fs_hi)
    isnothing(a) &&
        throw(ArgumentError("No known transition between '$fs_lo' and '$fs_hi'"))
    ω = transition_frequency(species, fs_lo => fs_hi)
    prefactor = multipole(species, fs_lo => fs_hi) == :E2 ? 20.0 : 6.0
    uconvert(u"rad/µs", sqrt(prefactor * π * u"c"^2 * intensity * a / (u"ħ" * ω^3)))
end

# The normalised geometric channel amplitudes of an E1, E2 or M1 coupling,
# indexed by q + rank + 1. An M1 coupling sees the magnetic polarisation
# n̂ × ε̂ of the wave, B = (n̂ × E)/c.
function channel_geometry(multipole, ε, direction)
    ε_scale = sqrt(sum(abs2, ε))
    iszero(ε_scale) && throw(ArgumentError("Polarisation must be non-zero"))
    ε̂ = ε ./ ε_scale
    multipole == :E1 && return dipole_geometry(ε̂)
    isnothing(direction) && throw(
        ArgumentError(
            "The beam direction is required for an $multipole coupling, which " *
            "depends on it",
        ),
    )
    n_scale = sqrt(sum(abs2, direction))
    iszero(n_scale) && throw(ArgumentError("Beam direction must be non-zero"))
    n̂ = direction ./ n_scale
    multipole == :E2 && return quadrupole_geometry(ε̂, n̂)
    multipole == :M1 ||
        throw(ArgumentError("Unknown multipole :$multipole; expected :E1, :E2 or :M1"))
    b̂ = SVector{3,ComplexF64}(
        n̂[2] * ε̂[3] - n̂[3] * ε̂[2],
        n̂[3] * ε̂[1] - n̂[1] * ε̂[3],
        n̂[1] * ε̂[2] - n̂[2] * ε̂[1],
    )
    if !isapprox(sum(abs2, b̂), 1; atol=1e-9)
        throw(
            ArgumentError(
                "Polarisation must be transverse to the beam direction " *
                "(ε ⋅ direction = 0) for a magnetic-dipole coupling, whose " *
                "magnetic polarisation direction × ε it fixes",
            ),
        )
    end
    dipole_geometry(b̂)
end

# Returns the complex coupling matrix over `basis` in the upper⟨row|lower⟩⟨col|
# block of the given level lists, as scale × amplitude(lower, upper) ×
# geometry[q]; `scale` fixes the units (1 for a relative matrix).
function assemble_couplings(
    basis::StateBasis,
    lower_levels,
    upper_levels,
    amplitude,
    rank,
    geometry,
    scale,
)
    num_states = length(basis)
    C = zeros(typeof(complex(1.0) * scale), num_states, num_states)
    for (i, lo) in enumerate(basis)
        lo.level in lower_levels || continue
        for (k, hi) in enumerate(basis)
            hi.level in upper_levels || continue
            q = hi.m - lo.m
            abs(q) <= rank || continue
            a = amplitude(lo, hi)
            iszero(a) && continue
            C[k, i] = scale * a * geometry[Int(q)+rank+1]
        end
    end
    C
end

"""
    rabi_frequency(species, lower::StateSpec => upper::StateSpec, intensity, ε, direction)
    rabi_frequency(s::Eigenstructure, lower::StateSpec => upper::StateSpec, intensity, ε, direction)

Returns the complex angular Rabi frequency ``Ω`` of a running wave with the
given intensity, polarisation `ε` and propagation direction. The Hamiltonian
coupling is ``(ħ/2)(Ω |u⟩⟨l| + \\mathrm{h.c.})``; an isolated resonant pair
has excitation probability ``\\sin^2(|Ω|t/2)``.

Following [James1998; Eqs. (5.11)–(5.13) in the arXiv numbering](@cite),

```math
Ω_\\mathrm{E1} = \\sqrt{\\frac{6π c^2 I A}{ħ ω^3}} \\, a_{ul} \\, d_{Δm}(\\hat ε),
\\qquad
Ω_\\mathrm{E2} = \\sqrt{\\frac{20π c^2 I A}{ħ ω^3}} \\, a_{ul} \\, Γ_{Δm},
```

where ``a_{ul}`` is [`transition_amplitude`](@ref), and ``d_q`` and ``Γ_q``
come from [`dipole_geometry`](@ref) and [`quadrupole_geometry`](@ref). A
magnetic-dipole (M1) pair (cf. [`Levels.multipole`](@ref)) couples to the
wave's magnetic field ``\\vec B = \\hat n × \\vec E / c``, so that

```math
Ω_\\mathrm{M1} = \\sqrt{\\frac{6π c^2 I A}{ħ ω^3}} \\, a_{ul} \\,
    d_{Δm}(\\hat n × \\hat ε),
```

the E1 prefactor following from ``A_\\mathrm{M1} = ω^3 |⟨J'‖μ‖J⟩|^2 /
(3π ε_0 ħ c^5 (2J'+1))``. The species form uses zero-field amplitudes; the
structure form includes field-induced ``F`` mixing. Both use the fine-structure
pair's Einstein A coefficient and centroid frequency ``ω`` for the physical
scale.

`ε` and `direction` are normalised internally; `direction` enters only for
E2 and M1, and an M1 coupling requires `ε` transverse to it. The complex
phase, including the phase of `ε`, is retained so matrices built from these
values preserve interference between components. `abs(Ω)` gives the
magnitude. [`dipole_components`](@ref) defines the interaction sign.
"""
function rabi_frequency(
    species::OneElectronSpecies,
    transition::Pair{<:StateSpec,<:StateSpec},
    intensity,
    ε,
    direction,
)
    lower = validate_state_kind(species, transition.first)
    upper = validate_state_kind(species, transition.second)
    rabi_from_amplitude(
        species,
        lower,
        upper,
        intensity,
        ε,
        direction,
        transition_amplitude(species, lower => upper),
    )
end

function rabi_frequency(
    s::Eigenstructure,
    transition::Pair{<:StateSpec,<:StateSpec},
    intensity,
    ε,
    direction,
)
    m_lo, i = resolve_state(s, transition.first)
    m_hi, k = resolve_state(s, transition.second)
    rabi_from_amplitude(
        s.species,
        m_lo.basis[i],
        m_hi.basis[k],
        intensity,
        ε,
        direction,
        transition_amplitude(s, transition),
    )
end

# The James-formula core shared by the rabi_frequency methods: scales the given
# relative angular amplitude to the absolute Rabi frequency for the beam
# intensity and geometry.
function rabi_from_amplitude(
    species,
    lower::StateSpec,
    upper::StateSpec,
    intensity,
    ε,
    direction,
    amplitude,
)
    fs_lo = fine_structure(lower.level)
    fs_hi = fine_structure(upper.level)
    m = multipole(species, fs_lo => fs_hi)
    rank = multipole_rank(m)
    scale = rabi_scale(species, fs_lo, fs_hi, intensity)
    geometry = channel_geometry(m, ε, direction)
    Δm = upper.m - lower.m
    g = abs(Δm) <= rank ? geometry[Int(Δm)+rank+1] : 0.0im
    uconvert(u"rad/µs", scale * (amplitude * g))
end

"""
    coupling_matrix(species, basis::StateBasis, lower => upper, intensity, ε, direction)
    coupling_matrix(s::Eigenstructure, lower => upper, intensity, ε, direction)

Returns the complex angular Rabi frequencies of [`rabi_frequency`](@ref) in
the upper-row, lower-column block, zero elsewhere. The Hamiltonian in
angular-frequency units is `(C + C') / 2`.

The species form uses the supplied coupled basis. The structure form uses
`basis(s)` in the field eigenbasis, with complete-manifold amplitudes rotated
before any selection of states. For hyperfine species each level argument
may name a whole fine-structure manifold or select one ``F`` label after
rotation. Partial and reordered bases retain the corresponding matrix elements.
"""
function coupling_matrix(
    species::OneElectronSpecies,
    basis::StateBasis,
    pair::Pair,
    intensity,
    ε,
    direction,
)
    validate_basis_kind(species, basis)
    p = level_pair(species, pair)
    assemble_couplings(
        basis,
        p.lower_levels,
        p.upper_levels,
        (lo, hi) -> zero_field_amplitude(species, lo, hi, p.rank),
        p.rank,
        channel_geometry(p.multipole, ε, direction),
        rabi_scale(species, p.fs_lo, p.fs_hi, intensity),
    )
end

function coupling_matrix(s::Eigenstructure, pair::Pair, intensity, ε, direction)
    p = level_pair(s.species, pair)
    assemble_couplings(
        s.basis,
        p.lower_levels,
        p.upper_levels,
        amplitude_lookup(s, p.fs_lo, p.fs_hi, p.rank),
        p.rank,
        channel_geometry(p.multipole, ε, direction),
        rabi_scale(s.species, p.fs_lo, p.fs_hi, intensity),
    )
end

"""
    relative_couplings(species, basis::StateBasis, lower => upper, ε, direction)
    relative_couplings(s::Eigenstructure, lower => upper, ε, direction)

Returns the dimensionless complex E1, E2 or M1 coupling amplitudes
``c_{m→m'} = ⟨\\mathrm{up}|T^k_q|\\mathrm{lo}⟩_\\mathrm{rel} \\, g_{Δm}`` — the
[`transition_amplitude`](@ref) times the normalised geometric channel
amplitude ([`dipole_geometry`](@ref) or [`quadrupole_geometry`](@ref)) — in
the upper⟨row|lower⟩⟨col| block, zero elsewhere: [`coupling_matrix`](@ref)
without its physical scale (an M1 pair takes the dipole geometry of the magnetic
polarisation `direction × ε`, cf. [`rabi_frequency`](@ref)). The species form is
the zero-field matrix in the coupled basis, the structure form the at-field
matrix in the field eigenbasis of `s` indexed by `basis(s)` (rotated as a whole
fine-structure pair, then selected, like `coupling_matrix`). The overall scale
is typically fixed by normalising one probed component to a known carrier Rabi
frequency with [`rabi_normalised`](@ref).
"""
function relative_couplings(
    species::OneElectronSpecies,
    basis::StateBasis,
    pair::Pair,
    ε,
    direction,
)
    validate_basis_kind(species, basis)
    p = level_pair(species, pair)
    assemble_couplings(
        basis,
        p.lower_levels,
        p.upper_levels,
        (lo, hi) -> zero_field_amplitude(species, lo, hi, p.rank),
        p.rank,
        channel_geometry(p.multipole, ε, direction),
        1.0,
    )
end

function relative_couplings(s::Eigenstructure, pair::Pair, ε, direction)
    p = level_pair(s.species, pair)
    assemble_couplings(
        s.basis,
        p.lower_levels,
        p.upper_levels,
        amplitude_lookup(s, p.fs_lo, p.fs_hi, p.rank),
        p.rank,
        channel_geometry(p.multipole, ε, direction),
        1.0,
    )
end

"""
    rabi_normalised(couplings, basis::StateBasis, transition::Pair, Ω0)

Scales the given relative coupling matrix such that the entry for the
`lower => upper` `transition` has magnitude `abs(Ω0)` (its carrier Rabi
frequency, an angular frequency; a complex [`rabi_frequency`](@ref) contributes its
magnitude only), fixing the physical scale of all the couplings. A pure index
operation: the matrix keeps its representation (pass `basis(s)` for a
structure-form matrix), and its internally consistent relative phases.

An error is raised if the amplitude for the given transition (nearly) vanishes in
`couplings`, as the normalisation would then be ill-defined.
"""
function rabi_normalised(
    couplings::AbstractMatrix,
    basis::StateBasis,
    transition::Pair,
    Ω0,
)
    c = couplings[
        stateindex(basis, transition.second),
        stateindex(basis, transition.first),
    ]
    if abs(c) <= 1e-9 * maximum(abs, couplings)
        throw(
            ArgumentError(
                "Coupling for transition $(transition.first) => $(transition.second) " *
                "(nearly) vanishes for the given geometry",
            ),
        )
    end
    (abs(Ω0) / abs(c)) .* couplings
end

"""
    intensity_for_rabi(species, lower::StateSpec => upper::StateSpec, Ω, ε, direction)
    intensity_for_rabi(s::Eigenstructure, lower::StateSpec => upper::StateSpec, Ω, ε, direction)

Returns the intensity at which a running wave of polarisation `ε` and
propagation direction `direction` drives the `lower` → `upper` component with the
Rabi frequency magnitude `abs(Ω)` (an angular frequency) — the inverse of
[`rabi_frequency`](@ref), whose arguments and species/structure forms it
shares, since the Rabi frequency scales as the square root of the intensity.
Raises if the geometry does not drive the component at all.
"""
function intensity_for_rabi(
    atom::Union{OneElectronSpecies,Eigenstructure},
    transition::Pair{<:StateSpec,<:StateSpec},
    Ω,
    ε,
    direction,
)
    lower, upper = transition
    reference = 1.0u"W/m^2"
    Ω_ref = abs(rabi_frequency(atom, lower => upper, reference, ε, direction))
    iszero(Ω_ref) && throw(
        ArgumentError(
            "The component $lower → $upper is not driven by this polarisation " *
            "and beam direction",
        ),
    )
    uconvert(u"W/m^2", reference * abs2(uconvert(NoUnits, abs(Ω) / Ω_ref)))
end

# The reduced electric-dipole matrix element per unit Clebsch–Gordan amplitude,
# ⟨J'‖d‖J⟩/√(2J'+1) = √(3π ε₀ ħ c³ A / ω³), from the Einstein A coefficient
# and the centroid frequency of the pair (the inverse of the ElectricDipole
# conversion of the species constructors).
function reduced_dipole_scale(species, fs_lo, fs_hi)
    multipole(species, fs_lo => fs_hi) == :E1 || throw(
        ArgumentError(
            "'$fs_lo' → '$fs_hi' is not an electric-dipole transition; " *
            "dipole_components covers E1 only",
        ),
    )
    a = einstein_a(species, fs_lo => fs_hi)
    isnothing(a) &&
        throw(ArgumentError("No known transition between '$fs_lo' and '$fs_hi'"))
    ω = transition_frequency(species, fs_lo => fs_hi)
    uconvert(u"C*m", sqrt(3π * u"ε0" * u"ħ" * u"c"^3 * a / ω^3))
end

"""
    dipole_components(species, basis::StateBasis, lower => upper)
    dipole_components(s::Eigenstructure, lower => upper)

Returns three E1 electric-dipole operator matrices ``d_q`` in coulomb metres,
indexed by `q + 2` for ``q = -1, 0, +1``. Entries occupy the upper-row,
lower-column block; all others are zero. Basis and ``F``-label selection
follow [`coupling_matrix`](@ref). An E2 or M1 pair raises an error.

The reduced element is defined by

```math
⟨J' m'|d_q|J m⟩ = ⟨J m; 1 q|J' m'⟩
    \\frac{⟨J'‖d‖J⟩}{\\sqrt{2J'+1}}, \\qquad
\\frac{⟨J'‖d‖J⟩}{\\sqrt{2J'+1}}
    = \\sqrt{\\frac{3π ε_0 ħ c^3 A}{ω^3}},
```

with the fine-structure Einstein A coefficient and centroid frequency ``ω``.
The reduced scale is taken positive, and each entry is that scale times
[`transition_amplitude`](@ref). Rate data supplies no relative sign between
the reduced elements of different pairs; coherent sums through several
intermediate manifolds require a consistent choice of those signs.

For ``\\vec E(t) = \\mathrm{Re}[\\vec ε E_0 e^{-iωt}]``, with unit ``\\vec ε`` and
``E_0 = \\sqrt{2I/(ε_0 c)}``, the interaction ``H_\\mathrm{int} = -\\vec d ⋅ \\vec E``
gives the rotating-wave amplitude

```math
Ω_\\mathrm{int} = -\\frac{E_0}{ħ}\\sum_q g_q(ε)⟨u|d_q|l⟩,
\\qquad g_q(ε) = (-1)^q ε_{-q}.
```

[`rabi_frequency`](@ref) and [`coupling_matrix`](@ref) use the contraction
with the overall minus sign omitted, a common phase shift of the field.
For unit `ε`, their E1 entries therefore equal
`E0/ħ * sum(dipole_geometry(ε)[q+2] * d[q+2][u, l] for q in -1:1)`.
"""
function dipole_components(species::OneElectronSpecies, basis::StateBasis, pair::Pair)
    validate_basis_kind(species, basis)
    p = level_pair(species, pair)
    scale = reduced_dipole_scale(species, p.fs_lo, p.fs_hi)
    amplitude = (lo, hi) -> zero_field_amplitude(species, lo, hi, 1)
    assemble_components(basis, p.lower_levels, p.upper_levels, amplitude, scale)
end

function dipole_components(s::Eigenstructure, pair::Pair)
    p = level_pair(s.species, pair)
    scale = reduced_dipole_scale(s.species, p.fs_lo, p.fs_hi)
    amplitude = amplitude_lookup(s, p.fs_lo, p.fs_hi, 1)
    assemble_components(s.basis, p.lower_levels, p.upper_levels, amplitude, scale)
end

function assemble_components(basis, lower_levels, upper_levels, amplitude, scale)
    num_states = length(basis)
    components = [zeros(num_states, num_states) .* scale for _ in -1:1]
    for (i, lo) in enumerate(basis)
        lo.level in lower_levels || continue
        for (k, hi) in enumerate(basis)
            hi.level in upper_levels || continue
            q = hi.m - lo.m
            abs(q) <= 1 || continue
            components[Int(q)+2][k, i] = scale * amplitude(lo, hi)
        end
    end
    components
end

export transition_amplitude,
    rabi_frequency,
    coupling_matrix,
    rabi_normalised,
    intensity_for_rabi,
    relative_couplings,
    dipole_components
public clebsch_gordan,
    multipole,
    multipole_rank,
    parity_multipole,
    hyperfine_reduction,
    relative_amplitudes
