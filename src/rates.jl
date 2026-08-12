using Unitful
using WignerSymbols

"""
    multipole_rank(lower, upper)

Returns the electric-multipole rank of the transition between the two given
levels — 1 (E1) if the level parities differ, 2 (E2) if they are equal — as
determined from the orbital angular momenta. Magnetic multipoles and electric
orders beyond E2 are not considered.
"""
function multipole_rank(lower, upper)
    lo = convert(NoHyperfineNumberSpec, lower)
    hi = convert(NoHyperfineNumberSpec, upper)
    isodd(lo.l + hi.l) ? 1 : 2
end

"""
    clebsch_gordan(lower::StateSpec, upper::StateSpec)

Returns the bare Clebsch–Gordan coefficient ``⟨j m; R Δm | j' m'⟩`` between two
fine-structure states, with the rank ``R`` the electric-multipole order from
[`multipole_rank`](@ref).

This is the purely geometric coupling factor; the general transition amplitude
relative to the reduced matrix element — species-aware, covering hyperfine
``F`` levels and their exact mixing at a static field — is
[`transition_amplitude`](@ref).
"""
function clebsch_gordan(lower::StateSpec, upper::StateSpec)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    rank = multipole_rank(lo, hi)
    Δm = upper.m - lower.m
    abs(Δm) <= rank || return 0.0
    Float64(clebschgordan(lo.j, lower.m, rank, Δm, hi.j, upper.m))
end

# Canonicalises a state specification and checks it names a hyperfine (F) level.
function hyperfine_state(state::StateSpec)
    spec = parse_level(state.level)
    if !(spec isa HyperfineNumberSpec)
        throw(
            ArgumentError(
                "State '$state' must specify a hyperfine (F) level " *
                "for a hyperfine species",
            ),
        )
    end
    StateSpec(spec, state.m)
end

"""
    transition_amplitude(species, lower::StateSpec, upper::StateSpec[, B]; rank)
    transition_amplitude(m_lower::HyperfineManifold, m_upper::HyperfineManifold,
                         lower::StateSpec, upper::StateSpec; rank)

Returns the amplitude of the rank-`rank` multipole transition between the two
given states relative to the reduced matrix element of the *fine-structure*
transition, ``⟨\\mathrm{up}|T^k_q|\\mathrm{lo}⟩ / (⟨J'‖T^k‖J⟩/\\sqrt{2J'+1})``
with ``q = Δm``. `rank` defaults to the electric-multipole order of
[`multipole_rank`](@ref); pass it explicitly for other operators (e.g. `1` for
M1 within one fine-structure level).

For fine-structure states this is the bare Clebsch–Gordan coefficient
``⟨j m; k Δm | j' m'⟩`` (cf. [`Levels.clebsch_gordan`](@ref)); for hyperfine
states of a [`HyperfineOneElectronSpecies`](@ref) at zero field it is
``⟨F m; k Δm | F' m'⟩ β^{(k)}(F → F')`` with the
[`Levels.hyperfine_reduction`](@ref) factor ``β``, which in the ``I → 0``
limit degenerates exactly (including sign) to the fine-structure form.

With the static flux density `B` (along the quantisation axis ẑ) as the final
positional argument, hyperfine states denote the **adiabatically-labelled
eigenstates** of the hyperfine + Zeeman Hamiltonian (cf.
[`hyperfine_manifold`](@ref)) and the amplitude is exact at that field,
including the ``F`` mixing within the manifolds (which modifies the ⁴³Ca⁺
D``_{5/2}`` component amplitudes at the few-percent level already at 0.5 mT).
``m_F`` stays exact, so each state pair still couples through the single
spherical component ``q = Δm``. For a fine-structure species a static field
along ẑ leaves the ``m_J`` eigenstates unchanged, so `B` is accepted and
ignored — the zero-field form is already exact — keeping generic code uniform
across species.

The manifold-pair form evaluates the same at-field amplitude from pre-solved
[`HyperfineManifold`](@ref)s (which must share one species and field; pass the
same manifold twice for a transition within one manifold), avoiding the
re-diagonalisation in loops over many components; the matrix counterpart is
the [`eigenbasis_transform`](@ref) rotation. At-field amplitudes connect
continuously (including sign) to the zero-field ones at low field, but the
individual signs inherit the eigenvector convention of
[`hyperfine_manifold`](@ref); magnitudes, and relative phases within one
consistently-rotated matrix, are the physically meaningful quantities.
"""
function transition_amplitude(species, lower::StateSpec, upper::StateSpec; rank=nothing)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    k = something(rank, multipole_rank(lo, hi))
    Δm = upper.m - lower.m
    (abs(Δm) <= k && is_triangle(lo.j, k, hi.j)) || return 0.0
    Float64(clebschgordan(lo.j, lower.m, k, Δm, hi.j, upper.m))
end

transition_amplitude(species, lower::StateSpec, upper::StateSpec, B; rank=nothing) =
    transition_amplitude(species, lower, upper; rank)

function transition_amplitude(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec,
    upper::StateSpec;
    rank=nothing,
)
    lo = hyperfine_state(lower)
    hi = hyperfine_state(upper)
    validate_hyperfine(species, lo.level)
    validate_hyperfine(species, hi.level)
    k = something(
        rank,
        multipole_rank(fine_structure(lo.level), fine_structure(hi.level)),
    )
    Δm = hi.m - lo.m
    abs(Δm) <= k || return 0.0
    β = hyperfine_reduction(species.nuclear_spin, lo.level, hi.level; rank=k)
    iszero(β) && return 0.0
    Float64(clebschgordan(lo.level.f, lo.m, k, Δm, hi.level.f, hi.m)) * β
end

function transition_amplitude(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec,
    upper::StateSpec,
    B;
    rank=nothing,
)
    iszero(B) && return transition_amplitude(species, lower, upper; rank)
    lo = hyperfine_state(lower)
    hi = hyperfine_state(upper)
    validate_hyperfine(species, lo.level)
    validate_hyperfine(species, hi.level)
    fs_lo = fine_structure(lo.level)
    fs_hi = fine_structure(hi.level)
    m_lower = hyperfine_manifold(species, fs_lo, B)
    m_upper = fs_hi == fs_lo ? m_lower : hyperfine_manifold(species, fs_hi, B)
    transition_amplitude(m_lower, m_upper, lo, hi; rank)
end

function transition_amplitude(
    m_lower::HyperfineManifold,
    m_upper::HyperfineManifold,
    lower::StateSpec,
    upper::StateSpec;
    rank=nothing,
)
    if m_lower.species !== m_upper.species
        throw(ArgumentError("Manifold solutions must belong to one species"))
    end
    if m_lower.field != m_upper.field
        throw(
            ArgumentError(
                "Manifold solutions must share one static field, " *
                "got $(m_lower.field) and $(m_upper.field)",
            ),
        )
    end
    lo = hyperfine_state(lower)
    hi = hyperfine_state(upper)
    k = something(rank, multipole_rank(m_lower.level, m_upper.level))
    q = hi.m - lo.m
    abs(q) <= k || return 0.0

    # m_F is exact, so both eigenvectors live in single-m_F blocks and only the
    # F decompositions are summed over.
    v_lo = view(m_lower.states, :, stateindex(m_lower.basis, lo))
    v_hi = view(m_upper.states, :, stateindex(m_upper.basis, hi))
    amplitude = 0.0
    for (i, sl) in enumerate(m_lower.basis)
        (sl.m == lo.m && !iszero(v_lo[i])) || continue
        for (j, su) in enumerate(m_upper.basis)
            (su.m == hi.m && !iszero(v_hi[j])) || continue
            β = hyperfine_reduction(
                m_lower.species.nuclear_spin,
                sl.level,
                su.level;
                rank=k,
            )
            iszero(β) && continue
            amplitude +=
                v_hi[j] *
                v_lo[i] *
                Float64(clebschgordan(sl.level.f, sl.m, k, q, su.level.f, su.m)) *
                β
        end
    end
    amplitude
end

# Triangle condition (|a - b| ≤ c ≤ a + b with integer perimeter) for angular
# momentum coupling.
is_triangle(a, b, c) = abs(a - b) <= c <= a + b && isinteger(a + b - c)

"""
    hyperfine_reduction(nuclear_spin, lower, upper; rank)
    hyperfine_reduction(species::HyperfineOneElectronSpecies, lower, upper; rank)

Returns the hyperfine reduction factor ``β^{(k)}(F → F')`` relating the
amplitude of a rank-`k` electronic multipole transition between hyperfine
levels to that between their fine-structure levels,

```math
β^{(k)}(F → F') = (-1)^{J + I + F' + k}
    \\sqrt{(2F + 1)(2J' + 1)}
    \\begin{Bmatrix} J & J' & k \\\\ F' & F & I \\end{Bmatrix},
```

such that ``⟨F' m'|T^k_q|F m⟩ / (⟨J'‖T‖J⟩/\\sqrt{2J'+1}) =
⟨F m; k q|F' m'⟩ β^{(k)}``. The phase convention matches the coupled basis of
[`Levels.coupling_transform`](@ref) (``⟨I m_I; J m_J | F m_F⟩``, nuclear spin
first); in the ``I → 0`` limit ``β ≡ +1`` exactly.

The squares are the relative line strengths: ``\\sum_F β^2 = 1`` for fixed
``F'`` (each hyperfine sublevel decays at the full fine-structure rate, with
``β^2`` the branching fractions), which is why F-resolved rates must never be
entered as separate Einstein A coefficients.

`rank` defaults to the electric-multipole order of [`multipole_rank`](@ref);
pass it explicitly for other operators (e.g. `1` for M1 within one
fine-structure level).
"""
function hyperfine_reduction(nuclear_spin, lower, upper; rank=nothing)
    lo = parse_level(lower)
    hi = parse_level(upper)
    if !(lo isa HyperfineNumberSpec && hi isa HyperfineNumberSpec)
        throw(ArgumentError("Levels must specify hyperfine (F) levels"))
    end
    k = something(rank, multipole_rank(fine_structure(lo), fine_structure(hi)))
    if !is_triangle(lo.f, k, hi.f) || !is_triangle(lo.j, k, hi.j)
        return 0.0
    end
    Float64(
        (-1)^Int(lo.j + nuclear_spin + hi.f + k) *
        sqrt((2 * lo.f + 1) * (2 * hi.j + 1)) *
        wigner6j(lo.j, hi.j, k, hi.f, lo.f, nuclear_spin),
    )
end

function hyperfine_reduction(
    species::HyperfineOneElectronSpecies,
    lower,
    upper;
    rank=nothing,
)
    lo = parse_level(lower)
    hi = parse_level(upper)
    lo isa HyperfineNumberSpec && validate_hyperfine(species, lo)
    hi isa HyperfineNumberSpec && validate_hyperfine(species, hi)
    hyperfine_reduction(species.nuclear_spin, lo, hi; rank)
end

"""
    rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity, ε, n)
    rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity, ε, n, B)

Returns the Rabi frequency (in angular units, such that the excitation
probability oscillates as ``\\sin^2(Ω t / 2)`` on resonance) for driving the
given transition with a running wave of the given intensity, polarisation `ε`
and propagation direction `n` (cf. [`beam_vectors`](@ref); `n` only enters for
quadrupole transitions).

The coupling strength is obtained from the Einstein A coefficient of the
transition following `[James1998]`,

``Ω = \\sqrt{6 π c^2 I A / (ħ ω^3)} \\, |⟨j m; 1 Δm | j' m'⟩ \\, d_{Δm}|``

for electric-dipole and

``Ω = \\sqrt{20 π c^2 I A / (ħ ω^3)} \\, |⟨j m; 2 Δm | j' m'⟩ \\, Γ_{Δm}|``

for electric-quadrupole transitions, with ``d_q`` and ``Γ_q`` the geometric
channel amplitudes from [`dipole_geometry`](@ref) and
[`quadrupole_geometry`](@ref).

Only the directions of `ε` and `n` matter — both are normalised internally, as
the field amplitude is fixed by the intensity.

For hyperfine states of a [`HyperfineOneElectronSpecies`](@ref), the Einstein A
coefficient is that of the fine-structure transition, the angular factor is
the zero-field [`transition_amplitude`](@ref) (``F``-basis Clebsch–Gordan
coefficient times the [`Levels.hyperfine_reduction`](@ref) factor), and the
transition frequency includes the zero-field hyperfine shifts. At finite field
the ``F`` mixing within the manifolds modifies the component amplitudes (at
the few-percent level for the ⁴³Ca⁺ D``_{5/2}`` manifold at 0.5 mT); passing
the static flux density `B` (along the quantisation axis ẑ) as the final
positional argument evaluates the angular factor exactly at that field, with
the states denoting the adiabatically-labelled eigenstates (cf.
[`transition_amplitude`](@ref)). The ``ω^3`` prefactor keeps the zero-field
transition frequency, whose Zeeman corrections are fractionally ``~10^{-9}``
on an optical transition. For a fine-structure species the `B` form returns
the plain result, which is already exact.

# References

- `[James1998]`: D. F. V. James, "Quantum dynamics of cold trapped ions with
  application to quantum computation", Appl. Phys. B **66**, 181 (1998),
  [doi:10.1007/s003400050373](https://doi.org/10.1007/s003400050373); Eqs.
  (5.11)–(5.13) in the
  [arXiv:quant-ph/9702053](https://arxiv.org/abs/quant-ph/9702053) numbering.
"""
function rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity, ε, n)
    rabi_from_amplitude(
        species,
        lower,
        upper,
        intensity,
        ε,
        n,
        transition_amplitude(species, lower, upper),
    )
end

function rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity, ε, n, B)
    rabi_from_amplitude(
        species,
        lower,
        upper,
        intensity,
        ε,
        n,
        transition_amplitude(species, lower, upper, B),
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
    n,
    angular,
)
    lo = parse_level(lower.level)
    hi = parse_level(upper.level)
    fs_lo = fine_structure(lo)
    fs_hi = fine_structure(hi)
    a = einstein_a(species, fs_lo, fs_hi)
    isnothing(a) &&
        throw(ArgumentError("No known transition between '$fs_lo' and '$fs_hi'"))
    ω = transition_frequency(species, lo, hi)
    rank = multipole_rank(fs_lo, fs_hi)
    Δm = upper.m - lower.m
    ε_scale = sqrt(sum(abs2, ε))
    n_scale = sqrt(sum(abs2, n))
    if iszero(ε_scale) || (rank == 2 && iszero(n_scale))
        throw(ArgumentError("Polarisation and beam direction must be non-zero"))
    end
    prefactor, geometry = if abs(Δm) > rank
        0.0, 0.0im
    elseif rank == 1
        6.0, dipole_geometry(ε)[Int(Δm)+2] / ε_scale
    else
        20.0, quadrupole_geometry(ε, n)[Int(Δm)+3] / (ε_scale * n_scale)
    end
    scale = prefactor * π * u"c"^2 * intensity * a / (u"ħ" * ω^3)
    uconvert(u"µs^-1", sqrt(scale) * abs(angular * geometry))
end

"""
    rabi_normalised(couplings, basis::StateBasis, transition::Pair, Ω0)

Scales the given relative coupling matrix such that the entry for the
`lower => upper` `transition` has magnitude `Ω0` (its carrier Rabi frequency, in
angular units), fixing the physical scale of all the couplings.

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
    (Ω0 / abs(c)) .* couplings
end

export transition_amplitude, rabi_frequency, rabi_normalised
public clebsch_gordan, multipole_rank, hyperfine_reduction
