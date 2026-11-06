using Unitful
using WignerSymbols

"""
Returns ``⟨j m; 1 q | j' m + q⟩``, or zero where that is not a valid
electric-dipole channel.
"""
function dipole_cg(j, m, q, j_upper)
    abs(j - j_upper) <= 1 || return 0.0
    abs(m + q) <= j_upper || return 0.0
    Float64(clebschgordan(j, m, 1, q, j_upper, m + q))
end

# Diagonal electronic rank-two operator in the |J, m_J⟩ basis (J_z in ħ).
tensor_projection(j, m) = j > 1//2 ? (3m^2 - j * (j + 1)) / (j * (2j - 1)) : 0//1

"""
Returns one E1 channel's scalar polarisability at photon energy `ħω`:
``2Δ d^2 / [3(2j+1)(Δ^2-(ħω)^2)]``, for source angular momentum `j`, reduced
dipole `d` and signed intermediate-level energy difference `Δ`.

Also supplies the static contribution subtracted when resolving an
[`ImplicitPolarisability`](@ref), keeping the anchor and evaluation consistent.
"""
channel_scalar_polarisability(d, j, Δ, ħω) =
    uconvert(u"C*m^2/V", d^2 / (3 * (2j + 1)) * 2Δ / (Δ^2 - ħω^2))

"""
Returns one E1 channel's tensor polarisability for source angular momentum
`j > 1/2` and intermediate angular momentum `j_channel`. Other arguments are
as for [`Levels.channel_scalar_polarisability`](@ref).

Subtracts ``α_0`` from the stretched state's π-polarised response
``α_0 + α_2``.
"""
function channel_tensor_polarisability(d, j, j_channel, Δ, ħω)
    weight = dipole_cg(j, j, 0, j_channel)^2 / (2 * j_channel + 1) - 1 / (3 * (2j + 1))
    uconvert(u"C*m^2/V", d^2 * weight * 2Δ / (Δ^2 - ħω^2))
end

"""
    polarisation_weights(ε)

Returns the fractions ``w_q = |ε_{-q}|^2`` of the beam intensity driving each of
the ``Δm = q`` electric-dipole channels, as a length-3 tuple indexed by
`q + 2`.

These are the squared [`dipole_geometry`](@ref) amplitudes, normalised to
sum to one. `ε` must be nonzero; its scale does not affect the result.
"""
function polarisation_weights(ε)
    w = ntuple(i -> abs2(spherical_component(ε, 2 - i)), 3)
    total = sum(w)
    if total <= 0
        throw(ArgumentError("Polarisation vector must be non-zero"))
    end
    w ./ total
end

"""
    quadrupole_weights(ε, direction)

Returns the weights ``|Γ_q|^2`` with which the beam intensity drives each of the
``Δm = q`` electric-quadrupole channels, as a length-5 tuple indexed by
`q + 3`.

Uses [`quadrupole_geometry`](@ref) with normalised `ε` and `direction`. Both must be
nonzero and transverse, ``ε ⋅ direction = 0`` (tolerance ``10^{-3}`` after
normalisation). The weights need not sum to one: only the rank-two part of
``ε ⊗ direction`` couples.
"""
function quadrupole_weights(ε, direction)
    ε_scale = sqrt(sum(abs2, ε))
    n_scale = sqrt(sum(abs2, direction))
    if iszero(ε_scale) || iszero(n_scale)
        throw(ArgumentError("Polarisation and beam direction must be non-zero"))
    end
    ε_unit = ε ./ ε_scale
    n_unit = direction ./ n_scale
    if abs(sum(ε_unit .* n_unit)) > 1e-3
        throw(
            ArgumentError(
                "Polarisation must be transverse to the beam direction (ε ⋅ direction = 0)",
            ),
        )
    end
    Γ = quadrupole_geometry(ε_unit, n_unit)
    ntuple(i -> abs2(Γ[i]), 5)
end

"""
Returns a cache of manifold solutions for the intermediate levels of the
polarisability sums, seeded with the structure's own solutions.
"""
manifold_cache(s::Eigenstructure{S,L}) where {S,L} =
    Dict{NoHyperfineNumberSpec,ManifoldSolution{L}}(s.manifolds)

"""
Returns the solution of an intermediate manifold at the structure's field,
solving and caching it on first use. A hyperfine intermediate without
hyperfine constants stays degenerate at its centroid in the coupled basis (no
hyperfine or Zeeman structure) — the documented far-detuned approximation,
which [`Levels.validate_channel_manifold`](@ref) refuses for near-resonant
partners.
"""
function intermediate_solution(s::Eigenstructure, fs, cache)
    get!(cache, fs) do
        solve_manifold(s.species, fs, s.B; require_constants=false)
    end
end

"""
Requires the hyperfine constants of a near-resonant channel partner, whose
field solution must be exact.
"""
validate_channel_manifold(species, fs) = nothing

function validate_channel_manifold(species::HyperfineOneElectronSpecies, fs)
    if !haskey(species.hyperfine, fs)
        throw(
            ArgumentError(
                "No hyperfine constants known for level '$fs', whose " *
                "eigenstates at a static field therefore cannot be resolved",
            ),
        )
    end
end

"""
    state_polarisabilities(s::Eigenstructure, state::StateSpec, ħω)

Returns the electric-dipole polarisability of the labelled eigenstate `state`
of `s` for each of the three driving polarisations at photon energy `ħω`, as
a length-3 vector indexed by `q + 2`. Contracts with
[`polarisation_weights`](@ref); cross terms vanish because every state has
definite ``m``. The public interface is [`polarisability`](@ref).
"""
function state_polarisabilities(s::Eigenstructure, state::StateSpec, ħω)
    m, index = resolve_state(s, state)
    if isnothing(level_polarisability(s.species, m.level))
        throw(ArgumentError("No light-shift data known for level '$(m.level)'"))
    end
    α = manifold_polarisabilities(s, m, ħω, nothing, manifold_cache(s))
    α[index, :]
end

"""
Returns the static remainder's tensor factor in each eigenstate:
``⟨(3 J_z^2 - J(J+1)) / (J(2J-1))⟩``, with ``J_z`` in units of ``ħ``;
zero for ``J ≤ 1/2``.

The operator is diagonal in the ``|m_J, m_I⟩`` product basis. Evaluating it
there retains coherences between ``F`` components that a weighted average of
coupled-state expectations would lose.
"""
tensor_expectations(species, m::ManifoldSolution) =
    [Float64(tensor_projection(m.level.j, s.m)) for s in m.basis]

function tensor_expectations(species::HyperfineOneElectronSpecies, m::ManifoldSolution)
    j, i_nuc = m.level.j, species.nuclear_spin
    diagonal =
        [Float64(tensor_projection(j, m_j)) for m_j in (-j):j for _ in (-i_nuc):i_nuc]
    products = coupling_transform(species, m.level) * m.states
    vec(diagonal' * abs2.(products))
end

"""
Returns a `length(source.basis) × 3` matrix of E1 polarisabilities of the
eigenstates of the `source` manifold of `s`, indexed by `[state, q + 2]` like
[`polarisation_weights`](@ref). Each intermediate state contributes
``|d_{ts}|^2 [w_q/(E_t-E_s-ħω) + w_{-q}/(E_t-E_s+ħω)]``, with dipole matrix
element ``d_{ts}`` and energies at the structure's field. Adds the static
scalar remainder and the tensor remainder weighted by
[`Levels.tensor_expectations`](@ref).

`exclude` names an explicit intermediate level whose near-resonant term
moves into [`LightShiftCoefficients`](@ref)'s channel sum. The
counter-rotating term stays here; a lumped remainder cannot be split.
"""
function manifold_polarisabilities(
    s::Eigenstructure,
    source::ManifoldSolution,
    ħω,
    exclude,
    cache,
)
    species = s.species
    fs = source.level
    data = level_polarisability(species, fs)
    if !isnothing(exclude) && !haskey(data.reduced_dipoles, exclude)
        throw(
            ArgumentError(
                "The '$fs' → '$exclude' channel is lumped into the static " *
                "remainder of the polarisability data, so its near-resonant part " *
                "cannot be separated out",
            ),
        )
    end
    num_states = length(source.basis)
    # Uniform units avoid conversions inside the channel sum.
    α = fill(0.0u"C^2*m^2/J", num_states, 3)
    ħω_J = uconvert(u"J", ħω)
    e_source = uconvert.(u"J", u"ħ" .* source.energies)
    for (upper, d) in data.reduced_dipoles
        target = intermediate_solution(s, upper, cache)
        amplitudes = relative_amplitudes(species, source, target, 1)
        strength = uconvert(u"C^2*m^2", d^2 / (2 * upper.j + 1))
        centroid = uconvert(u"J", species.energies[upper] - species.energies[fs])
        e_target = uconvert.(u"J", u"ħ" .* target.energies)
        # Only the Δm = q = −1, 0, 1 partners couple; within a level's block
        # the canonical bases run over m in ascending order.
        for level in target.basis.levels
            block = staterange(target.basis, level)
            f_target = momentum(level)
            for (i, state) in enumerate(source.basis), q in -1:1
                m_target = state.m + q
                abs(m_target) <= f_target || continue
                k = first(block) + Int(m_target + f_target)
                @assert target.basis[k].m == m_target
                weight = strength * abs2(amplitudes[k, i])
                iszero(weight) && continue
                Δ = centroid + e_target[k] - e_source[i]
                # For a lower intermediate level (Δ < 0), the near-resonant
                # denominator is Δ + ħω.
                drop_rotating = upper == exclude && Δ > zero(Δ)
                drop_counter = upper == exclude && !drop_rotating
                if !drop_rotating
                    α[i, q+2] += weight / (Δ - ħω_J)
                end
                if !drop_counter
                    α[i, 2-q] += weight / (Δ + ħω_J)
                end
            end
        end
    end
    α = uconvert.(u"C*m^2/V", α)
    tensor = tensor_expectations(species, source)
    for i in 1:num_states, q in -1:1
        α[i, q+2] +=
            data.static_scalar +
            data.static_tensor * ((3 * (q == 0) - 1) / 2) * tensor[i]
    end
    α
end

"""
Returns the background polarisabilities of the basis states of `s` as
length-3 rows, evaluating each manifold once. Levels without data have `NaN`
rows. `exclude(level)` supplies each level's excluded intermediate partner.
"""
function polarisability_rows(s::Eigenstructure, ħω, exclude, cache)
    species = s.species
    blocks = Dict{NoHyperfineNumberSpec,Matrix{typeof(1.0u"C*m^2/V")}}()
    for fs in manifold_levels(s)
        isnothing(level_polarisability(species, fs)) && continue
        blocks[fs] =
            manifold_polarisabilities(s, s.manifolds[fs], ħω, exclude(fs), cache)
    end
    [
        begin
            fs = fine_structure(state.level)
            haskey(blocks, fs) ?
            blocks[fs][stateindex(s.manifolds[fs].basis, state), :] :
            fill(NaN * u"C*m^2/V", 3)
        end for state in s.basis
    ]
end

"""
One near-resonant coupling channel of a basis state: a Zeeman component
connecting it to `partner` in the other manifold of a level pair the laser sits
close to.

The signed `weight` is such that the second-order shift of the state from this
channel is `weight × intensity × w / (δ − Δ)`, with `w` the geometric weight of
its ``Δm = q`` component ([`polarisation_weights`](@ref)`[q + 2]` for rank 1,
[`quadrupole_weights`](@ref)`[q + 3]` for rank 2) and `δ` the laser offset from
the same reference `Δ` is measured against. It is positive for a state in the
lower manifold of the pair and negative for one in the upper (level repulsion
away from the driving photon); `4 × |weight| × intensity × w` is the squared
Rabi frequency of the component, and `γ` the total decay rate of the upper
level — both used by the validity warnings at evaluation time.
"""
struct ResonantChannel{S<:StateSpec}
    "The state at the other end (adiabatic labels at the structure's field)."
    partner::S

    "The multipole rank of the level pair (1 = E1, 2 = E2)."
    rank::Int8

    "The ``Δm`` of the component, indexing the geometric weights."
    q::Int8

    "Signed shift weight (see above)."
    weight::typeof(1.0u"rad^2*m^2/(W*µs^2)")

    "Position of the component relative to the pair's zero-field reference."
    Δ::typeof(1.0u"rad/µs")

    "Total decay rate of the upper level of the pair."
    γ::typeof(1.0u"µs^-1")
end

channel_type(basis::StateBasis{L}) where {L} = ResonantChannel{StateSpec{L}}

# Ordered (lower, upper) pairs among the fine-structure `levels` connected by
# an electric-quadrupole transition with a known Einstein A coefficient.
function quadrupole_pairs(species, levels)
    pairs = Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec}[]
    for lo in levels, hi in levels
        species.energies[lo] < species.energies[hi] || continue
        multipole(species, lo => hi) == :E2 || continue
        abs(lo.j - hi.j) <= 2 <= lo.j + hi.j || continue
        isnothing(einstein_a(species, lo => hi)) && continue
        push!(pairs, (lo, hi))
    end
    pairs
end

# Ω²/(4 × intensity) for a unit-amplitude component; see rabi_frequency.
function channel_scale(species, lo, hi, rank)
    a = einstein_a(species, lo => hi)
    ω = transition_frequency(species, lo => hi)
    prefactor = rank == 1 ? 6.0 : 20.0
    uconvert(u"m^2/(W*µs^2)", prefactor * π * u"c"^2 * a / (4 * u"ħ" * ω^3))
end

"""
Validates that the reference pair of a [`RelativeFrequency`](@ref) laser is a
transition the channel machinery can describe — an E1 or E2 one; near-resonant M1
channels are not supported — returning its ordered fine-structure levels and
multipole rank.
"""
function validate_reference(species, laser::RelativeFrequency)
    lo = fine_structure(laser.lower)
    hi = fine_structure(laser.upper)
    m = multipole(species, lo => hi)
    m == :M1 && throw_m1_channel(lo, hi)
    rank = multipole_rank(m)
    if !(abs(lo.j - hi.j) <= rank <= lo.j + hi.j)
        throw(
            ArgumentError(
                "'$(laser.lower)' and '$(laser.upper)' are not connected by an " *
                "$m transition",
            ),
        )
    end
    if isnothing(einstein_a(species, lo => hi))
        if !isnothing(einstein_a(species, hi => lo))
            throw(
                ArgumentError(
                    "'$(laser.lower)' is of higher energy than '$(laser.upper)'; " *
                    "give the reference as lower => upper",
                ),
            )
        end
        throw(
            ArgumentError(
                "No known Einstein A coefficient between '$lo' and '$hi', which " *
                "the near-resonant channel amplitudes require",
            ),
        )
    end
    lo, hi, rank
end

# Near-resonant light-shift channels cover E1 and E2 pairs only.
throw_m1_channel(lo, hi) = throw(
    ArgumentError(
        "'$lo' → '$hi' is a magnetic-dipole transition; near-resonant M1 " *
        "light-shift channels are not supported",
    ),
)

# Maps each E1 reference level to the partner excluded from its background.
function exclusion_partner(species, laser)
    laser isa RelativeFrequency || return Returns(nothing)
    lo = fine_structure(laser.lower)
    hi = fine_structure(laser.upper)
    multipole(species, lo => hi) == :E1 || return Returns(nothing)
    function (level)
        fs = fine_structure(parse_level(level))
        fs == lo ? hi : fs == hi ? lo : nothing
    end
end

# The other fine-structure manifold of the laser's reference pair for the
# given level, or nothing if the level belongs to neither side.
function reference_partner_manifold(laser::RelativeFrequency, level)
    lo = fine_structure(laser.lower)
    hi = fine_structure(laser.upper)
    fs = fine_structure(level)
    fs == lo ? hi : fs == hi ? lo : nothing
end

"""
Builds channels for each E2 pair within the structure's manifolds and for the
[`RelativeFrequency`](@ref) reference pair, if present, from the manifold
solutions at the structure's field (partners outside the basis are solved
into `cache`). Component positions are relative to the named zero-field
interval for the reference pair, or to the centroid interval for other pairs
(which cancels in driven mode).
"""
function resonant_channels(s::Eigenstructure, laser, cache)
    species = s.species
    basis = s.basis
    channels = [channel_type(basis)[] for _ in basis]
    pairs = quadrupole_pairs(species, manifold_levels(s))
    named = nothing
    if laser isa RelativeFrequency
        named = (fine_structure(laser.lower), fine_structure(laser.upper))
        named in pairs || push!(pairs, named)
    end
    for (lo, hi) in pairs
        validate_channel_manifold(species, lo)
        validate_channel_manifold(species, hi)
        rank = multipole_rank(species, lo => hi)
        scale = channel_scale(species, lo, hi, rank)
        γ = uconvert(u"µs^-1", 1 / lifetime(species, hi))
        m_lo = get!(cache, lo) do
            solve_manifold(species, lo, s.B)
        end
        m_hi = get!(cache, hi) do
            solve_manifold(species, hi, s.B)
        end
        rotated = relative_amplitudes(species, m_lo, m_hi, rank)
        # Centroid intervals cancel; subtract the named F levels' Casimir
        # offsets to express component positions relative to the reference.
        ref = 0.0u"rad/µs"
        if (lo, hi) == named
            ref =
                level_offset(species, laser.upper) - level_offset(species, laser.lower)
        end
        for (i, state) in enumerate(basis)
            fs = fine_structure(state.level)
            if fs == lo
                li = stateindex(m_lo.basis, state)
                for (k, partner) in enumerate(m_hi.basis)
                    q = partner.m - state.m
                    abs(q) <= rank || continue
                    amp = rotated[k, li]
                    iszero(amp) && continue
                    Δ = uconvert(u"rad/µs", m_hi.energies[k] - m_lo.energies[li] - ref)
                    push!(
                        channels[i],
                        ResonantChannel(
                            partner,
                            Int8(rank),
                            Int8(Int(q)),
                            uconvert(u"rad^2*m^2/(W*µs^2)", scale * amp^2),
                            Δ,
                            γ,
                        ),
                    )
                end
            elseif fs == hi
                ui = stateindex(m_hi.basis, state)
                for (k, partner) in enumerate(m_lo.basis)
                    q = state.m - partner.m
                    abs(q) <= rank || continue
                    amp = rotated[ui, k]
                    iszero(amp) && continue
                    Δ = uconvert(u"rad/µs", m_hi.energies[ui] - m_lo.energies[k] - ref)
                    push!(
                        channels[i],
                        ResonantChannel(
                            partner,
                            Int8(rank),
                            Int8(Int(q)),
                            uconvert(u"rad^2*m^2/(W*µs^2)", -scale * amp^2),
                            Δ,
                            γ,
                        ),
                    )
                end
            end
        end
    end
    channels
end

# The zero-field offset of a reference level from its fine-structure centroid:
# the Casimir hyperfine shift of an F level, nothing for a fine-structure one.
level_offset(species, level::HyperfineNumberSpec) = hyperfine_shift(species, level)
level_offset(species, level::NoHyperfineNumberSpec) = 0.0u"rad/µs"

"""
A bare laser frequency closer than this to an explicit electric-dipole channel
requires a [`RelativeFrequency`](@ref) reference: its channel detunings must
cancel the uncertain tabulated line centre.
"""
const BACKGROUND_GUARD_WINDOW = u"ħ" * 2π * 100.0u"GHz"

"""
Raises an error if the laser photon energy `ħω` falls within
[`Levels.BACKGROUND_GUARD_WINDOW`](@ref) of an explicit electric-dipole channel
of any of the given levels, except the excluded (reference-pair) partner.
"""
function validate_background_detunings(species, levels, ħω, exclude)
    for level in levels
        fs = fine_structure(parse_level(level))
        data = level_polarisability(species, fs)
        isnothing(data) && continue
        e_level = species.energies[fs]
        for (upper, _) in data.reduced_dipoles
            upper == exclude(level) && continue
            detuning = abs(abs(species.energies[upper] - e_level) - abs(ħω))
            if detuning < BACKGROUND_GUARD_WINDOW
                cyclic = round(u"GHz", detuning / u"ħ" / 2π; digits=1)
                throw(
                    ArgumentError(
                        "The laser is within $cyclic of the '$fs' → '$upper' " *
                        "resonance, where precise component detunings require " *
                        "a reference to the line centre; give the laser as a " *
                        "RelativeFrequency naming that transition instead",
                    ),
                )
            end
        end
    end
end

"""
Light-shift data for the states of an [`Eigenstructure`](@ref), precomputed
for one laser.

Construct with [`LightShiftCoefficients`](@ref)`(s, laser)` and evaluate with
[`light_shift`](@ref) (parked beam, single states) or
[`driven_light_shift`](@ref) (resonantly driven component); see there for the
sign and unit conventions.
"""
struct LightShiftCoefficients{
    S<:Eigenstructure,
    E<:Quantity,
    T<:Quantity,
    R<:Union{Nothing,RelativeFrequency},
}
    "The structure — species, static field and basis — the data was computed for."
    structure::S

    "The photon energy the background polarisabilities were computed for."
    photon_energy::E

    """
    Background (far-detuned electric-dipole) polarisability of each state for
    each ``Δm`` channel, as a `length(basis) × 3` matrix indexed by
    `[stateindex, q + 2]`.

    States of levels without [`LevelPolarisability`](@ref) data carry `NaN`
    entries; evaluating a background shift involving them raises an error. For
    an electric-dipole reference pair, the co-rotating part of that one channel
    is excluded here — it lives in `channels`, Zeeman-resolved, instead.
    """
    polarisabilities::Matrix{T}

    """
    Near-resonant channels of each basis state (cf.
    [`Levels.ResonantChannel`](@ref)), by basis index.
    """
    channels::Vector{Vector{ResonantChannel}}

    "The [`RelativeFrequency`](@ref) the laser was given as, if it was."
    laser::R
end

"""
    LightShiftCoefficients(s::Eigenstructure, laser)

Precomputes light-shift data for every state in `basis(s)`, for a laser given
as a wavelength, angular frequency or [`RelativeFrequency`](@ref). Reuses the
structure's field eigenstates for repeated intensity, polarisation and
detuning evaluations.

The E1 background uses [`LevelPolarisability`](@ref) data; rows without data
are stored as `NaN` and raise an error if evaluated. Near-resonant channels
use Einstein A coefficients and are built for E2 pairs among the structure's
manifolds and for the `RelativeFrequency` pair (E1 or E2). Neither part includes
magnetic-dipole (M1) couplings to the beam's magnetic field; an M1 reference pair
raises an error.

For an E1 reference pair, its co-rotating term is moved from each available
background into the channel sum; the partner must be explicit in that
background's dipole data. A bare frequency within ~100 GHz of an explicit
E1 line is refused; use `RelativeFrequency` to cancel the tabulated line
centre from precise detunings.

[`light_shift`](@ref) evaluates a parked beam and [`driven_light_shift`](@ref)
a resonantly driven component. See [Light shifts and polarisabilities](@ref
reference-polarisability) for the background/channel split and validity limits.
"""
function LightShiftCoefficients(s::Eigenstructure, laser)
    species = s.species
    basis = s.basis
    if laser isa RelativeFrequency
        if !(laser.lower isa eltype(basis.levels))
            throw(
                ArgumentError(
                    "The RelativeFrequency reference levels must be of the same " *
                    "kind as the basis levels",
                ),
            )
        end
        validate_reference(species, laser)
    end
    ħω = photon_energy(species, laser)
    exclude = exclusion_partner(species, laser)
    validate_background_detunings(species, unique(basis.levels), ħω, exclude)
    cache = manifold_cache(s)
    rows = polarisability_rows(s, ħω, exclude, cache)
    channels = resonant_channels(s, laser, cache)
    LightShiftCoefficients(
        s,
        ħω,
        permutedims(reduce(hcat, rows)),
        Vector{Vector{ResonantChannel}}(channels),
        laser isa RelativeFrequency ? laser : nothing,
    )
end

"""
Returns the polarisability of the basis state with the given index for the
[`polarisation_weights`](@ref) `w`.
"""
function state_polarisability(c::LightShiftCoefficients, index::Integer, w)
    if isnan(c.polarisabilities[index, 1])
        throw(
            ArgumentError(
                "No light-shift data known for level '$(c.structure.basis[index].level)'",
            ),
        )
    end
    w[1] * c.polarisabilities[index, 1] +
    w[2] * c.polarisabilities[index, 2] +
    w[3] * c.polarisabilities[index, 3]
end

"""
Converts a polarisability into the corresponding shift at the given intensity.
"""
shift_at_intensity(α, intensity) =
    uconvert(u"rad/µs", -intensity * α / (2 * u"c" * u"ε0" * u"ħ"))

"""
Raises an error unless `parts` is a valid part selector.
"""
function validate_parts(parts)
    if !(parts in (:total, :background, :resonant))
        throw(
            ArgumentError(
                "parts must be :total, :background or :resonant, " *
                "not $(repr(parts))",
            ),
        )
    end
end

# Perturbation-theory validity thresholds of the channel sum: warn when a
# channel Rabi frequency exceeds a tenth of its detuning ((Ω/2d)² > 2.5e-3),
# or the detuning comes within five linewidths of the pole.
const CHANNEL_RABI_RATIO = 0.1
const CHANNEL_LINEWIDTH_RATIO = 5.0

"""
Sums the signed near-resonant shift of the basis state `index` over its
channels to the `other` fine-structure manifold, for a laser at `δ` from the
channels' reference and the given geometric weights (`w_pol` indexed `q + 2`
for rank 1, `w_quad` indexed `q + 3` for rank 2); `exclude` drops the channel
to one partner state (the resonantly driven component of
[`driven_light_shift`](@ref)). Channels the geometry does not drive are
skipped before any division; a driven channel at exactly zero detuning — a
degenerate spectator, for which second-order perturbation theory has nothing
to expand in — raises an error.
"""
function resonant_state_shift(c, index, other, δ, intensity, w_pol, w_quad, exclude)
    shift = 0.0u"rad/µs"
    state = c.structure.basis[index]
    # Channels the geometry does not drive — including numerical residue such
    # as the 1e-17 Δm = -1 weight of a σ⁺ beam — are skipped before any
    # division.
    floor_pol = 1e-12 * sum(w_pol)
    floor_quad = isnothing(w_quad) ? 0.0 : 1e-12 * sum(w_quad)
    for ch in c.channels[index]
        fine_structure(ch.partner.level) == other || continue
        ch.partner == exclude && continue
        w = ch.rank == 1 ? w_pol[Int(ch.q)+2] : w_quad[Int(ch.q)+3]
        w <= (ch.rank == 1 ? floor_pol : floor_quad) && continue
        d = δ - ch.Δ
        if iszero(d)
            throw(
                ArgumentError(
                    "The laser is exactly resonant with the $state → $(ch.partner) " *
                    "component, which the beam geometry drives: its second-order " *
                    "shift is undefined (degenerate spectator, e.g. at zero field)",
                ),
            )
        end
        rabi² = 4 * abs(ch.weight) * intensity * w
        if rabi² > (CHANNEL_RABI_RATIO * d)^2
            @warn "Channel Rabi frequency is comparable to its detuning; " *
                  "second-order perturbation theory may be inaccurate." state = state partner =
                ch.partner rabi = uconvert(u"rad/µs", sqrt(rabi²)) detuning =
                uconvert(u"rad/µs", d)
        end
        if abs(d) < CHANNEL_LINEWIDTH_RATIO * ch.γ
            @warn "Laser within a few linewidths of a channel resonance, where " *
                  "the pole approximation of the shift breaks down." state = state partner =
                ch.partner detuning = uconvert(u"rad/µs", d) linewidth = ch.γ
        end
        shift += uconvert(u"rad/µs", ch.weight * intensity * w / d)
    end
    shift
end

"""
    light_shift(coefficients::LightShiftCoefficients, state::StateSpec, intensity, ε; direction, δ, parts)
    light_shift(s::Eigenstructure, state::StateSpec, laser, intensity, ε; direction, δ, parts)

Returns the angular-frequency shift of `state` under a beam parked at the
coefficients' laser frequency. The state must belong to the structure's
selected basis. A parked transition shift is the upper state's shift minus
the lower state's; [`driven_light_shift`](@ref) instead describes a
resonantly driven component.

`parts` selects `:total` (default), `:background` or `:resonant`. The E1
background depends on polarisation `ε`. For a [`RelativeFrequency`](@ref)
laser, the resonant sum adds channels to the other fine-structure manifold
of the reference pair. `direction` is required for an E2 reference pair and
unused for E1. With a bare frequency only the background is evaluated;
requesting `:resonant` raises an error.

`δ` overrides the reference offset in the channel detunings. The background
stays at the construction frequency; rebuild the coefficients when its
variation matters. The structure form rebuilds them on every call.

The result is a diagonal second-order response, with the degeneracy caveat
of [`polarisability`](@ref). Channel couplings must be small compared with
their detunings. Evaluation warns near strong coupling or linewidth-scale
detunings; a coupled channel at exactly zero detuning raises an error.
"""
function light_shift(
    c::LightShiftCoefficients,
    state::StateSpec,
    intensity,
    ε;
    direction=nothing,
    δ=nothing,
    parts=:total,
)
    validate_parts(parts)
    basis = c.structure.basis
    index = stateindex(basis, state)
    w_pol = polarisation_weights(ε)
    background() = shift_at_intensity(state_polarisability(c, index, w_pol), intensity)
    parts == :background && return background()
    if isnothing(c.laser)
        isnothing(δ) || throw(
            ArgumentError("δ requires a RelativeFrequency laser for its reference"),
        )
        if parts == :resonant
            throw(
                ArgumentError(
                    "The near-resonant shift of a single state requires the " *
                    "laser given as a RelativeFrequency: a bare frequency does " *
                    "not locate it within the Zeeman manifold",
                ),
            )
        end
        # A bare laser has no reference pair for a near-resonant channel sum.
        return background()
    end
    δ_eval = isnothing(δ) ? c.laser.offset : δ
    if !(δ_eval isa Unitful.Frequency)
        throw(ArgumentError("δ must be an angular frequency offset"))
    end
    other = reference_partner_manifold(c.laser, basis[index].level)
    resonant = if isnothing(other)
        0.0u"rad/µs"
    else
        rank = multipole_rank(c.structure.species, c.laser.lower => c.laser.upper)
        w_quad = nothing
        if rank == 2
            isnothing(direction) && throw(
                ArgumentError(
                    "The beam direction is required: the reference pair is an " *
                    "electric-quadrupole transition, whose coupling depends on it",
                ),
            )
            w_quad = quadrupole_weights(ε, direction)
        end
        resonant_state_shift(c, index, other, δ_eval, intensity, w_pol, w_quad, nothing)
    end
    uconvert(u"rad/µs", parts == :resonant ? resonant : background() + resonant)
end

function light_shift(
    s::Eigenstructure,
    state::StateSpec,
    laser,
    intensity,
    ε;
    direction=nothing,
    δ=nothing,
    parts=:total,
)
    c = LightShiftCoefficients(s, laser)
    light_shift(c, state, intensity, ε; direction, δ, parts)
end

function throw_undriveable(c, lower, upper)
    species = c.structure.species
    fs_lo = fine_structure(lower.level)
    fs_hi = fine_structure(upper.level)
    m = multipole(species, fs_lo => fs_hi)
    rank = multipole_rank(m)
    connected =
        abs(fs_lo.j - fs_hi.j) <= rank <= fs_lo.j + fs_hi.j &&
        !isnothing(einstein_a(species, fs_lo => fs_hi))
    if connected
        m == :M1 && throw_m1_channel(fs_lo, fs_hi)
        named =
            !isnothing(c.laser) &&
            (fine_structure(c.laser.lower), fine_structure(c.laser.upper)) ==
            (fs_lo, fs_hi)
        if rank == 1 && !named
            throw(
                ArgumentError(
                    "Driving the electric-dipole '$fs_lo' → '$fs_hi' line " *
                    "requires naming it as the RelativeFrequency reference at " *
                    "construction",
                ),
            )
        end
        Δm = upper.m - lower.m
        Δm = isinteger(Δm) ? Int(Δm) : Δm
        if abs(Δm) > rank
            throw(
                ArgumentError(
                    "The Δm = $Δm component '$lower' => '$upper' of an $m " *
                    "transition cannot be driven (|Δm| ≤ $rank), so there is no " *
                    "resonance whose shift could be observed",
                ),
            )
        end
        throw(
            ArgumentError(
                "The component '$lower' => '$upper' has vanishing amplitude, so " *
                "it cannot be resonantly driven",
            ),
        )
    end
    if !isnothing(einstein_a(species, fs_hi => fs_lo))
        throw(
            ArgumentError(
                "'$lower' is of higher energy than '$upper'; give the transition " *
                "as lower => upper",
            ),
        )
    end
    throw(
        ArgumentError(
            "'$lower' and '$upper' are not connected by a transition with a " *
            "known Einstein A coefficient",
        ),
    )
end

"""
    driven_light_shift(coefficients::LightShiftCoefficients, transition::Pair, intensity, ε; direction, parts)
    driven_light_shift(s::Eigenstructure, lower::StateSpec => upper::StateSpec, intensity, ε; direction, parts)

Returns the angular-frequency displacement of the `lower => upper` resonance
under a laser tuned to that component. A positive result raises the resonance
frequency; both states must be in the structure's selected basis.

`parts` selects `:total` (default), `:background` or `:resonant`. The
background is the difference of the states' E1 shifts. The resonant part sums
spectator channels sharing either state, with detunings relative to the
probed component. The probed channel itself is excluded as the coherent drive.
For E2 this is the ac Stark shift of [Lindvall2025; Sec. III F 2](@cite); its
symmetries are described in [Light shifts and polarisabilities](@ref
reference-polarisability).

An E2 pair needs channels in the coefficients and a beam `direction`. An E1
pair must be the [`RelativeFrequency`](@ref) reference pair; `direction` is
unused. Evaluating the background also requires its construction frequency
to match the probed line within relative tolerance ``10^{-3}``.

The spectator sum uses second-order perturbation theory, with the warnings
of [`light_shift`](@ref). A coupled spectator at the probed frequency raises
an error; a closed cycle without spectators is valid at zero field.
`:background` skips the channels and does not check drivability.
The structure form constructs coefficients with the pair's zero-field
interval as the reference on every call.
"""
function driven_light_shift(
    c::LightShiftCoefficients,
    transition::Pair,
    intensity,
    ε;
    direction=nothing,
    parts=:total,
)
    validate_parts(parts)
    basis = c.structure.basis
    species = c.structure.species
    il = stateindex(basis, transition.first)
    iu = stateindex(basis, transition.second)
    lower, upper = basis[il], basis[iu]
    w_pol = polarisation_weights(ε)
    background = 0.0u"rad/µs"
    if parts != :resonant
        resonance = uconvert(
            u"J",
            u"ħ" * transition_frequency(species, lower.level => upper.level),
        )
        if !isapprox(c.photon_energy, resonance; rtol=1e-3)
            wavelength(ħω) = round(u"nm", 2π * u"ħ" * u"c" / ħω; digits=3)
            throw(
                ArgumentError(
                    "Driven-mode evaluation presumes the laser to be on the " *
                    "probed transition, but the coefficients were computed for " *
                    "$(wavelength(c.photon_energy)) rather than the " *
                    "$(wavelength(resonance)) of '$lower' => '$upper'",
                ),
            )
        end
        background = shift_at_intensity(
            state_polarisability(c, iu, w_pol) - state_polarisability(c, il, w_pol),
            intensity,
        )
        parts == :background && return background
    end
    probed = findfirst(ch -> ch.partner == upper, c.channels[il])
    isnothing(probed) && throw_undriveable(c, lower, upper)
    ch = c.channels[il][probed]
    if ch.weight < zero(ch.weight)
        throw(
            ArgumentError(
                "'$lower' is of higher energy than '$upper'; give the transition " *
                "as lower => upper",
            ),
        )
    end
    w_quad = nothing
    if ch.rank == 2
        isnothing(direction) && throw(
            ArgumentError(
                "The beam direction is required for an electric-quadrupole " *
                "transition, whose coupling depends on it",
            ),
        )
        w_quad = quadrupole_weights(ε, direction)
    end
    fs_lo = fine_structure(lower.level)
    fs_hi = fine_structure(upper.level)
    resonant =
        resonant_state_shift(c, iu, fs_lo, ch.Δ, intensity, w_pol, w_quad, lower) -
        resonant_state_shift(c, il, fs_hi, ch.Δ, intensity, w_pol, w_quad, upper)
    uconvert(u"rad/µs", parts == :resonant ? resonant : background + resonant)
end

function driven_light_shift(
    s::Eigenstructure,
    transition::Pair{<:StateSpec,<:StateSpec},
    intensity,
    ε;
    direction=nothing,
    parts=:total,
)
    lower, upper = transition
    laser = RelativeFrequency(lower.level => upper.level, 0.0u"rad/µs")
    c = LightShiftCoefficients(s, laser)
    driven_light_shift(c, lower => upper, intensity, ε; direction, parts)
end

"""
Returns the level and its light-shift data, or raises an error if there is none.
"""
function polarisability_data(species, level)
    spec = convert(NoHyperfineNumberSpec, level)
    data = level_polarisability(species, spec)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$level'"))
    end
    spec, data
end

"""
    polarisability(species, state::StateSpec, frequency[, ε])
    polarisability(s::Eigenstructure, state::StateSpec, frequency[, ε])

Returns a state's E1 polarisability in SI units (``C m^2/V``). Without `ε`,
returns a vector for pure ``σ^-``, ``π`` and ``σ^+`` fields, indexed by `q + 2`.
With `ε`, contracts it with [`polarisation_weights`](@ref), ignoring the scale
of `ε`. For a linear field, the ``π`` entry applies along ẑ and the mean of
the ``σ`` entries applies perpendicular to it.

`frequency` accepts a wavelength or angular frequency; `0u"rad/s"` gives
the static response. Frequencies within ~100 GHz of an explicit E1 line are
refused. For that regime use [`LightShiftCoefficients`](@ref) with a
[`RelativeFrequency`](@ref), which this accessor does not accept.

The species form uses zero-field states. The structure form accepts any
label in `coupled_basis(s)` and uses its field eigenstate, including ``F``
mixing. Explicit dipole channels resolve source and intermediate hyperfine
and Zeeman energies. A hyperfine intermediate without constants is instead
approximated as degenerate at its centroid, without Zeeman structure.
A static scalar/tensor remainder accounts for unresolved channels.

The result is a diagonal second-order response, including co- and
counter-rotating E1 terms. Within a degenerate subspace the effective
light-shift operator can mix states; its eigenvalues then give the physical
shifts, rather than these diagonal entries. For a sinusoidal field of peak
amplitude ``E``, the cycle-averaged energy shift is ``-α E^2/4``; for a static
field it is ``-α E^2/2``. For light intensity ``I`` the angular shift is
``-I α/(2c ε_0 ħ)``.
"""
function polarisability(s::Eigenstructure, state::StateSpec, frequency)
    ħω = photon_energy(frequency)
    m, index = resolve_state(s, state)
    validate_background_detunings(s.species, [m.level], ħω, Returns(nothing))
    state_polarisabilities(s, state, ħω)
end

polarisability(s::Eigenstructure, state::StateSpec, frequency, ε) =
    sum(polarisation_weights(ε) .* polarisability(s, state, frequency))

function polarisability(species::OneElectronSpecies, state::StateSpec, frequency)
    canonical = validate_state_kind(species, state)
    polarisability(
        Eigenstructure(species, StateBasis([canonical]); B=0.0u"mT"),
        canonical,
        frequency,
    )
end

polarisability(species::OneElectronSpecies, state::StateSpec, frequency, ε) =
    sum(polarisation_weights(ε) .* polarisability(species, state, frequency))

function throw_relative_polarisability()
    throw(
        ArgumentError(
            "polarisability evaluates the plain second-order sum at an absolute " *
            "frequency; a laser placed relative to a line belongs to " *
            "LightShiftCoefficients and light_shift",
        ),
    )
end
polarisability(::OneElectronSpecies, ::StateSpec, ::RelativeFrequency) =
    throw_relative_polarisability()
polarisability(::OneElectronSpecies, ::StateSpec, ::RelativeFrequency, ε) =
    throw_relative_polarisability()
polarisability(::Eigenstructure, ::StateSpec, ::RelativeFrequency) =
    throw_relative_polarisability()
polarisability(::Eigenstructure, ::StateSpec, ::RelativeFrequency, ε) =
    throw_relative_polarisability()

# The zero-field stretched-state response of an F level, from which the
# level's scalar/vector/tensor components are extracted.
function stretched_polarisabilities(species, level::HyperfineNumberSpec, ħω)
    validate_level_kind(species, level)
    state = StateSpec(level, level.f)
    s = Eigenstructure(species, StateBasis([state]); B=0.0u"mT")
    state_polarisabilities(s, state, ħω)
end

"""
    scalar_polarisability(species, level, frequency)

Returns the zero-field scalar polarisability ``α_0(ω)`` at a wavelength or
angular frequency; zero frequency gives the static limit. Accepts a
fine-structure level of either species kind, or an ``F`` level of a hyperfine
species. These level curves have no near-resonance guard.

Together with [`vector_polarisability`](@ref) and
[`tensor_polarisability`](@ref), this decomposes the polarisability of the
sublevel ``m`` for polarisation `ε` as

``α = α_0 + 𝒜 \\frac{m}{2j} α_1 +
\\frac{3 |ε_0|^2 - 1}{2} \\frac{3 m^2 - j(j+1)}{j(2j-1)} α_2``,

where `ε` is normalised, ``𝒜 = |ε_{-1}|^2 - |ε_{+1}|^2`` and
``ε_0 = ε · ẑ``. The vector term is absent for ``j = 0``, the tensor term
for ``j ≤ 1/2``. [`polarisability`](@ref) evaluates the state response directly.

For an ``F`` level, replaces ``j`` by ``F`` and extracts the components from
the zero-field stretched-state responses:
``α_{σ^-} = α_0 - α_1/2 - α_2/2``, ``α_π = α_0 + α_2``,
``α_{σ^+} = α_0 + α_1/2 - α_2/2``. These single-``F`` factors do not apply
to a mixed field eigenstate.
"""
function scalar_polarisability(species, level, frequency)
    ħω = photon_energy(frequency)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        α = stretched_polarisabilities(species, parsed, ħω)
        return (α[1] + α[2] + α[3]) / 3
    end
    spec, data = polarisability_data(species, level)
    e_level = species.energies[spec]
    α = data.static_scalar
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        α += channel_scalar_polarisability(d, spec.j, Δ, ħω)
    end
    α
end

"""
    vector_polarisability(species, level, frequency)

Returns the zero-field vector polarisability ``α_1(ω)`` at a wavelength or
angular frequency; see [`scalar_polarisability`](@ref) for the convention.

Vanishes in the static limit. The model assigns no vector remainder to the
unresolved channels, so only explicit dipoles contribute.
"""
function vector_polarisability(species, level, frequency)
    ħω = photon_energy(frequency)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f > 0 || return 0.0u"C*m^2/V"
        α = stretched_polarisabilities(species, parsed, ħω)
        return α[3] - α[1]
    end
    spec, data = polarisability_data(species, level)
    e_level = species.energies[spec]
    j = spec.j
    α = 0.0u"C*m^2/V"
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        weight = dipole_cg(j, j, 1, upper.j)^2 - dipole_cg(j, j, -1, upper.j)^2
        α += uconvert(u"C*m^2/V", d^2 / (2 * upper.j + 1) * weight * 2ħω / (Δ^2 - ħω^2))
    end
    α
end

"""
    tensor_polarisability(species, level, frequency)

Returns the zero-field tensor polarisability ``α_2(ω)`` at a wavelength or
angular frequency; see [`scalar_polarisability`](@ref) for the convention.

Zero when the queried angular momentum (``J`` or ``F``) is at most ``1/2``.
An ``F ≥ 1`` level of a ``J = 1/2`` parent can have a tensor response from
resolved intermediate hyperfine splittings.
"""
function tensor_polarisability(species, level, frequency)
    ħω = photon_energy(frequency)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f >= 1 || return 0.0u"C*m^2/V"
        α = stretched_polarisabilities(species, parsed, ħω)
        return (2α[2] - α[1] - α[3]) / 3
    end
    spec, data = polarisability_data(species, level)
    j = spec.j
    j > 1//2 || return 0.0u"C*m^2/V"
    e_level = species.energies[spec]
    α = data.static_tensor
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        α += channel_tensor_polarisability(d, j, upper.j, Δ, ħω)
    end
    α
end

Base.show(io::IO, coefficients::LightShiftCoefficients) = print(
    io,
    "LightShiftCoefficients(",
    length(coefficients.structure.basis),
    " states, B = ",
    coefficients.structure.B,
    ")",
)

export LightShiftCoefficients,
    driven_light_shift,
    light_shift,
    polarisability,
    scalar_polarisability,
    tensor_polarisability,
    vector_polarisability
public polarisation_weights, quadrupole_weights
