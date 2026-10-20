using LinearAlgebra: I
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
``2Δ d^2 / [3(2j+1)(Δ^2-ħω^2)]``, for source angular momentum `j`, reduced
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
    quadrupole_weights(ε, n)

Returns the weights ``|Γ_q|^2`` with which the beam intensity drives each of the
``Δm = q`` electric-quadrupole channels, as a length-5 tuple indexed by
`q + 3`.

Uses [`quadrupole_geometry`](@ref) with normalised `ε` and `n`. Both must be
nonzero and transverse, ``ε ⋅ n = 0`` (tolerance ``10^{-3}`` after
normalisation). The weights need not sum to one: only the rank-two part of
``ε ⊗ n`` couples.
"""
function quadrupole_weights(ε, n)
    ε_scale = sqrt(sum(abs2, ε))
    n_scale = sqrt(sum(abs2, n))
    if iszero(ε_scale) || iszero(n_scale)
        throw(ArgumentError("Polarisation and beam direction must be non-zero"))
    end
    ε_unit = ε ./ ε_scale
    n_unit = n ./ n_scale
    if abs(sum(ε_unit .* n_unit)) > 1e-3
        throw(
            ArgumentError(
                "Polarisation must be transverse to the beam direction (ε ⋅ n = 0)",
            ),
        )
    end
    Γ = quadrupole_geometry(ε_unit, n_unit)
    ntuple(i -> abs2(Γ[i]), 5)
end

"""
    state_polarisabilities(species, state::StateSpec, ħω[, B])

Returns the electric-dipole polarisability of `state` for each of the three
driving polarisations at photon energy `ħω`, as a length-3 vector indexed by
`q + 2`. Contracts with [`polarisation_weights`](@ref); cross terms vanish
because every state has definite ``m``.

Uses [`Levels.manifold_polarisabilities`](@ref) at field `B`, or at zero
field when `B = nothing`. The public interface is [`polarisability`](@ref).
"""
function state_polarisabilities(species, state::StateSpec, ħω, B=nothing)
    basis = StateBasis([state])
    validate_basis_kind(species, basis)
    isnothing(B) || validate_shift_field(B)
    if isnothing(level_polarisability(species, state.level))
        throw(ArgumentError("No light-shift data known for level '$(state.level)'"))
    end
    cache = manifold_cache(basis)
    only(polarisability_rows(species, basis, ħω, B, Returns(nothing), cache))
end

"""
Raises an error unless `B` is a real, finite, nonzero scalar flux density.
"""
function validate_shift_field(B)
    if !(B isa Unitful.BField)
        throw(
            ArgumentError(
                "B must be the signed scalar magnetic flux density along the " *
                "quantisation axis, not e.g. the Cartesian field vector " *
                "zeeman_hamiltonian() accepts",
            ),
        )
    end
    if !isreal(B) || !isfinite(B)
        throw(ArgumentError("B must be a real, finite scalar flux density"))
    end
    if iszero(B)
        throw(
            ArgumentError(
                "Field-resolved states and channels require nonzero B; " *
                "omit B for the zero-field states",
            ),
        )
    end
end

"""
Requires ``F`` levels for a hyperfine species and ``J`` levels for a
fine-structure species, so that `m` retains its meaning.
"""
function validate_basis_kind(species, basis::StateBasis)
    if eltype(basis.levels) <: HyperfineNumberSpec
        throw(
            ArgumentError(
                "Level '$(first(basis.levels))' carries hyperfine structure, " *
                "but the species has none",
            ),
        )
    end
end

function validate_basis_kind(species::HyperfineOneElectronSpecies, basis::StateBasis)
    if !(eltype(basis.levels) <: HyperfineNumberSpec)
        throw(
            ArgumentError(
                "States of a hyperfine species must specify hyperfine (F) " *
                "levels, e.g. 'S_1/2 F=4', not '$(first(basis.levels))'",
            ),
        )
    end
end

"""
Requires hyperfine constants for a field-resolved source or resonant partner.
At zero field (`B = nothing`), the coupled states need none.
"""
validate_field_solution(species, fs, B) = nothing

function validate_field_solution(species::HyperfineOneElectronSpecies, fs, B)
    if !isnothing(B) && !haskey(species.hyperfine, fs)
        throw(
            ArgumentError(
                "No hyperfine constants known for level '$fs', whose " *
                "eigenstates at a static field therefore cannot be resolved",
            ),
        )
    end
end

"""
A fine-structure manifold's canonical basis, energies relative to its
centroid (angular units), and eigenvector columns in that basis. An operator
``X`` rotates into the eigenbasis as `states' * X * states`.

Built by [`Levels.manifold_solution`](@ref) for either species kind.
"""
struct ManifoldSolution{L<:LevelSpec}
    "The fine-structure level."
    level::NoHyperfineNumberSpec

    "The canonical basis of the level's states."
    basis::StateBasis{L}

    "The state energies relative to the level centroid, aligned to `basis`."
    energies::Vector{typeof(1.0u"µs^-1")}

    "Eigenvector columns in the coupled basis, aligned to `basis`."
    states::Matrix{Float64}
end

identity_states(num_states) = Matrix{Float64}(I, num_states, num_states)

manifold_cache(::StateBasis{L}) where {L} =
    Dict{NoHyperfineNumberSpec,ManifoldSolution{L}}()

"""
Returns the solution for fine-structure level `fs` at field `B` (`nothing`
for zero field), cached once per coefficient construction.

Fine-structure states retain their basis vectors and acquire Zeeman energies.
Hyperfine states use [`hyperfine_manifold`](@ref) with `B`, or coupled states
with Casimir energies at zero field. An intermediate without hyperfine
constants stays degenerate at its centroid, with no Zeeman shift; callers
use [`Levels.validate_field_solution`](@ref) where this fallback is forbidden.
"""
function manifold_solution(species, fs, B, cache)
    get!(cache, fs) do
        basis = StateBasis([fs])
        energies = [
            isnothing(B) ? 0.0u"µs^-1" :
            uconvert(u"µs^-1", zeeman_shift(species, s, B)) for s in basis
        ]
        ManifoldSolution(fs, basis, energies, identity_states(length(basis)))
    end
end

function manifold_solution(species::HyperfineOneElectronSpecies, fs, B, cache)
    get!(cache, fs) do
        known = haskey(species.hyperfine, fs)
        if known && !isnothing(B)
            m = hyperfine_manifold(species, fs, B)
            return ManifoldSolution(fs, m.basis, m.energies, m.states)
        end
        basis = StateBasis(species, fs)
        energies =
            [known ? hyperfine_shift(species, s.level) : 0.0u"µs^-1" for s in basis]
        ManifoldSolution(fs, basis, energies, identity_states(length(basis)))
    end
end

"""
Returns the `[target, source]` matrix of [`transition_amplitude`](@ref)
values for rank `rank`, rotated into both field eigenbases. Evaluates the
hyperfine reduction once per level pair. Definite ``m`` makes each entry
belong to a single component ``q = m_{target} - m_{source}``.
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

"""
Returns the multiplier of the rank-`rank` Clebsch–Gordan coefficient:
the triangle condition (0 or 1) for fine-structure levels, or the
[`Levels.hyperfine_reduction`](@ref) factor ``β`` for hyperfine levels.
"""
reduction_factor(species, lo::NoHyperfineNumberSpec, hi::NoHyperfineNumberSpec, rank) =
    is_triangle(lo.j, rank, hi.j) ? 1.0 : 0.0

reduction_factor(
    species::HyperfineOneElectronSpecies,
    lo::HyperfineNumberSpec,
    hi::HyperfineNumberSpec,
    rank,
) = hyperfine_reduction(species.nuclear_spin, lo, hi; rank)

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
Returns a `length(source.basis) × 3` matrix of E1 polarisabilities, indexed
by `[state, q + 2]` like [`polarisation_weights`](@ref). Each intermediate
state contributes
``|d_{ts}|^2 [w_q/(E_t-E_s-ħω) + w_{-q}/(E_t-E_s+ħω)]``, with dipole matrix
element ``d_{ts}`` and energies at the common field `B`. Adds the static
scalar remainder and the tensor remainder weighted by
[`Levels.tensor_expectations`](@ref).

`exclude` names an explicit intermediate level whose near-resonant term
moves into [`LightShiftCoefficients`](@ref)'s channel sum. The
counter-rotating term stays here; a lumped remainder cannot be split.
"""
function manifold_polarisabilities(
    species,
    source::ManifoldSolution,
    ħω,
    B,
    exclude,
    cache,
)
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
        target = manifold_solution(species, upper, B, cache)
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
Returns background polarisabilities in `basis` order as length-3 rows,
evaluating each manifold once. Levels without data have `NaN` rows.
`exclude(level)` supplies each level's excluded intermediate partner.
"""
function polarisability_rows(species, basis::StateBasis, ħω, B, exclude, cache)
    blocks = Dict{NoHyperfineNumberSpec,Matrix{typeof(1.0u"C*m^2/V")}}()
    for level in basis.levels
        fs = fine_structure(level)
        (haskey(blocks, fs) || isnothing(level_polarisability(species, fs))) && continue
        validate_field_solution(species, fs, B)
        source = manifold_solution(species, fs, B, cache)
        blocks[fs] =
            manifold_polarisabilities(species, source, ħω, B, exclude(fs), cache)
    end
    [
        begin
            fs = fine_structure(state.level)
            haskey(blocks, fs) ? blocks[fs][stateindex(cache[fs].basis, state), :] :
            fill(NaN * u"C*m^2/V", 3)
        end for state in basis
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
    "The state at the other end (adiabatic labels at the construction field)."
    partner::S

    "The multipole rank of the level pair (1 = E1, 2 = E2)."
    rank::Int8

    "The ``Δm`` of the component, indexing the geometric weights."
    q::Int8

    "Signed shift weight (see above)."
    weight::typeof(1.0u"m^2/(W*µs^2)")

    "Position of the component relative to the pair's zero-field reference."
    Δ::typeof(1.0u"µs^-1")

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
        multipole_rank(lo, hi) == 2 || continue
        abs(lo.j - hi.j) <= 2 <= lo.j + hi.j || continue
        isnothing(einstein_a(species, lo, hi)) && continue
        push!(pairs, (lo, hi))
    end
    pairs
end

# Ω²/(4 × intensity) for a unit-amplitude component; see rabi_frequency.
function channel_scale(species, lo, hi, rank)
    a = einstein_a(species, lo, hi)
    ω = transition_frequency(species, lo, hi)
    prefactor = rank == 1 ? 6.0 : 20.0
    uconvert(u"m^2/(W*µs^2)", prefactor * π * u"c"^2 * a / (4 * u"ħ" * ω^3))
end

"""
Validates that the reference pair of a [`RelativeFrequency`](@ref) laser is a
transition the channel machinery can describe, returning its ordered
fine-structure levels and multipole rank.
"""
function validate_reference(species, laser::RelativeFrequency)
    lo = fine_structure(laser.lower)
    hi = fine_structure(laser.upper)
    rank = multipole_rank(lo, hi)
    if !(abs(lo.j - hi.j) <= rank <= lo.j + hi.j)
        throw(
            ArgumentError(
                "'$(laser.lower)' and '$(laser.upper)' are not connected by an " *
                "E$rank transition",
            ),
        )
    end
    if isnothing(einstein_a(species, lo, hi))
        if !isnothing(einstein_a(species, hi, lo))
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

# Maps each E1 reference level to the partner excluded from its background.
function exclusion_partner(laser)
    laser isa RelativeFrequency || return Returns(nothing)
    lo = fine_structure(laser.lower)
    hi = fine_structure(laser.upper)
    multipole_rank(lo, hi) == 1 || return Returns(nothing)
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
Builds channels for each E2 pair within `basis` and for the
[`RelativeFrequency`](@ref) reference pair, if present. Uses the cached
manifold solutions at `B` for amplitudes and energies. Component positions
are relative to the named zero-field interval for the reference pair, or
to the centroid interval for other pairs (which cancels in driven mode).
"""
function resonant_channels(species, basis::StateBasis, laser, B, cache)
    channels = [channel_type(basis)[] for _ in basis]
    pairs = quadrupole_pairs(species, unique!(fine_structure.(basis.levels)))
    named = nothing
    if laser isa RelativeFrequency
        named = (fine_structure(laser.lower), fine_structure(laser.upper))
        named in pairs || push!(pairs, named)
    end
    for (lo, hi) in pairs
        validate_field_solution(species, lo, B)
        validate_field_solution(species, hi, B)
        rank = multipole_rank(lo, hi)
        scale = channel_scale(species, lo, hi, rank)
        γ = uconvert(u"µs^-1", 1 / lifetime(species, hi))
        m_lo = manifold_solution(species, lo, B, cache)
        m_hi = manifold_solution(species, hi, B, cache)
        rotated = relative_amplitudes(species, m_lo, m_hi, rank)
        # Centroid intervals cancel; subtract the named F levels' Casimir
        # offsets to express component positions relative to the reference.
        ref = 0.0u"µs^-1"
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
                    Δ = uconvert(u"µs^-1", m_hi.energies[k] - m_lo.energies[li] - ref)
                    push!(
                        channels[i],
                        ResonantChannel(
                            partner,
                            Int8(rank),
                            Int8(Int(q)),
                            scale * amp^2,
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
                    Δ = uconvert(u"µs^-1", m_hi.energies[ui] - m_lo.energies[k] - ref)
                    push!(
                        channels[i],
                        ResonantChannel(
                            partner,
                            Int8(rank),
                            Int8(Int(q)),
                            -scale * amp^2,
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
level_offset(species, level::NoHyperfineNumberSpec) = 0.0u"µs^-1"

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
Light-shift data for a set of states, precomputed for one laser.

Construct with [`LightShiftCoefficients`](@ref)`(species, basis, laser[; B])`
and evaluate with [`light_shift`](@ref) (parked beam, single states) or
[`driven_light_shift`](@ref) (resonantly driven component); see there for the
sign and unit conventions.
"""
struct LightShiftCoefficients{
    S,
    L<:LevelSpec,
    E<:Quantity,
    T<:Quantity,
    R<:Union{Nothing,RelativeFrequency},
    F<:Union{Nothing,Quantity},
}
    "The species the data was computed for."
    species::S

    "The states the coefficients refer to, fixing the row order."
    basis::StateBasis{L}

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
    [`Levels.ResonantChannel`](@ref)), by basis index; empty unless the static
    flux density `B` was given at construction.
    """
    channels::Vector{Vector{ResonantChannel{StateSpec{L}}}}

    "The [`RelativeFrequency`](@ref) the laser was given as, if it was."
    laser::R

    "The static flux density the background and channels were computed at, if any."
    field::F
end

"""
    LightShiftCoefficients(species, basis::StateBasis, laser[; B])
    LightShiftCoefficients(species, levels_or_states::AbstractVector, laser[; B])

Precomputes the ac Stark shift data of every state in the given basis, for a
laser given as a wavelength, an angular frequency, or a
[`RelativeFrequency`](@ref).

Caches the atomic structure for repeated intensity, polarisation and detuning
evaluations, for example in laser-parameter fits.

The E1 **background** requires [`LevelPolarisability`](@ref) data. Levels
without it are admitted, but evaluating their background raises. With `B`,
builds **near-resonant channels** for each E2 pair within the basis and the
`RelativeFrequency` reference pair (E1 or E2). Channels need only Einstein A
coefficients. Evaluates with [`light_shift`](@ref) for a parked beam, or
[`driven_light_shift`](@ref) for a resonantly driven component.

A `RelativeFrequency` laser requires `B`. For an E1 reference pair, moves
the near-resonant term from the background into the channels; that pair
must be explicit in the polarisability data. Bare frequencies within
~100 GHz of an explicit E1 line are refused because precise detunings
require cancellation of the uncertain tabulated line centre.

With `B`, the background describes field eigenstates, including hyperfine
``F`` mixing and Zeeman energies; without it, zero-field states. See
[`polarisability`](@ref) for the conventions and intermediate-state fallback.
"""
LightShiftCoefficients(species, basis::StateBasis, laser; B=nothing) =
    build_coefficients(species, basis, laser, B)

LightShiftCoefficients(species, levels_or_states::AbstractVector, laser; kwargs...) =
    LightShiftCoefficients(species, StateBasis(levels_or_states), laser; kwargs...)

function build_coefficients(species, basis, laser, B)
    validate_basis_kind(species, basis)
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
        if isnothing(B)
            throw(
                ArgumentError(
                    "A RelativeFrequency laser requires the static flux density " *
                    "B, as the near-resonant channels are field-resolved",
                ),
            )
        end
    end
    isnothing(B) || validate_shift_field(B)
    ħω = photon_energy(species, laser)
    exclude = exclusion_partner(laser)
    validate_background_detunings(species, unique(basis.levels), ħω, exclude)
    cache = manifold_cache(basis)
    rows = polarisability_rows(species, basis, ħω, B, exclude, cache)
    channels =
        isnothing(B) ? [channel_type(basis)[] for _ in basis] :
        resonant_channels(species, basis, laser, B, cache)
    LightShiftCoefficients(
        species,
        basis,
        ħω,
        permutedims(reduce(hcat, rows)),
        channels,
        laser isa RelativeFrequency ? laser : nothing,
        isnothing(B) ? nothing : uconvert(u"mT", B),
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
                "No light-shift data known for level '$(c.basis[index].level)'",
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
    uconvert(u"µs^-1", -intensity * α / (2 * u"c" * u"ε0" * u"ħ"))

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
[`driven_light_shift`](@ref)).
"""
function resonant_state_shift(c, index, other, δ, intensity, w_pol, w_quad, exclude)
    shift = 0.0u"µs^-1"
    state = c.basis[index]
    for ch in c.channels[index]
        fine_structure(ch.partner.level) == other || continue
        ch.partner == exclude && continue
        w = ch.rank == 1 ? w_pol[Int(ch.q)+2] : w_quad[Int(ch.q)+3]
        iszero(w) && continue
        d = δ - ch.Δ
        rabi² = 4 * abs(ch.weight) * intensity * w
        if rabi² > (CHANNEL_RABI_RATIO * d)^2
            @warn "Channel Rabi frequency is comparable to its detuning; " *
                  "second-order perturbation theory may be inaccurate." state = state partner =
                ch.partner rabi = uconvert(u"µs^-1", sqrt(rabi²)) detuning =
                uconvert(u"µs^-1", d)
        end
        if abs(d) < CHANNEL_LINEWIDTH_RATIO * ch.γ
            @warn "Laser within a few linewidths of a channel resonance, where " *
                  "the pole approximation of the shift breaks down." state = state partner =
                ch.partner detuning = uconvert(u"µs^-1", d) linewidth = ch.γ
        end
        shift += ch.weight * intensity * w / d
    end
    shift
end

"""
    light_shift(coefficients::LightShiftCoefficients, state::StateSpec, intensity, ε; n, δ, parts)
    light_shift(species, state::StateSpec, laser, intensity, ε; n, B, δ, parts)

Returns the ac Stark (light) shift of the given state (in angular units) for a
beam of the given intensity and polarisation `ε`, parked at the laser frequency
the coefficients were built for.

The shift is the far-detuned electric-dipole **background** — which depends on
the polarisation alone, never on the beam direction — plus, for a
[`RelativeFrequency`](@ref) laser (necessarily built with `B`), the
**near-resonant** sum over the Zeeman components connecting the state to the
other manifold of the reference pair, at their exact at-field detunings from
the stated laser offset. For an electric-quadrupole reference pair the beam
direction `n` is required, as the quadrupole coupling depends on it; for an
electric-dipole one it is not used. `parts` selects `:total` (the default),
`:background` or `:resonant`.

`δ` overrides the [`RelativeFrequency`](@ref) offset in the channel sum.
The smooth background remains at the construction frequency; rebuild the
coefficients if its variation matters over a sweep.

The change a parked beam makes to the splitting of a `lower => upper` pair — a
Ramsey-type measurement, with no component resonantly driven — is the
difference of the two states' shifts; the displacement of the resonance
observed when *driving* a component is [`driven_light_shift`](@ref) instead.

The near-resonant sum is second-order perturbation theory: every channel Rabi
frequency must stay small against its detuning, and detunings within a few
linewidths of a pole are outside the model; both conditions warn when
violated. The one-shot species form is a convenience wrapper that rebuilds the
coefficients on every call.

`B` also selects field eigenstates for the background, including for a bare
frequency or wavelength. Without `B`, it describes zero-field states; see
[`polarisability`](@ref).
"""
function light_shift(
    c::LightShiftCoefficients,
    state::StateSpec,
    intensity,
    ε;
    n=nothing,
    δ=nothing,
    parts=:total,
)
    validate_parts(parts)
    index = stateindex(c.basis, state)
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
                    "laser given as a RelativeFrequency (with B): a bare " *
                    "frequency does not locate it within the Zeeman manifold",
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
    other = reference_partner_manifold(c.laser, c.basis[index].level)
    resonant = if isnothing(other)
        0.0u"µs^-1"
    else
        rank = multipole_rank(fine_structure(c.laser.lower), fine_structure(c.laser.upper))
        w_quad = nothing
        if rank == 2
            isnothing(n) && throw(
                ArgumentError(
                    "The beam direction n is required: the reference pair is an " *
                    "electric-quadrupole transition, whose coupling depends on it",
                ),
            )
            w_quad = quadrupole_weights(ε, n)
        end
        resonant_state_shift(c, index, other, δ_eval, intensity, w_pol, w_quad, nothing)
    end
    parts == :resonant ? resonant : background() + resonant
end

function light_shift(
    species,
    state::StateSpec,
    laser,
    intensity,
    ε;
    n=nothing,
    B=nothing,
    δ=nothing,
    parts=:total,
)
    c = LightShiftCoefficients(species, [state.level], laser; B)
    light_shift(c, state, intensity, ε; n, δ, parts)
end

# The pair forms of light_shift were replaced by the parked/driven split;
# fail with guidance rather than a MethodError.
light_shift(c::LightShiftCoefficients, transition::Pair, intensity, ε; kwargs...) =
    throw(
        ArgumentError(
            "The light shift of a pair is driven_light_shift (laser locked to " *
            "the driven component) or, for a parked beam, the difference of the " *
            "two states' light_shift results",
        ),
    )
light_shift(
    species,
    lower::StateSpec,
    upper::StateSpec,
    laser,
    intensity,
    ε;
    kwargs...,
) = throw(
    ArgumentError(
        "The light shift of a pair is driven_light_shift (laser locked to " *
        "the driven component; no laser argument — the probed transition fixes " *
        "it) or, for a parked beam, the difference of the two states' " *
        "light_shift results",
    ),
)

# Explains why no channel connects the probed pair, mirroring the checks the
# construction applies.
function throw_undriveable(c, lower, upper)
    fs_lo = fine_structure(lower.level)
    fs_hi = fine_structure(upper.level)
    rank = multipole_rank(fs_lo, fs_hi)
    connected =
        abs(fs_lo.j - fs_hi.j) <= rank <= fs_lo.j + fs_hi.j &&
        !isnothing(einstein_a(c.species, fs_lo, fs_hi))
    if connected
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
                    "The Δm = $Δm component '$lower' => '$upper' of an E$rank " *
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
    if !isnothing(einstein_a(c.species, fs_hi, fs_lo))
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
    driven_light_shift(coefficients::LightShiftCoefficients, transition::Pair, intensity, ε; n, parts)
    driven_light_shift(species, lower::StateSpec, upper::StateSpec, intensity, ε; n, B, parts)

Returns the angular-frequency displacement of the `lower => upper` resonance
under a laser tuned to that component. A positive result raises the resonance
frequency.

The **background** part is the difference of the two states' far-detuned
electric-dipole shifts. The **resonant** part sums, for both states, the
couplings to every *other* Zeeman component sharing one of them — the driven
channel itself is the coherent drive, not a shift — with the laser pinned to
the probed component, so the laser frequency drops out and only exact at-field
splittings enter the detunings. `parts` selects `:total` (the default),
`:background` or `:resonant`.

For an E2 line, `parts = :resonant` is the "E2 ac Stark shift" of
`[Lindvall2025]`, Sec. III F 2. In a fine-structure species, linear
polarisation makes this part odd under reversal of both projections, so it
cancels in a ``±m`` pair average. Ellipticity breaks this cancellation. For
hyperfine states, spectators in other ``F`` levels generally prevent it;
with linear polarisation the field-reversal identity is instead
``shift(-m, B) = shift(+m, -B)``, with both projections reversed.

The probed pair may be any electric-quadrupole transition with channels in the
coefficients (which requires construction with `B`); an electric-dipole pair
must additionally be the [`RelativeFrequency`](@ref) reference pair. The beam
direction `n` is required for electric-quadrupole pairs and unused for
electric-dipole ones. For the background part, the coefficients' laser
frequency must be consistent with the probed transition (rtol ``10^{-3}``) —
in driven mode the laser *is* on that line. Like the parked
[`light_shift`](@ref), the resonant sum is second-order perturbation theory
and warns when a spectator Rabi frequency approaches its detuning (which
breaks down first for near-degenerate components at very low fields).

# References

- `[Lindvall2025]`: T. Lindvall, A. E. Wallin, K. J. Hanhijärvi, and T. Fordell,
  "⁸⁸Sr⁺ Optical Clock with 7.9 × 10⁻¹⁹ Systematic Uncertainty and Measurement of
  Its Absolute Frequency", Phys. Rev. Applied **24**, 044082 (2025),
  [doi:10.1103/cztf-bfvp](https://doi.org/10.1103/cztf-bfvp).
"""
function driven_light_shift(
    c::LightShiftCoefficients,
    transition::Pair,
    intensity,
    ε;
    n=nothing,
    parts=:total,
)
    validate_parts(parts)
    il = stateindex(c.basis, transition.first)
    iu = stateindex(c.basis, transition.second)
    lower, upper = c.basis[il], c.basis[iu]
    w_pol = polarisation_weights(ε)
    background = 0.0u"µs^-1"
    if parts != :resonant
        resonance = uconvert(
            u"J",
            u"ħ" * transition_frequency(c.species, lower.level, upper.level),
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
    if isnothing(c.field)
        throw(
            ArgumentError(
                "The near-resonant part requires the channels: construct the " *
                "LightShiftCoefficients with the static flux density B",
            ),
        )
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
        isnothing(n) && throw(
            ArgumentError(
                "The beam direction n is required for an electric-quadrupole " *
                "transition, whose coupling depends on it",
            ),
        )
        w_quad = quadrupole_weights(ε, n)
    end
    fs_lo = fine_structure(lower.level)
    fs_hi = fine_structure(upper.level)
    resonant =
        resonant_state_shift(c, iu, fs_lo, ch.Δ, intensity, w_pol, w_quad, lower) -
        resonant_state_shift(c, il, fs_hi, ch.Δ, intensity, w_pol, w_quad, upper)
    parts == :resonant ? resonant : background + resonant
end

function driven_light_shift(
    species,
    lower::StateSpec,
    upper::StateSpec,
    intensity,
    ε;
    n=nothing,
    B=nothing,
    parts=:total,
)
    levels = [lower.level, upper.level]
    if isnothing(B)
        if parts != :background
            throw(
                ArgumentError(
                    "The near-resonant part requires the static flux density B",
                ),
            )
        end
        fs_lo = fine_structure(parse_level(lower.level))
        fs_hi = fine_structure(parse_level(upper.level))
        if multipole_rank(fs_lo, fs_hi) == 1
            throw(
                ArgumentError(
                    "The background of a driven electric-dipole line requires B: " *
                    "separating out its resonant channel is field-resolved",
                ),
            )
        end
        # For an electric-quadrupole pair the background needs no channels, so
        # a bare construction at the line's own frequency suffices.
        c = LightShiftCoefficients(
            species,
            levels,
            transition_frequency(species, lower.level, upper.level),
        )
        return driven_light_shift(c, lower => upper, intensity, ε; n, parts)
    end
    laser = RelativeFrequency(lower.level => upper.level, 0.0u"µs^-1")
    c = LightShiftCoefficients(species, levels, laser; B)
    driven_light_shift(c, lower => upper, intensity, ε; n, parts)
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
    polarisability(species, state::StateSpec, frequency[, ε]; B=nothing)

Returns one state's E1 polarisability in SI units (``C m^2/V``). Without `ε`,
returns a length-3 vector for pure ``σ^-``, ``π`` and ``σ^+`` driving fields,
indexed by `q + 2`. With `ε`, contracts that vector with
[`polarisation_weights`](@ref); the scale of `ε` is ignored. A real `ε`
specifies a linearly oscillating field: the ``π`` entry applies along ẑ,
and the mean of the ``σ`` entries applies perpendicular to it.

`frequency` is the wavelength or angular frequency of the driving field — a
laser or a trap rf drive; `0u"s^-1"` gives the static response. Frequencies
within ~100 GHz of an explicit E1 line are refused because precise detunings
require cancellation of the uncertain line centre. That regime is handled
by [`LightShiftCoefficients`](@ref) and [`light_shift`](@ref); this accessor
rejects [`RelativeFrequency`](@ref).

Without `B`, the state is a zero-field ``|J, m_J⟩`` or ``|F, m_F⟩`` state.
With `B` (signed, finite, nonzero scalar flux density along ẑ), a hyperfine
state denotes the adiabatically labelled eigenstate, including ``F`` mixing.
Fine-structure eigenvectors are unchanged. The sum resolves known hyperfine
splittings and, with `B`, Zeeman energies in both manifolds. An intermediate
without hyperfine constants stays degenerate at its centroid, with no Zeeman
shift. The static remainder is evaluated as an electronic operator in the
source state.

Returns the diagonal second-order E1 response, including rotating and
counter-rotating terms and a static remainder; excludes E2 channels. For a
sinusoidal field of peak amplitude `E`, the cycle-averaged energy shift is
``-α E^2/4``; for a static field it is ``-α E^2/2``. For light intensity `I`,
the angular shift is ``-I α/(2c ε_0 ħ)``.
"""
function polarisability(species, state::StateSpec, frequency; B=nothing)
    ħω = photon_energy(frequency)
    validate_background_detunings(species, [state.level], ħω, Returns(nothing))
    state_polarisabilities(species, state, ħω, B)
end

polarisability(species, state::StateSpec, frequency, ε; B=nothing) =
    sum(polarisation_weights(ε) .* polarisability(species, state, frequency; B))

polarisability(
    species,
    state::StateSpec,
    laser::RelativeFrequency,
    args...;
    kwargs...,
) = throw(
    ArgumentError(
        "polarisability evaluates the plain second-order sum at an absolute " *
        "frequency; a laser placed relative to a line belongs to " *
        "LightShiftCoefficients and light_shift",
    ),
)

"""
    scalar_polarisability(species, level, laser)

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
function scalar_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        α = state_polarisabilities(species, StateSpec(parsed, parsed.f), ħω)
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
    vector_polarisability(species, level, laser)

Returns the zero-field vector polarisability ``α_1(ω)`` at a wavelength or
angular frequency; see [`scalar_polarisability`](@ref) for the convention.

Vanishes in the static limit. The model assigns no vector remainder to the
unresolved channels, so only explicit dipoles contribute.
"""
function vector_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f > 0 || return 0.0u"C*m^2/V"
        α = state_polarisabilities(species, StateSpec(parsed, parsed.f), ħω)
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
    tensor_polarisability(species, level, laser)

Returns the zero-field tensor polarisability ``α_2(ω)`` at a wavelength or
angular frequency; see [`scalar_polarisability`](@ref) for the convention.

Zero for ``j ≤ 1/2``, which cannot support a rank-two tensor response.
"""
function tensor_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f >= 1 || return 0.0u"C*m^2/V"
        α = state_polarisabilities(species, StateSpec(parsed, parsed.f), ħω)
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

export LightShiftCoefficients,
    driven_light_shift,
    light_shift,
    polarisability,
    scalar_polarisability,
    tensor_polarisability,
    vector_polarisability
public polarisation_weights, quadrupole_weights
