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

"""
    polarisation_weights(ε)

Returns the fractions ``w_q = |ε_{-q}|^2`` of the beam intensity driving each of
the ``Δm = q`` electric-dipole channels, as a length-3 tuple indexed by
`q + 2`.

These are the squared magnitudes of the [`dipole_geometry`](@ref) amplitudes,
normalised to sum to one; only the direction of `ε` matters, as the overall
scale is set by the intensity.
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

These are the squared magnitudes of the [`quadrupole_geometry`](@ref) amplitudes
for a unit polarisation and beam direction, so only the directions of `ε` and
`n` matter; the overall scale is set by the intensity. Unlike the
electric-dipole [`polarisation_weights`](@ref) they do not sum to one, as only
part of ``ε ⊗ n`` is of rank two.

The polarisation must be transverse to the beam direction, ``ε ⋅ n = 0``, as it
is for any propagating field; pairs that do not belong to one physical beam are
rejected.
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
Returns the polarisability of `state` for each of the three ``Δm`` channels
separately, as a length-3 vector indexed by `q + 2`.

For hyperfine states of a [`HyperfineOneElectronSpecies`](@ref), the explicit
channels are summed over the hyperfine levels of each intermediate
fine-structure level, with the ``F``-basis angular factors
(``\\mathrm{CG} × β``, cf. [`Levels.hyperfine_reduction`](@ref)) and the
intermediate-state hyperfine splittings resolved in the detunings (levels
without known hyperfine constants keep degenerate ``F`` levels at their
centroid). The hyperfine-resolved detunings are what gives e.g. the ⁴³Ca⁺
S``_{1/2}`` ``F`` levels their small tensor polarisability (``∝ A_P/Δ``),
which the fine-structure ``J = 1/2`` level lacks; with the splittings zeroed,
the result collapses exactly to the standard 6-j re-projection of the
fine-structure polarisability. The tensor part of the lumped static remainder
is re-projected with the expectation value of the ``J``-basis tensor operator
in the coupled state (diagonal-in-``F`` approximation; the ``F``-off-diagonal
elements are second order in the hyperfine mixing).

The `exclude` keyword names one intermediate fine-structure level whose
*co-rotating* term is left out of the sum — used when that channel is
near-resonant and handled Zeeman-resolved by the channel machinery of
[`LightShiftCoefficients`](@ref) instead; its smooth counter-rotating term
stays in the background. The excluded channel must be explicit in the
[`LevelPolarisability`](@ref) data, not lumped into the static remainder,
which cannot be split.
"""
function state_polarisabilities(species, state::StateSpec, ħω; exclude=nothing)
    level = convert(NoHyperfineNumberSpec, state.level)
    data = level_polarisability(species, level)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$(state.level)'"))
    end
    if !isnothing(exclude) && !haskey(data.reduced_dipoles, exclude)
        throw(
            ArgumentError(
                "The '$level' → '$exclude' channel is lumped into the static " *
                "remainder of the polarisability data, so its near-resonant part " *
                "cannot be separated out",
            ),
        )
    end
    j, m = level.j, state.m
    if abs(m) > j || !isinteger(j - m)
        throw(ArgumentError("Invalid projection m = $m for level with j = $j"))
    end
    e_level = species.energies[level]

    α = fill(0.0u"C*m^2/V", 3)
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        r = d^2 / (2 * upper.j + 1)
        # For a channel to a level below, the near-resonant term is the one at
        # Δ + ħω (the labels swap), so the exclusion picks by the sign of Δ.
        drop_rotating = upper == exclude && Δ > zero(Δ)
        drop_counter = upper == exclude && !drop_rotating
        for q in -1:1
            rotating = drop_rotating ? 0.0 : dipole_cg(j, m, q, upper.j)^2
            counter = drop_counter ? 0.0 : dipole_cg(j, m, -q, upper.j)^2
            α[q+2] +=
                uconvert(u"C*m^2/V", r * (rotating / (Δ - ħω) + counter / (Δ + ħω)))
        end
    end

    # The tensor part of the lumped remainder scales with the same angular factors
    # as an explicit channel would; the scalar part is channel-independent.
    tensor = j > 1//2 ? (3m^2 - j * (j + 1)) / (j * (2j - 1)) : 0//1
    for q in -1:1
        α[q+2] +=
            data.static_scalar + data.static_tensor * ((3 * (q == 0) - 1) / 2) * tensor
    end
    α
end

function state_polarisabilities(
    species::HyperfineOneElectronSpecies,
    state::StateSpec{HyperfineNumberSpec},
    ħω;
    exclude=nothing,
)
    spec = validate_hyperfine(species, state.level)
    fs = fine_structure(spec)
    data = level_polarisability(species, fs)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$(state.level)'"))
    end
    if !isnothing(exclude) && !haskey(data.reduced_dipoles, exclude)
        throw(
            ArgumentError(
                "The '$fs' → '$exclude' channel is lumped into the static " *
                "remainder of the polarisability data, so its near-resonant part " *
                "cannot be separated out",
            ),
        )
    end
    f, m = spec.f, state.m
    if abs(m) > f || !isinteger(f - m)
        throw(ArgumentError("Invalid projection m = $m for level with F = $f"))
    end
    e_level = level_energy(species, spec)

    # Zero-field hyperfine shift where constants are known; intermediate levels
    # without them keep their F levels degenerate at the centroid.
    shift_of(level) =
        haskey(species.hyperfine, fine_structure(level)) ?
        u"ħ" * hyperfine_shift(species, level) : zero(1.0u"J")

    α = fill(0.0u"C*m^2/V", 3)
    for (upper, d) in data.reduced_dipoles
        r = d^2 / (2 * upper.j + 1)
        for upper_f in hyperfine_levels(species, upper)
            β = hyperfine_reduction(species.nuclear_spin, spec, upper_f; rank=1)
            iszero(β) && continue
            Δ = species.energies[upper] + shift_of(upper_f) - e_level
            drop_rotating = upper == exclude && Δ > zero(Δ)
            drop_counter = upper == exclude && !drop_rotating
            for q in -1:1
                rotating = drop_rotating ? 0.0 : dipole_cg(f, m, q, upper_f.f)^2
                counter = drop_counter ? 0.0 : dipole_cg(f, m, -q, upper_f.f)^2
                α[q+2] += uconvert(
                    u"C*m^2/V",
                    r * β^2 * (rotating / (Δ - ħω) + counter / (Δ + ħω)),
                )
            end
        end
    end

    # The J-basis tensor operator of the lumped remainder, evaluated in the
    # coupled state (diagonal-in-F approximation): its angular factor is the
    # CG-weighted average of (3m_J² − J(J+1))/(J(2J−1)) over the state's
    # |m_I, m_J⟩ decomposition. Zero for J = 1/2, where the operator vanishes
    # identically — the genuinely-new hyperfine tensor term of such levels
    # comes from the resolved detunings above instead.
    j = fs.j
    i_nuc = species.nuclear_spin
    tensor = if j > 1//2
        sum(
            Float64(clebschgordan(i_nuc, m - m_j, j, m_j, f, m))^2 *
            (3m_j^2 - j * (j + 1)) / (j * (2j - 1)) for
            m_j in (-j):j if abs(m - m_j) <= i_nuc;
            init=0.0,
        )
    else
        0.0
    end
    for q in -1:1
        α[q+2] +=
            data.static_scalar + data.static_tensor * ((3 * (q == 0) - 1) / 2) * tensor
    end
    α
end

"""
Raises an error unless `B` is a sensible signed scalar flux density for the
near-resonant channel construction.
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
    if iszero(B)
        throw(
            ArgumentError(
                "The near-resonant channels are undefined at zero magnetic " *
                "field, where the Zeeman components are degenerate",
            ),
        )
    end
end

# Relative rank-`rank` multipole amplitudes between the adiabatically-labelled
# eigenstates of two solved manifolds: the coupled-basis CG × β amplitudes
# conjugated with the eigenvector matrices ([upper eigenstate, lower
# eigenstate], real).
function rotated_multipole_amplitudes(
    species::HyperfineOneElectronSpecies,
    m_lower::HyperfineManifold,
    m_upper::HyperfineManifold,
    rank::Integer,
)
    c = zeros(length(m_upper.basis), length(m_lower.basis))
    for (i, ls) in enumerate(m_lower.basis), (k, us) in enumerate(m_upper.basis)
        q = us.m - ls.m
        abs(q) <= rank || continue
        β = hyperfine_reduction(species.nuclear_spin, ls.level, us.level; rank)
        iszero(β) && continue
        c[k, i] =
            Float64(clebschgordan(ls.level.f, ls.m, rank, q, us.level.f, us.m)) * β
    end
    m_upper.states' * c * m_lower.states
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

# Common scale of the channel weights of one (lower, upper) level pair: the
# James-formula Ω²-per-intensity prefactor over four (cf. rabi_frequency), so
# that weight × intensity × w is the Ω²/4 of a unit-amplitude component.
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

# The background-exclusion partner function of the given laser: maps a basis
# level to the other member of an electric-dipole reference pair (whose
# co-rotating term moves into the channel sum), or to nothing.
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
Builds the per-state near-resonant channel lists over `basis` at the static
flux density `B`: for every state, one entry per Zeeman component connecting it
to the other manifold of each relevant level pair — the electric-quadrupole
pairs within the basis, plus the reference pair of a
[`RelativeFrequency`](@ref) `laser` (of either rank). Channel positions are
exact at-field offsets from the pair's zero-field interval; for the reference
pair that is the very reference the laser offset is measured against, for the
others the reference cancels out of the driven-mode differences it is used in.
"""
function resonant_channels(species, basis::StateBasis{NoHyperfineNumberSpec}, laser, B)
    channels = [channel_type(basis)[] for _ in basis]
    pairs = quadrupole_pairs(species, unique(basis.levels))
    if laser isa RelativeFrequency
        named = (fine_structure(laser.lower), fine_structure(laser.upper))
        named in pairs || push!(pairs, named)
    end
    for (lo, hi) in pairs
        rank = multipole_rank(lo, hi)
        scale = channel_scale(species, lo, hi, rank)
        γ = uconvert(u"µs^-1", 1 / lifetime(species, hi))
        # Relative amplitude of the m_lo → m_hi component (zero also covers the
        # dipole-forbidden j combinations).
        amplitude(m_lo, q, m_hi) =
            rank == 1 ? dipole_cg(lo.j, m_lo, q, hi.j) :
            Float64(clebschgordan(lo.j, m_lo, 2, q, hi.j, m_hi))
        for (i, state) in enumerate(basis)
            if state.level == lo
                for m in (-hi.j):hi.j
                    q = m - state.m
                    abs(q) <= rank || continue
                    amp = amplitude(state.m, q, m)
                    iszero(amp) && continue
                    partner = StateSpec(hi, m)
                    Δ = uconvert(
                        u"µs^-1",
                        zeeman_shift(species, partner, B) -
                        zeeman_shift(species, state, B),
                    )
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
            elseif state.level == hi
                for m in (-lo.j):lo.j
                    q = state.m - m
                    abs(q) <= rank || continue
                    amp = amplitude(m, q, state.m)
                    iszero(amp) && continue
                    partner = StateSpec(lo, m)
                    Δ = uconvert(
                        u"µs^-1",
                        zeeman_shift(species, state, B) -
                        zeeman_shift(species, partner, B),
                    )
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

function resonant_channels(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    laser,
    B,
)
    channels = [channel_type(basis)[] for _ in basis]
    pairs = quadrupole_pairs(species, unique!(fine_structure.(basis.levels)))
    named = nothing
    if laser isa RelativeFrequency
        named = (fine_structure(laser.lower), fine_structure(laser.upper))
        named in pairs || push!(pairs, named)
    end
    needed = unique!([level for pair in pairs for level in pair])
    manifolds = Dict(fs => hyperfine_manifold(species, fs, B) for fs in needed)
    for (lo, hi) in pairs
        rank = multipole_rank(lo, hi)
        scale = channel_scale(species, lo, hi, rank)
        γ = uconvert(u"µs^-1", 1 / lifetime(species, hi))
        m_lo, m_hi = manifolds[lo], manifolds[hi]
        rotated = rotated_multipole_amplitudes(species, m_lo, m_hi, rank)
        # For the reference pair, positions are measured from the named
        # zero-field F-pair interval; the manifold eigen-energies are relative
        # to their centroids, so only the Casimir offsets of the named levels
        # remain (the centroid interval cancels).
        ref = 0.0u"µs^-1"
        if (lo, hi) == named
            ref = uconvert(
                u"µs^-1",
                hyperfine_shift(species, laser.upper) -
                hyperfine_shift(species, laser.lower),
            )
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
                    Δ = uconvert(u"µs^-1", (m_hi.energies[k] - m_lo.energies[li]) - ref)
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
                    Δ = uconvert(u"µs^-1", (m_hi.energies[ui] - m_lo.energies[k]) - ref)
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

"""
A bare laser frequency closer than this to an explicit electric-dipole channel
would make the background silently miss the Zeeman/hyperfine structure of the
line (which only a [`RelativeFrequency`](@ref) reference resolves); beyond it,
neglecting that structure is a ``≲ 10^{-3}`` relative error.
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
            detuning = abs(abs(species.energies[upper] - e_level) - ħω)
            if detuning < BACKGROUND_GUARD_WINDOW
                cyclic = round(u"GHz", detuning / u"ħ" / 2π; digits=1)
                throw(
                    ArgumentError(
                        "The laser is within $cyclic of the '$fs' → '$upper' " *
                        "resonance, where the electric-dipole background misses " *
                        "the Zeeman structure of the line; give the laser as a " *
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

    "The static flux density the channels were computed at, if any."
    field::F
end

"""
    LightShiftCoefficients(species, basis::StateBasis, laser[; B])
    LightShiftCoefficients(species, levels_or_states::AbstractVector, laser[; B])

Precomputes the ac Stark shift data of every state in the given basis, for a
laser given as a wavelength, an angular frequency, or a
[`RelativeFrequency`](@ref).

All the atomic structure enters here, so that evaluating the shift for a
particular intensity and polarisation afterwards costs only a handful of
arithmetic operations — the intended use when fitting laser parameters against
many observed transition frequencies.

Two kinds of data are held. The far-detuned electric-dipole **background** —
the per-state polarisabilities at the laser frequency — requires
[`LevelPolarisability`](@ref) data; levels without it are admitted, but
evaluating their background raises. The **near-resonant channels** — per-state
couplings to the individual Zeeman components of a level pair the laser sits
close to, at their exact at-field positions — are built when the static flux
density `B` is given: for every electric-quadrupole pair within the basis
(enabling [`driven_light_shift`](@ref)), and for the reference pair of a
`RelativeFrequency` laser of either rank (additionally enabling the parked
single-state [`light_shift`](@ref) at the laser's stated offset). The channels
need only the Einstein A coefficient, not polarisability data.

A `RelativeFrequency` laser requires `B`. If its reference pair is an
electric-dipole one, the co-rotating part of that channel is moved from the
background into the channels (which requires it to be explicit in the
[`LevelPolarisability`](@ref) data rather than lumped into the static
remainder). A bare wavelength that falls within ~100 GHz of an explicit
electric-dipole channel is refused for the same reason: the background alone
would silently miss the Zeeman structure of the line.
"""
LightShiftCoefficients(
    species,
    basis::StateBasis{NoHyperfineNumberSpec},
    laser;
    B=nothing,
) = build_coefficients(species, basis, laser, B)

LightShiftCoefficients(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    laser;
    B=nothing,
) = build_coefficients(species, basis, laser, B)

LightShiftCoefficients(species, levels_or_states::AbstractVector, laser; kwargs...) =
    LightShiftCoefficients(species, StateBasis(levels_or_states), laser; kwargs...)

function build_coefficients(species, basis, laser, B)
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
    rows = [
        isnothing(level_polarisability(species, state.level)) ?
        fill(NaN * u"C*m^2/V", 3) :
        state_polarisabilities(species, state, ħω; exclude=exclude(state.level)) for
        state in basis
    ]
    validate_background_detunings(species, unique(basis.levels), ħω, exclude)
    channels =
        isnothing(B) ? [channel_type(basis)[] for _ in basis] :
        resonant_channels(species, basis, laser, B)
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

`δ` overrides the laser offset (relative to the same reference interval), so a
frequency sweep or fit can reuse one set of coefficients.

The change a parked beam makes to the splitting of a `lower => upper` pair — a
Ramsey-type measurement, with no component resonantly driven — is the
difference of the two states' shifts; the displacement of the resonance
observed when *driving* a component is [`driven_light_shift`](@ref) instead.

The near-resonant sum is second-order perturbation theory: every channel Rabi
frequency must stay small against its detuning, and detunings within a few
linewidths of a pole are outside the model; both conditions warn when
violated. The one-shot species form is a convenience wrapper that rebuilds the
coefficients on every call.
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
        # A bare laser is guaranteed far from every explicit channel (cf.
        # LightShiftCoefficients), so the background is the whole story.
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
    if !(laser isa RelativeFrequency) && !isnothing(B)
        throw(
            ArgumentError(
                "B only enters through the near-resonant channels of a " *
                "RelativeFrequency laser; the shift of a single state from a " *
                "bare laser frequency does not depend on it",
            ),
        )
    end
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

Returns the light shift of the observed resonance when the `lower => upper`
Zeeman component is resonantly driven (in angular units): the amount by which
the frequency at which resonant Rabi flopping is observed exceeds the
unperturbed transition frequency, with the laser servo-locked to the probed
component.

The **background** part is the difference of the two states' far-detuned
electric-dipole shifts. The **resonant** part sums, for both states, the
couplings to every *other* Zeeman component sharing one of them — the driven
channel itself is the coherent drive, not a shift — with the laser pinned to
the probed component, so the laser frequency drops out and only exact at-field
splittings enter the detunings. `parts` selects `:total` (the default),
`:background` or `:resonant`.

This is the shift ⁸⁸Sr⁺ clock evaluations quote as the "E2 ac Stark shift"
(`[Lindvall2025]`, Sec. III F 2): with a perfectly linear polarisation the
shifts of the two components of a ``±m`` Zeeman pair are equal and opposite
and cancel in the pair average, whereas an elliptical polarisation leaves a
net shift. (For a hyperfine species the cancellation survives only
approximately — spectators in other ``F`` levels sit at hyperfine-interval
detunings that are even under ``m → −m``; the rigorous mirror identity is
shift(−m pair, B) = shift(+m pair, −B).)

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
Returns the per-channel polarisabilities of the stretched state ``m = F`` of a
hyperfine level, from which the ``α_0``/``α_1``/``α_2`` decomposition is
extracted (``α_{σ^-} = α_0 - α_1/2 - α_2/2``, ``α_π = α_0 + α_2``,
``α_{σ^+} = α_0 + α_1/2 - α_2/2`` there, as the tensor angular factor is unity
at ``m = F``).
"""
stretched_channels(
    species::HyperfineOneElectronSpecies,
    spec::HyperfineNumberSpec,
    ħω,
) = state_polarisabilities(species, StateSpec(spec, spec.f), ħω)

function stretched_channels(species, spec, ħω)
    throw(
        ArgumentError(
            "Level '$spec' carries hyperfine structure, but the species has none",
        ),
    )
end

"""
    scalar_polarisability(species, level, laser)

Returns the dynamic scalar polarisability ``α_0(ω)`` of the given level at the
given laser wavelength or angular frequency.

Together with [`vector_polarisability`](@ref) and
[`tensor_polarisability`](@ref), this decomposes the polarisability of the
sublevel ``m`` for polarisation `ε` as

``α = α_0 + 𝒜 \\frac{m}{2j} α_1 +
\\frac{3 |ε_0|^2 - 1}{2} \\frac{3 m^2 - j(j+1)}{j(2j-1)} α_2``,

where ``𝒜 = |ε_{-1}|^2 - |ε_{+1}|^2`` is the degree of circular polarisation
about the quantisation axis and ``ε_0 = ε · ẑ``. [`light_shift`](@ref) evaluates
the equivalent sum directly, without going through this decomposition.

For a hyperfine ``F`` level (with ``j`` replaced by ``F`` in the decomposition),
the components are extracted from the per-channel polarisabilities of the
stretched state; cf. [`Levels.state_polarisabilities`](@ref) for the
approximations involved.
"""
function scalar_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        α = stretched_channels(species, parsed, ħω)
        return (α[1] + α[2] + α[3]) / 3
    end
    spec, data = polarisability_data(species, level)
    e_level = species.energies[spec]
    α = data.static_scalar
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        α += uconvert(u"C*m^2/V", d^2 / (3 * (2 * spec.j + 1)) * 2Δ / (Δ^2 - ħω^2))
    end
    α
end

"""
    vector_polarisability(species, level, laser)

Returns the dynamic vector polarisability ``α_1(ω)`` of the given level at the
given laser wavelength or angular frequency (see
[`scalar_polarisability`](@ref) for the convention).

The vector polarisability vanishes in the static limit, so the far-detuned
remainder lumped into [`LevelPolarisability`](@ref)`.static_scalar` contributes
nothing to it.
"""
function vector_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f > 0 || return 0.0u"C*m^2/V"
        α = stretched_channels(species, parsed, ħω)
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

Returns the dynamic tensor polarisability ``α_2(ω)`` of the given level at the
given laser wavelength or angular frequency (see
[`scalar_polarisability`](@ref) for the convention).

Zero for ``j ≤ 1/2``, which has no oriented sublevels to distinguish.
"""
function tensor_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f >= 1 || return 0.0u"C*m^2/V"
        α = stretched_channels(species, parsed, ħω)
        return (2α[2] - α[1] - α[3]) / 3
    end
    spec, data = polarisability_data(species, level)
    j = spec.j
    j > 1//2 || return 0.0u"C*m^2/V"
    e_level = species.energies[spec]
    α = data.static_tensor
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        # The π-polarised shift of the stretched state m = j is exactly
        # α_0 + α_2, as the tensor angular factor is unity there.
        weight = dipole_cg(j, j, 0, upper.j)^2 / (2 * upper.j + 1) - 1 / (3 * (2j + 1))
        α += uconvert(u"C*m^2/V", d^2 * weight * 2Δ / (Δ^2 - ħω^2))
    end
    α
end

export LightShiftCoefficients,
    driven_light_shift,
    light_shift,
    scalar_polarisability,
    tensor_polarisability,
    vector_polarisability
public polarisation_weights, quadrupole_weights
