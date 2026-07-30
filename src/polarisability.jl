using Unitful
using WignerSymbols

"""
    photon_energy(laser)

Returns the photon energy ``ħ ω`` of a laser given either by its wavelength or
by its angular frequency.
"""
photon_energy(wavelength::Unitful.Length) =
    uconvert(u"J", 2π * u"ħ" * u"c" / wavelength)
photon_energy(frequency::Unitful.Frequency) = uconvert(u"J", u"ħ" * frequency)

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
Returns the polarisability of `state` for each of the three ``Δm`` channels
separately, as a length-3 vector indexed by `q + 2`.
"""
function state_polarisabilities(species, state::StateSpec, ħω)
    level = convert(NoHyperfineNumberSpec, state.level)
    data = level_polarisability(species, level)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$(state.level)'"))
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
        for q in -1:1
            rotating = dipole_cg(j, m, q, upper.j)^2
            counter = dipole_cg(j, m, -q, upper.j)^2
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

"""
Light-shift data for a set of states, precomputed for one laser frequency.

Construct with [`LightShiftCoefficients`](@ref)`(species, basis, laser)` and
evaluate with [`light_shift`](@ref); see there for the sign and unit
conventions.
"""
struct LightShiftCoefficients{L<:LevelSpec,E<:Quantity,T<:Quantity}
    "The states the coefficients refer to, fixing the row order."
    basis::StateBasis{L}

    "The photon energy the coefficients were computed for."
    photon_energy::E

    """
    Polarisability of each state for each ``Δm`` channel, as a
    `length(basis) × 3` matrix indexed by `[stateindex, q + 2]`.
    """
    polarisabilities::Matrix{T}
end

"""
    LightShiftCoefficients(species, basis::StateBasis, laser)
    LightShiftCoefficients(species, levels_or_states::AbstractVector, laser)

Precomputes the ac Stark shift of every state in the given basis for a laser of
the given wavelength or angular frequency.

All the atomic structure enters here, so that evaluating the shift for a
particular intensity and polarisation afterwards costs only a handful of
arithmetic operations — the intended use when fitting laser parameters against
many observed transition frequencies.

An error is raised if the species has no [`LevelPolarisability`](@ref) data for
one of the levels involved.
"""
function LightShiftCoefficients(species, basis::StateBasis, laser)
    ħω = photon_energy(laser)
    rows = [state_polarisabilities(species, state, ħω) for state in basis]
    LightShiftCoefficients(basis, ħω, permutedims(reduce(hcat, rows)))
end

LightShiftCoefficients(species, levels_or_states::AbstractVector, laser) =
    LightShiftCoefficients(species, StateBasis(levels_or_states), laser)

"""
Returns the polarisability of the basis state with the given index for the
[`polarisation_weights`](@ref) `w`.
"""
state_polarisability(c::LightShiftCoefficients, index::Integer, w) =
    w[1] * c.polarisabilities[index, 1] +
    w[2] * c.polarisabilities[index, 2] +
    w[3] * c.polarisabilities[index, 3]

"""
Converts a polarisability into the corresponding shift at the given intensity.
"""
shift_at_intensity(α, intensity) =
    uconvert(u"µs^-1", -intensity * α / (2 * u"c" * u"ε0" * u"ħ"))

"""
    light_shift(coefficients::LightShiftCoefficients, state::StateSpec, intensity, ε)
    light_shift(coefficients::LightShiftCoefficients, transition::Pair, intensity, ε)
    light_shift(coefficients::LightShiftCoefficients, transitions::AbstractVector, intensity, ε)

Returns the ac Stark (light) shift of a state, or of a `lower => upper`
transition, from the precomputed `coefficients` (in angular units).

For a transition this is the shift of the upper state minus that of the lower
one, i.e. the amount by which the frequency at which resonant Rabi flopping is
observed exceeds the unperturbed transition frequency. Given a vector of
transitions, the shifts are returned in the same order, sharing the work of
resolving the polarisation.

The beam direction does not enter: an electric-dipole light shift depends only
on the polarisation `ε`, whose overall scale is irrelevant as the intensity is
given separately.
"""
light_shift(c::LightShiftCoefficients, state::StateSpec, intensity, ε) =
    shift_at_intensity(
        state_polarisability(c, stateindex(c.basis, state), polarisation_weights(ε)),
        intensity,
    )

function light_shift(c::LightShiftCoefficients, transition::Pair, intensity, ε)
    w = polarisation_weights(ε)
    lower = stateindex(c.basis, transition.first)
    upper = stateindex(c.basis, transition.second)
    shift_at_intensity(
        state_polarisability(c, upper, w) - state_polarisability(c, lower, w),
        intensity,
    )
end

function light_shift(
    c::LightShiftCoefficients,
    transitions::AbstractVector{<:Pair},
    intensity,
    ε,
)
    w = polarisation_weights(ε)
    map(transitions) do transition
        lower = stateindex(c.basis, transition.first)
        upper = stateindex(c.basis, transition.second)
        shift_at_intensity(
            state_polarisability(c, upper, w) - state_polarisability(c, lower, w),
            intensity,
        )
    end
end

"""
    light_shift(species, state::StateSpec, laser, intensity, ε)
    light_shift(species, lower::StateSpec, upper::StateSpec, laser, intensity, ε)

Returns the ac Stark (light) shift of a single state or transition for a laser
of the given wavelength or angular frequency, intensity and polarisation.

This is a convenience wrapper that does the full sum over intermediate levels on
every call; use [`LightShiftCoefficients`](@ref) to hoist that work out of a
loop over intensities or polarisations.
"""
light_shift(species, state::StateSpec, laser, intensity, ε) = shift_at_intensity(
    sum(
        polarisation_weights(ε) .*
        state_polarisabilities(species, state, photon_energy(laser)),
    ),
    intensity,
)

function light_shift(species, lower::StateSpec, upper::StateSpec, laser, intensity, ε)
    ħω = photon_energy(laser)
    w = polarisation_weights(ε)
    α_lower = sum(w .* state_polarisabilities(species, lower, ħω))
    α_upper = sum(w .* state_polarisabilities(species, upper, ħω))
    shift_at_intensity(α_upper - α_lower, intensity)
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
"""
function scalar_polarisability(species, level, laser)
    ħω = photon_energy(laser)
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
    light_shift, scalar_polarisability, tensor_polarisability, vector_polarisability
public photon_energy, polarisation_weights
