using Unitful

"""
Generic type representing an atomic species.
"""
abstract type Species end

"""
Atomic species with only one (relevant) electron – all configurations spin-1/2 –
and no hyperfine structure.
"""
struct NoHyperfineOneElectronSpecies{M<:Quantity,E<:Quantity,A<:Quantity} <: Species
    """
    The atomic mass.

    Conventionally that of the neutral atom; the difference to the bare ion mass
    is negligible for the current uses (e.g. the reduced-mass correction in
    [`lande_g`](@ref) is only affected at the ``10^{-10}`` level).
    """
    mass::M

    """
    Energies for different levels.

    The point of reference is chosen arbitrarily.
    """
    energies::Dict{NoHyperfineNumberSpec,E}

    """
    Einstein A coefficients for (lower, upper) pair of levels.

    The sum of all A coefficients for a given upper level is the reciprocal
    level lifetime, so this is the linewidth contribution in angular units.
    """
    einstein_as::Dict{Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec},A}
end

"""
    transition_frequency(species::NoHyperfineOneElectronSpecies, lower, upper)

Returns the frequency of a given transition (in angular units).

An error is raised if the two levels are not connected by a known transition or
if the levels are incorrectly ordered.
"""
function transition_frequency(species::NoHyperfineOneElectronSpecies, lower, upper)
    # TODO: Nicer errors on wrong keys.
    e1 = species.energies[convert(NoHyperfineNumberSpec, lower)]
    e2 = species.energies[convert(NoHyperfineNumberSpec, upper)]

    if e1 > e2
        throw(ArgumentError("'$lower' is of higher energy than '$upper'"))
    end

    uconvert(u"THz", (e2 - e1) / u"ħ")
end

"""
    einstein_a(species::NoHyperfineOneElectronSpecies, lower, upper)

Returns the Einstein A coefficient for the given two levels.

The sum of all A coefficients for a given upper level is the reciprocal
level lifetime, so this is the linewidth contribution in angular units.

Returns `nothing` if there is no decay from `upper` to `lower`.
"""
function einstein_a(species::NoHyperfineOneElectronSpecies, lower, upper)
    get(
        species.einstein_as,
        (convert(NoHyperfineNumberSpec, lower), convert(NoHyperfineNumberSpec, upper)),
        nothing,
    )
end

"""
    lifetime(species::NoHyperfineOneElectronSpecies, level)

Returns the lifetime of the given level.

This is the reciprocal of the sum of all the decay rates from the level.

Returns `nothing` if there are no decay channels defined from the given level.
"""
function lifetime(species::NoHyperfineOneElectronSpecies, level)
    nhns = convert(NoHyperfineNumberSpec, level)
    as = [a for ((_, up), a) in species.einstein_as if up == nhns]
    isempty(as) && return nothing
    uconvert(u"s", 1 / sum(as))
end

"""
    saturation_intensity(species::NoHyperfineOneElectronSpecies, lower, upper)

Returns the saturation intensity for the transition between the two levels.

An error is raised if the two levels are not connected by a known transition.
"""
function saturation_intensity(species::NoHyperfineOneElectronSpecies, lower, upper)
    ω = transition_frequency(species, lower, upper)
    τ = lifetime(species, upper)
    uconvert(u"W/m^2", u"ħ" * ω^3 / (6 * π * τ * u"c"^2))
end

export Species,
    NoHyperfineOneElectronSpecies, einstein_a, lifetime, saturation_intensity
public transition_frequency
