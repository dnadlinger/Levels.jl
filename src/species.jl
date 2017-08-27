using Unitful

"""
Generic type representing an atomic species.
"""
abstract type Species end

"""
Atomic species with only one (relevant) electron – all configurations spin-1/2 –
and no hyperfine structure.
"""
type NoHyperfineOneElectronSpecies <: Species
    """
    Energies for different levels.

    The point of reference is chosen arbitrarily.
    """
    energies::Dict{NoHyperfineNumberSpec, Quantity}

    """
    Einstein A coefficients for (lower, upper) pair of levels.

    The sum of all A coefficients for a given upper level is the reciprocal
    level lifetime, so this is the linewidth contribution in angular units.
    """
    einstein_as::Dict{Tuple{NoHyperfineNumberSpec,
        NoHyperfineNumberSpec}, Quantity}
end

"""
Returns the frequency of a given transition (in angular units).

An error is raised if the two levels are not connected by a known transition or
if the levels are incorrectly ordered.
"""
function transition_frequency(species::NoHyperfineOneElectronSpecies, lower, upper)
    # TODO: Nicer errors on wrong keys.
    e1 = species.energies[NoHyperfineNumberSpec(lower)]
    e2 = species.energies[NoHyperfineNumberSpec(upper)]

    if e1 > e2
        throw(ArgumentError("'$lower' is of higher energy than '$upper'"))
    end

    uconvert(THz, (e2 - e1) / u"ħ")
end

"""
Returns the Einstein A coefficient for the given two levels.

The sum of all A coefficients for a given upper level is the reciprocal
level lifetime, so this is the linewidth contribution in angular units.

Returns null if there is no decay from upper to lower.
"""
function einstein_a(species::NoHyperfineOneElectronSpecies, lower, upper)::Nullable{Quantity}
    get(species.einstein_as, (NoHyperfineNumberSpec(lower), NoHyperfineNumberSpec(upper)), Nullable{Quantity}())
end


"""
Returns the lifetime of the given level.

This is the reciprocal of the sum for all the decay rates from the level.

Returns null if there are no decay channels defined from the given level.
"""
function lifetime(species::NoHyperfineOneElectronSpecies, level)::Nullable{Quantity}
    nhns = NoHyperfineNumberSpec(level)
    as = values(filter((l, v) -> l[2] == nhns, species.einstein_as))
    if isempty(as)
        return Nullable{Quantity}()
    end
    uconvert(s, 1 / sum(as))
end

"""
Returns the saturation intensity for the transition between the two levels.

An error is raised if the two levels are not connected by a known transition.
"""
function saturation_intensity(species::NoHyperfineOneElectronSpecies, lower, upper)
    ω = transition_frequency(species, lower, upper)
    τ = get(lifetime(species, upper))
    uconvert(W / m^2, u"ħ" * ω^3 / (6 * π * τ * u"c"^2))
end

function clebsch_gordan(species::NoHyperfineOneElectronSpecies, lower, upper)
end

export Species, NoHyperfineOneElectronSpecies, einstein_a, lifetime, saturation_intensity
