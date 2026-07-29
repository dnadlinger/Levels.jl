using Unitful
using WignerSymbols

"""
    clebsch_gordan(lower::StateSpec, upper::StateSpec)

Returns the Clebsch–Gordan coefficient coupling the two given states, relative to
the reduced matrix element of the transition.
"""
function clebsch_gordan(lower::StateSpec, upper::StateSpec)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    # TODO: Verify transition type.
    Float64(
        sqrt(2 * hi.j + 1) *
        sum(wigner3j(lo.j, 1, hi.j, -lower.m, q, upper.m) for q in -1:1),
    )
end

"""
    rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity)

Returns the Rabi frequency (in angular units) for driving the given transition at
the given optical intensity.
"""
function rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity)
    i0 = saturation_intensity(species, lower.level, upper.level)
    cg = clebsch_gordan(lower, upper)
    a = einstein_a(species, lower.level, upper.level)
    τ = lifetime(species, upper.level)
    uconvert(u"kHz", sqrt(intensity / i0 * cg^2 * a / τ))
end

export rabi_frequency
public clebsch_gordan
