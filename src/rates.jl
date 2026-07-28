using Unitful.DefaultSymbols
using WignerSymbols

function clebsch_gordan(lower::StateSpec, upper::StateSpec)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    # TODO: Verify transition type.
    sqrt(2 * hi.j + 1) *
    sum(wigner3j(lo.j, 1, hi.j, -lower.m, q, upper.m) for q = -1:1) |> Float64
end

function rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity)
    i0 = saturation_intensity(species, lower.level, upper.level)
    cg = clebsch_gordan(lower, upper)
    a = einstein_a(species, lower.level, upper.level)
    τ = lifetime(species, upper.level)
    uconvert(kHz, sqrt(intensity / i0 * cg^2 * a / τ))
end

export rabi_frequency
