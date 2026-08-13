using Unitful

"""
    RelativeFrequency(lower => upper, offset)

A laser frequency given relative to the zero-field `lower` → `upper`
transition frequency — the reference point at which all Zeeman components of
the two levels are degenerate — plus `offset` (an angular frequency, e.g.
`2π * 1.5u"MHz"`).

The point of the parametrisation is precision: consumers such as the
near-resonant machinery of [`light_shift`](@ref) only ever use `offset`
together with exactly-known Zeeman and hyperfine splittings, so the tabulated
transition energy — and its uncertainty, typically tens of MHz, far coarser
than the Zeeman structure — drops out. Only converting to an absolute
frequency via [`Levels.photon_energy`](@ref)`(species, laser)` reintroduces
it.

For a [`HyperfineOneElectronSpecies`](@ref) the levels must be ``F`` levels
(`"S_1/2 F=4"`), whose Zeeman components are degenerate at zero field; the
reference interval then includes the zero-field hyperfine shifts of the two
levels. The reference pair itself need not be a drivable component — it is a
frequency marker, not a transition to be driven.
"""
struct RelativeFrequency{L<:LevelSpec,D<:Unitful.Frequency}
    "The lower level of the reference transition."
    lower::L

    "The upper level of the reference transition."
    upper::L

    "The angular-frequency offset of the laser from the zero-field interval."
    offset::D
end

function RelativeFrequency(pair::Pair, offset)
    lower = parse_level(pair.first)
    upper = parse_level(pair.second)
    if typeof(lower) != typeof(upper)
        throw(
            ArgumentError(
                "Reference levels must be of the same kind, " *
                "got '$(pair.first)' and '$(pair.second)'",
            ),
        )
    end
    if lower == upper
        throw(ArgumentError("Reference levels must differ"))
    end
    if !(offset isa Unitful.Frequency)
        throw(
            ArgumentError(
                "The offset must be an angular frequency (e.g. 2π * 1.5u\"MHz\")",
            ),
        )
    end
    RelativeFrequency(lower, upper, offset)
end

"""
    photon_energy(laser)
    photon_energy(species, laser)

Returns the photon energy ``ħ ω`` of a laser given by its wavelength, its
angular frequency, or a [`RelativeFrequency`](@ref); the latter requires the
species, as it is resolved against the tabulated transition energies (whose
uncertainty — typically tens of MHz — the absolute value therefore inherits).
"""
photon_energy(wavelength::Unitful.Length) =
    uconvert(u"J", 2π * u"ħ" * u"c" / wavelength)
photon_energy(frequency::Unitful.Frequency) = uconvert(u"J", u"ħ" * frequency)
photon_energy(laser::RelativeFrequency) = throw(
    ArgumentError(
        "Resolving a RelativeFrequency to an absolute photon energy requires " *
        "the atomic data: call photon_energy(species, laser)",
    ),
)
photon_energy(species, laser) = photon_energy(laser)
photon_energy(species, laser::RelativeFrequency) = uconvert(
    u"J",
    u"ħ" * (transition_frequency(species, laser.lower, laser.upper) + laser.offset),
)

export RelativeFrequency
public photon_energy
