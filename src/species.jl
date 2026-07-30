using Unitful

"""
Generic type representing an atomic species.
"""
abstract type Species end

"""
Bohr radius ``a_0 = 4 π ε_0 ħ^2 / (m_e e^2)``.
"""
const BOHR_RADIUS = uconvert(u"m", 4π * u"ε0" * u"ħ"^2 / (u"me" * u"q"^2))

"""
Atomic unit of the electric dipole moment, ``e a_0``.

Reduced matrix elements are conventionally quoted in these units in the
literature; see [`LevelPolarisability`](@ref).
"""
const DIPOLE_AU = uconvert(u"C*m", u"q" * BOHR_RADIUS)

"""
Atomic unit of the electric polarisability, ``e^2 a_0^2 / E_h``.

Polarisabilities are conventionally quoted in these units in the literature; see
[`LevelPolarisability`](@ref).
"""
const POLARIZABILITY_AU = uconvert(u"C*m^2/V", DIPOLE_AU^2 / (2 * u"R∞" * u"h" * u"c"))

"""
Electric-dipole data for the ac Stark shift (light shift) of one level.

The shift is evaluated as a sum over the intermediate levels an electric-dipole
transition connects to (cf. [`light_shift`](@ref)). The few channels that
dominate are given explicitly by their reduced matrix elements, so that their
detuning from the driving laser is accounted for; everything else — the ionic
core, and the many weak transitions to high-lying levels — is lumped into a pair
of static polarisabilities, which is a good approximation as long as the laser
is far detuned from those transitions.
"""
struct LevelPolarisability
    """
    Reduced electric-dipole matrix elements ``|⟨k‖d‖a⟩|`` from this level ``a``
    to the intermediate levels ``k`` treated explicitly.

    The energies of the levels `k` are taken from
    [`NoHyperfineOneElectronSpecies`](@ref)`.energies`, so every key must appear
    there as well.
    """
    reduced_dipoles::Dict{NoHyperfineNumberSpec,typeof(1.0u"C*m")}

    """
    Scalar polarisability of all the channels not listed in `reduced_dipoles`,
    in the static limit.
    """
    static_scalar::typeof(1.0u"C*m^2/V")

    """
    Tensor polarisability of all the channels not listed in `reduced_dipoles`,
    in the static limit (in the convention of [`tensor_polarisability`](@ref)).
    """
    static_tensor::typeof(1.0u"C*m^2/V")
end

"""
    LevelPolarisability(reduced_dipoles; static_scalar, static_tensor)

Creates the light-shift data for one level from a collection of
`level => reduced matrix element` pairs and the static remainder.

Levels may be given in spectroscopic notation, and all quantities are converted
to the stored units, so literature values can be written as e.g.
`3.078 * Levels.DIPOLE_AU`.
"""
function LevelPolarisability(
    reduced_dipoles;
    static_scalar=0.0u"C*m^2/V",
    static_tensor=0.0u"C*m^2/V",
)
    LevelPolarisability(
        Dict(
            convert(NoHyperfineNumberSpec, level) => uconvert(u"C*m", d) for
            (level, d) in reduced_dipoles
        ),
        uconvert(u"C*m^2/V", static_scalar),
        uconvert(u"C*m^2/V", static_tensor),
    )
end

"""
Atomic species with only one (relevant) electron – all configurations spin-1/2 –
and no hyperfine structure.
"""
Base.@kwdef struct NoHyperfineOneElectronSpecies{M<:Quantity,E<:Quantity,A<:Quantity} <:
                   Species
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

    """
    Light-shift data for the levels it is known for, if any.

    Levels missing from this dictionary have no [`light_shift`](@ref) defined.
    """
    polarisabilities::Dict{NoHyperfineNumberSpec,LevelPolarisability} =
        Dict{NoHyperfineNumberSpec,LevelPolarisability}()
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

    uconvert(u"ps^-1", (e2 - e1) / u"ħ")
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
    level_polarisability(species::NoHyperfineOneElectronSpecies, level)

Returns the [`LevelPolarisability`](@ref) light-shift data for the given level.

Returns `nothing` if no such data is known for the level.
"""
function level_polarisability(species::NoHyperfineOneElectronSpecies, level)
    get(species.polarisabilities, convert(NoHyperfineNumberSpec, level), nothing)
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
    NoHyperfineOneElectronSpecies,
    LevelPolarisability,
    einstein_a,
    lifetime,
    level_polarisability,
    saturation_intensity
public transition_frequency, BOHR_RADIUS, DIPOLE_AU, POLARIZABILITY_AU
