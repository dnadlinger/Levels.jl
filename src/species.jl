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
Atomic species with only one (relevant) electron, i.e. all configurations
spin-1/2.

Concrete subtypes are [`NoHyperfineOneElectronSpecies`](@ref) and
[`HyperfineOneElectronSpecies`](@ref); the level-data accessors
([`Levels.transition_frequency`](@ref), [`einstein_a`](@ref), [`lifetime`](@ref),
[`saturation_intensity`](@ref), [`level_polarisability`](@ref)) are shared
between them, with fine-structure data keyed on the
[`NoHyperfineNumberSpec`](@ref) part of the queried level.
"""
abstract type OneElectronSpecies <: Species end

"""
Atomic species with only one (relevant) electron – all configurations spin-1/2 –
and no hyperfine structure.
"""
Base.@kwdef struct NoHyperfineOneElectronSpecies{M<:Quantity,E<:Quantity,A<:Quantity} <:
                   OneElectronSpecies
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
Hyperfine coupling constants of one fine-structure level, stored as energies
(i.e. ``h`` times the conventionally quoted frequency values).

The magnetic-dipole constant `a` and electric-quadrupole constant `b` enter the
level Hamiltonian as ``A \\, \\vec{I} ⋅ \\vec{J}`` plus the standard Casimir
quadrupole term (cf. [`hyperfine_shift`](@ref)); `b` must be zero unless both
``I > 1/2`` and ``J > 1/2``.
"""
struct HyperfineConstants{E<:Quantity}
    "Magnetic-dipole hyperfine constant ``A``, as an energy."
    a::E

    "Electric-quadrupole hyperfine constant ``B``, as an energy."
    b::E
end

"""
    HyperfineConstants(; a, b = zero(a))

Creates the hyperfine coupling constants of one level, promoting the two
energies to a common type.
"""
HyperfineConstants(; a, b=zero(a)) = HyperfineConstants(promote(a, b)...)

"""
Atomic species with only one (relevant) electron – all configurations spin-1/2 –
and hyperfine structure from a nuclear spin ``I``.

The level data (`energies`, `einstein_as`, `polarisabilities`) is keyed on the
fine-structure levels exactly as for [`NoHyperfineOneElectronSpecies`](@ref) —
energies are hyperfine centroids, and Einstein A coefficients are fine-structure
rates (each hyperfine sublevel decays at the full fine-structure rate, so there
are never ``F``-resolved entries). The hyperfine structure itself is described
by `nuclear_spin`, `nuclear_g` and the per-level `hyperfine` coupling constants.
"""
Base.@kwdef struct HyperfineOneElectronSpecies{
    M<:Quantity,
    E<:Quantity,
    A<:Quantity,
    H<:HyperfineConstants,
} <: OneElectronSpecies
    """
    The atomic mass.

    Conventionally that of the neutral atom, as for
    [`NoHyperfineOneElectronSpecies`](@ref).
    """
    mass::M

    "The nuclear spin quantum number ``I``."
    nuclear_spin::Rational{Int}

    """
    The nuclear g-factor ``g_I``, in Bohr magnetons and in the convention
    ``H_Z = μ_B B (g_J m_J + g_I m_I)``, i.e. with the electron-like sign
    absorbed such that ``μ_I = -g_I I μ_B``.

    This is the effective moment of the nucleus bound in the ion (not corrected
    for diamagnetic shielding), which is the appropriate value for the Zeeman
    Hamiltonian.
    """
    nuclear_g::Float64

    """
    Energies of the hyperfine centroids of the fine-structure levels.

    The point of reference is chosen arbitrarily.
    """
    energies::Dict{NoHyperfineNumberSpec,E}

    "Hyperfine coupling constants for each fine-structure level."
    hyperfine::Dict{NoHyperfineNumberSpec,H}

    """
    Einstein A coefficients for (lower, upper) pairs of fine-structure levels.

    The sum of all A coefficients for a given upper level is the reciprocal
    level lifetime, so this is the linewidth contribution in angular units.
    """
    einstein_as::Dict{Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec},A}

    """
    Measured electronic g-factors overriding the LS-coupling Landé formula in
    [`lande_g`](@ref), where available.
    """
    lande_g_overrides::Dict{NoHyperfineNumberSpec,Float64} =
        Dict{NoHyperfineNumberSpec,Float64}()

    """
    Light-shift data for the levels it is known for, if any.

    Levels missing from this dictionary have no [`light_shift`](@ref) defined.
    """
    polarisabilities::Dict{NoHyperfineNumberSpec,LevelPolarisability} =
        Dict{NoHyperfineNumberSpec,LevelPolarisability}()
end

"""
Checks that the ``F`` quantum number of the given hyperfine level is compatible
with the nuclear spin of the species, returning the level.
"""
function validate_hyperfine(
    species::HyperfineOneElectronSpecies,
    spec::HyperfineNumberSpec,
)
    i = species.nuclear_spin
    if spec.f < abs(i - spec.j) ||
       spec.f > i + spec.j ||
       !isinteger(spec.f - i - spec.j)
        throw(
            ArgumentError(
                "F = $(spec.f) is outside |I - J| … I + J for level '$spec' " *
                "with I = $i",
            ),
        )
    end
    spec
end

"""
Returns the energy of the given level relative to the species' reference point:
the stored (centroid) energy for a fine-structure level, plus the zero-field
hyperfine shift for a hyperfine level.
"""
level_energy(species::OneElectronSpecies, level::NoHyperfineNumberSpec) =
    species.energies[level]
level_energy(species::HyperfineOneElectronSpecies, level::HyperfineNumberSpec) =
    species.energies[fine_structure(level)] + u"ħ" * hyperfine_shift(species, level)
function level_energy(species::OneElectronSpecies, level::LevelSpec)
    throw(
        ArgumentError(
            "Level '$level' is not supported by a $(nameof(typeof(species)))",
        ),
    )
end

"""
    transition_frequency(species::OneElectronSpecies, lower, upper)

Returns the frequency of a given transition (in angular units).

For hyperfine levels (of a [`HyperfineOneElectronSpecies`](@ref)), the
zero-field hyperfine shifts of any ``F``-resolved levels are included on top of
the centroid splitting.

An error is raised if the two levels are not connected by a known transition or
if the levels are incorrectly ordered.
"""
function transition_frequency(species::OneElectronSpecies, lower, upper)
    # TODO: Nicer errors on wrong keys.
    e1 = level_energy(species, parse_level(lower))
    e2 = level_energy(species, parse_level(upper))

    if e1 > e2
        throw(ArgumentError("'$lower' is of higher energy than '$upper'"))
    end

    uconvert(u"ps^-1", (e2 - e1) / u"ħ")
end

"""
    einstein_a(species::OneElectronSpecies, lower, upper)

Returns the Einstein A coefficient for the given two levels.

The sum of all A coefficients for a given upper level is the reciprocal
level lifetime, so this is the linewidth contribution in angular units.
Hyperfine levels are resolved to their fine-structure part, as A coefficients
are fine-structure quantities (the ``F``-resolved rates are branching fractions
of them; cf. [`Levels.hyperfine_reduction`](@ref)).

Returns `nothing` if there is no decay from `upper` to `lower`.
"""
function einstein_a(species::OneElectronSpecies, lower, upper)
    get(species.einstein_as, (fine_structure(lower), fine_structure(upper)), nothing)
end

"""
    lifetime(species::OneElectronSpecies, level)

Returns the lifetime of the given level.

This is the reciprocal of the sum of all the decay rates from the level.
Hyperfine levels are resolved to their fine-structure part — every hyperfine
sublevel decays at the full fine-structure rate.

Returns `nothing` if there are no decay channels defined from the given level.
"""
function lifetime(species::OneElectronSpecies, level)
    nhns = fine_structure(level)
    as = [a for ((_, up), a) in species.einstein_as if up == nhns]
    isempty(as) && return nothing
    uconvert(u"s", 1 / sum(as))
end

"""
    level_polarisability(species::OneElectronSpecies, level)

Returns the [`LevelPolarisability`](@ref) light-shift data for the given level
(the fine-structure part for hyperfine levels).

Returns `nothing` if no such data is known for the level.
"""
function level_polarisability(species::OneElectronSpecies, level)
    get(species.polarisabilities, fine_structure(level), nothing)
end

"""
    saturation_intensity(species::OneElectronSpecies, lower, upper)

Returns the saturation intensity for the transition between the two levels.

An error is raised if the two levels are not connected by a known transition.
"""
function saturation_intensity(species::OneElectronSpecies, lower, upper)
    ω = transition_frequency(species, lower, upper)
    τ = lifetime(species, upper)
    uconvert(u"W/m^2", u"ħ" * ω^3 / (6 * π * τ * u"c"^2))
end

export Species,
    OneElectronSpecies,
    NoHyperfineOneElectronSpecies,
    HyperfineOneElectronSpecies,
    HyperfineConstants,
    LevelPolarisability,
    einstein_a,
    lifetime,
    level_polarisability,
    saturation_intensity
public transition_frequency, BOHR_RADIUS, DIPOLE_AU, POLARIZABILITY_AU, level_energy
