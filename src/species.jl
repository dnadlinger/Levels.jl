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
const POLARISABILITY_AU = uconvert(u"C*m^2/V", DIPOLE_AU^2 / (2 * u"R∞" * u"h" * u"c"))

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
Returns the Einstein A coefficient of the electric-dipole transition with
(angular) frequency `ω`, upper-level angular momentum `j_upper` and absolute
reduced dipole matrix element `d`,
``A = ω^3 d^2 / (3 π ε_0 ħ c^3 (2 j_\\mathrm{upper} + 1))``.
"""
einstein_a_from_dipole(d, ω, j_upper) =
    uconvert(u"µs^-1", ω^3 * d^2 / (3π * u"ε0" * u"ħ" * u"c"^3 * (2 * j_upper + 1)))

"""
Returns the absolute reduced dipole matrix element implied by the Einstein A
coefficient `a` (the inverse of [`Levels.einstein_a_from_dipole`](@ref)).
"""
dipole_from_einstein_a(a, ω, j_upper) =
    uconvert(u"C*m", sqrt(3π * u"ε0" * u"ħ" * u"c"^3 * (2 * j_upper + 1) * a / ω^3))

"""
The strength of an electric-dipole transition, specified through the absolute
reduced dipole matrix element ``|⟨\\mathrm{upper}‖d‖\\mathrm{lower}⟩|`` rather
than the Einstein A coefficient itself.

Accepted as a value in the `einstein_as` table by the species keyword
constructors, which convert it to the equivalent Einstein A coefficient
([`Levels.einstein_a_from_dipole`](@ref), with the frequency from the level
energies) — for pairs whose best-determined literature quantity is the matrix
element (typically from a many-body calculation) rather than a directly
measured rate. Downstream, only the converted rate exists; there is one
strength per transition, however it was specified.
"""
struct ReducedDipole
    "The absolute reduced dipole matrix element."
    d::typeof(1.0u"C*m")
end

"""
    ReducedDipole(d)

Creates the strength specification from a matrix element in any dipole-moment
units, so literature values can be written as e.g.
`ReducedDipole(4.187 * Levels.DIPOLE_AU)`.
"""
ReducedDipole(d::Quantity) = ReducedDipole(uconvert(u"C*m", d))

resolve_rate(rate, pair, energies) = rate

function resolve_rate(spec::ReducedDipole, (lower, upper)::Tuple, energies)
    Δ = energies[upper] - energies[lower]
    if Δ <= zero(Δ)
        throw(
            ArgumentError(
                "'$lower' is not below '$upper'; einstein_as keys are " *
                "(lower, upper) pairs",
            ),
        )
    end
    if multipole_rank(lower, upper) != 1 || abs(upper.j - lower.j) > 1
        throw(
            ArgumentError(
                "A ReducedDipole strength applies only to electric-dipole " *
                "pairs, which '$lower' → '$upper' is not",
            ),
        )
    end
    einstein_a_from_dipole(spec.d, Δ / u"ħ", upper.j)
end

"""
Resolves any [`ReducedDipole`](@ref) entries in an `einstein_as` table into
the equivalent Einstein A coefficients, passing rate entries through
unchanged.
"""
function resolve_einstein_as(einstein_as, energies)
    Dict(pair => resolve_rate(rate, pair, energies) for (pair, rate) in einstein_as)
end

"""
Electric-dipole light-shift data for one level, specified through the
quantities that are genuinely independent rather than the assembled
[`LevelPolarisability`](@ref) numbers: the intermediate levels whose channels
are treated explicitly, and the level's **total** static polarisabilities,
from which the lumped static remainder follows by subtracting the explicit
channels' static contributions.

A channel given as a bare level has its reduced dipole matrix element
**derived from the species' Einstein A coefficient** for the pair
([`Levels.dipole_from_einstein_a`](@ref)); a `level => dipole` pair carries an
explicit matrix element, and is accepted only for channels absent from the
`einstein_as` table (which must stay complete per upper level, as it defines
the level [`lifetime`](@ref)s) — a pair must never carry two independent
strengths.

The species constructors resolve this into a [`LevelPolarisability`](@ref) —
the one place the required `energies` and `einstein_as` tables are both in
hand — so a channel strength is stored exactly once per species, and an update
to a measured rate propagates into the light-shift background automatically.

Note that the resolved static remainder is a bookkeeping quantity rather than
a physical tail contribution: it is whatever the anchor totals exceed the
explicit channels by, and so also absorbs any tension between the channel
dipoles used here and the matrix elements underlying the anchors (cf. the
⁴³Ca⁺ S``_{1/2}`` remainder commentary in `species_data.jl`).
"""
struct ImplicitPolarisability
    """
    The intermediate levels whose channels are treated explicitly, each with
    an explicit reduced dipole matrix element, or with `nothing` for one
    derived from the species' Einstein A coefficient for the pair.
    """
    channels::Vector{Pair{NoHyperfineNumberSpec,Union{Nothing,typeof(1.0u"C*m")}}}

    """
    Total static scalar polarisability ``α_0(0)`` of the level, anchoring the
    lumped static remainder — or `nothing` if no total is known, leaving the
    explicit channels as the whole model (zero remainder).
    """
    static_scalar::Union{Nothing,typeof(1.0u"C*m^2/V")}

    """
    Total static tensor polarisability ``α_2(0)`` of the level (in the
    convention of [`tensor_polarisability`](@ref)), or `nothing` for no anchor
    as for `static_scalar`.
    """
    static_tensor::Union{Nothing,typeof(1.0u"C*m^2/V")}
end

"""
    ImplicitPolarisability(channels; static_scalar, static_tensor = nothing)

Creates the light-shift specification for one level from a collection of
intermediate levels — each a bare level (in spectroscopic notation or as a
[`NoHyperfineNumberSpec`](@ref); dipole derived from the Einstein A
coefficient) or a `level => dipole` pair — and the total static
polarisabilities.

All quantities are converted to the stored units, so literature values can be
written as e.g. `74.62 * Levels.POLARISABILITY_AU`. An anchor total given as
`nothing` keeps the corresponding static remainder at zero — for datasets
whose only known inputs are the explicit channels; `static_scalar` has no
default, so leaving out the anchor is always an explicit choice. (For a level
whose channels carry no tensor content — any ``j = 1/2`` one — the default
`static_tensor = nothing` is indistinguishable from a zero total.)
"""
function ImplicitPolarisability(channels; static_scalar, static_tensor=nothing)
    normalise(channel::Pair) =
        convert(NoHyperfineNumberSpec, channel.first) =>
            uconvert(u"C*m", channel.second)
    normalise(channel) = convert(NoHyperfineNumberSpec, channel) => nothing
    anchor(total) = isnothing(total) ? nothing : uconvert(u"C*m^2/V", total)
    ImplicitPolarisability(
        [normalise(c) for c in channels],
        anchor(static_scalar),
        anchor(static_tensor),
    )
end

resolve_polarisability(data::LevelPolarisability, level, energies, einstein_as) = data

"""
Resolves an [`ImplicitPolarisability`](@ref) specification for `level` into the
[`LevelPolarisability`](@ref) it implies, given the species' `energies` and
`einstein_as` tables (cf. the species keyword constructors).
"""
function resolve_polarisability(
    data::ImplicitPolarisability,
    level::NoHyperfineNumberSpec,
    energies,
    einstein_as,
)
    j = level.j
    e_level = energies[level]
    dipoles = Pair{NoHyperfineNumberSpec,typeof(1.0u"C*m")}[]
    channels_scalar = 0.0u"C*m^2/V"
    channels_tensor = 0.0u"C*m^2/V"
    for (channel, dipole) in data.channels
        if any(p -> p.first == channel, dipoles)
            throw(
                ArgumentError("Duplicate '$level' → '$channel' polarisability channel"),
            )
        end
        if !haskey(energies, channel)
            throw(
                ArgumentError(
                    "No known energy for the '$level' → '$channel' " *
                    "polarisability channel",
                ),
            )
        end
        if multipole_rank(level, channel) != 1 || abs(channel.j - j) > 1
            throw(
                ArgumentError(
                    "'$level' → '$channel' is not an electric-dipole pair, so " *
                    "it cannot be a polarisability channel",
                ),
            )
        end
        Δ = energies[channel] - e_level
        pair = Δ > zero(Δ) ? (level, channel) : (channel, level)
        d = if isnothing(dipole)
            a = get(einstein_as, pair, nothing)
            if isnothing(a)
                throw(
                    ArgumentError(
                        "No known Einstein A coefficient between '$level' and " *
                        "'$channel', which deriving that polarisability channel " *
                        "requires (a channel absent from einstein_as can be " *
                        "given as an explicit level => dipole pair)",
                    ),
                )
            end
            dipole_from_einstein_a(a, abs(Δ) / u"ħ", pair[2].j)
        else
            if haskey(einstein_as, pair)
                throw(
                    ArgumentError(
                        "An explicit dipole for '$level' → '$channel' would " *
                        "duplicate the strength its einstein_as entry already " *
                        "fixes; give the channel as a bare level to derive it " *
                        "from that",
                    ),
                )
            end
            dipole
        end
        push!(dipoles, channel => d)
        # The channels' static contributions, via the same helpers the
        # evaluation side uses (polarisability.jl — included later, but loaded
        # by species-construction time, cf. src/Levels.jl), so the resolved
        # data reproduces the anchor totals exactly.
        channels_scalar += channel_scalar_polarisability(d, j, Δ, zero(Δ))
        if j > 1//2
            channels_tensor +=
                channel_tensor_polarisability(d, j, channel.j, Δ, zero(Δ))
        end
    end
    # Whatever the anchor totals exceed the explicit channels by becomes the
    # lumped remainder; an unanchored (nothing) total leaves none at all.
    remainder(total, channels) = isnothing(total) ? zero(channels) : total - channels
    LevelPolarisability(
        dipoles;
        static_scalar=remainder(data.static_scalar, channels_scalar),
        static_tensor=remainder(data.static_tensor, channels_tensor),
    )
end

"""
Resolves any [`ImplicitPolarisability`](@ref) entries in a level →
polarisability-data dictionary, passing [`LevelPolarisability`](@ref) entries
through unchanged.
"""
function resolve_polarisabilities(polarisabilities, energies, einstein_as)
    resolved = Dict{NoHyperfineNumberSpec,LevelPolarisability}()
    for (level, data) in polarisabilities
        spec = convert(NoHyperfineNumberSpec, level)
        if !haskey(energies, spec)
            throw(
                ArgumentError(
                    "No known energy for level '$spec', for which " *
                    "polarisability data was given",
                ),
            )
        end
        resolved[spec] = resolve_polarisability(data, spec, energies, einstein_as)
    end
    resolved
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
struct NoHyperfineOneElectronSpecies{M<:Quantity,E<:Quantity,A<:Quantity} <:
       OneElectronSpecies
    """
    The mass of the ion.

    That of the actual charged species, i.e. the neutral-atom mass less one electron,
    plus the (tiny) mass equivalent of the ionisation energy. In practice, this is only
    a stylistic distinction; experimental or theoretical uncertainties in the derived
    quantities relevant to ion-trap AMO physics will dwarf the small difference in mass
    definitions.
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
    Measured electronic g-factors overriding the LS-coupling Landé formula in
    [`lande_g`](@ref), where available.
    """
    lande_g_overrides::Dict{NoHyperfineNumberSpec,Float64}

    """
    Light-shift data for the levels it is known for, if any.

    Levels missing from this dictionary have no [`light_shift`](@ref) defined.
    """
    polarisabilities::Dict{NoHyperfineNumberSpec,LevelPolarisability}
end

"""
    NoHyperfineOneElectronSpecies(; mass, energies, einstein_as,
        lande_g_overrides = Dict(), polarisabilities = Dict())

Creates the species from its level data.

`einstein_as` entries may be given as [`ReducedDipole`](@ref) matrix elements,
and `polarisabilities` may mix [`LevelPolarisability`](@ref) values with
[`ImplicitPolarisability`](@ref) specifications; both are resolved against the
level energies (and, for the latter, the resolved A coefficients) here.
"""
function NoHyperfineOneElectronSpecies(;
    mass,
    energies,
    einstein_as,
    lande_g_overrides=Dict{NoHyperfineNumberSpec,Float64}(),
    polarisabilities=Dict{NoHyperfineNumberSpec,LevelPolarisability}(),
)
    resolved_as = resolve_einstein_as(einstein_as, energies)
    NoHyperfineOneElectronSpecies(
        mass,
        energies,
        resolved_as,
        lande_g_overrides,
        resolve_polarisabilities(polarisabilities, energies, resolved_as),
    )
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
struct HyperfineOneElectronSpecies{
    M<:Quantity,
    E<:Quantity,
    A<:Quantity,
    H<:HyperfineConstants,
} <: OneElectronSpecies
    """
    The mass of the ion, as for [`NoHyperfineOneElectronSpecies`](@ref).
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
    lande_g_overrides::Dict{NoHyperfineNumberSpec,Float64}

    """
    Light-shift data for the levels it is known for, if any.

    Levels missing from this dictionary have no [`light_shift`](@ref) defined.
    """
    polarisabilities::Dict{NoHyperfineNumberSpec,LevelPolarisability}
end

"""
    HyperfineOneElectronSpecies(; mass, nuclear_spin, nuclear_g, energies, hyperfine,
        einstein_as, lande_g_overrides = Dict(), polarisabilities = Dict())

Creates the species from its level data.

As for [`NoHyperfineOneElectronSpecies`](@ref), `einstein_as` entries may be
given as [`ReducedDipole`](@ref) matrix elements, and `polarisabilities` may
mix [`LevelPolarisability`](@ref) values with [`ImplicitPolarisability`](@ref)
specifications; both are resolved here.
"""
function HyperfineOneElectronSpecies(;
    mass,
    nuclear_spin,
    nuclear_g,
    energies,
    hyperfine,
    einstein_as,
    lande_g_overrides=Dict{NoHyperfineNumberSpec,Float64}(),
    polarisabilities=Dict{NoHyperfineNumberSpec,LevelPolarisability}(),
)
    resolved_as = resolve_einstein_as(einstein_as, energies)
    HyperfineOneElectronSpecies(
        mass,
        nuclear_spin,
        nuclear_g,
        energies,
        hyperfine,
        resolved_as,
        lande_g_overrides,
        resolve_polarisabilities(polarisabilities, energies, resolved_as),
    )
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

Returns the saturation intensity ``I_0`` of the transition between the two levels, in
the Oxford ion-trap group's convention (follwoing A. Steane's `[SATSUM]` note).

``I_0`` is the intensity at which the resonant lowest-order excitation rate on an
electric-dipole Zeeman component equals ``C^2 A``, with ``C`` the Clebsch–Gordan factor
and ``A`` the Einstein coefficient of the transition:
```math
I_0 = \\frac{ħ ω^3}{6 π c^2 τ} = \\frac{4 π^2 ħ c}{3 λ^3 τ},
```
where ``ω`` is the transition frequency and ``τ`` the **total** [`lifetime`](@ref) of
the upper level (all decay channels summed).

Equivalently, the [`rabi_frequency`](@ref) of a component driven at unit geometric
amplitude obeys ``Ω^2 = C^2 (I/I_0) A/τ``. For a closed transition (``A = 1/τ``) this
makes ``I_0`` **twice** the two-level saturation intensity
``I_\\mathrm{sat} = π h c Γ/(3 λ^3)`` of e.g. `[Foot2005]` §7.6.1: in
the saturation parameter ``s = I/I_\\mathrm{sat} = 2 Ω^2/Γ^2`` convention,
``I = I_0`` corresponds to ``s = 2``, i.e. ``Ω = Γ`` on the stretch component.

An error is raised if a level is unknown to the species, the levels are given in the
wrong energy order, or the upper level has no known decay channels; a tabulated A
coefficient for the pair itself is not required, as it does not enter the definition.

# References

- `[SATSUM]`: A. M. Steane, "A summary of saturation intensities, or: How to find
  scattering rates and Rabi frequencies quickly", internal Oxford ion-trap group note
  (2001-12-03).
- `[Foot2005]`: C. J. Foot, "Atomic Physics", Oxford University Press (2005);
  same convention in Metcalf & van der Straten, "Laser Cooling and Trapping"
  (1999), §2.4.
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
    ImplicitPolarisability,
    ReducedDipole,
    einstein_a,
    lifetime,
    level_polarisability,
    saturation_intensity
public transition_frequency, BOHR_RADIUS, DIPOLE_AU, POLARISABILITY_AU, level_energy
