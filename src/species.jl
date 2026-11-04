using Unitful

"""
Generic type representing an atomic species.
"""
abstract type Species end

function normalise_level_table(entries, value=identity)
    pairs = [convert(NoHyperfineNumberSpec, k) => value(v) for (k, v) in entries]
    allunique(first.(pairs)) ||
        throw(ArgumentError("Duplicate canonical level in data"))
    Dict(pairs)
end

function validate_energy(value)
    (value isa Unitful.Energy && isreal(value) && isfinite(value)) ||
        throw(ArgumentError("An energy must be a real, finite energy, got $value"))
    float(uconvert(u"J", value))
end

function validate_mass(mass)
    (mass isa Unitful.Mass && isreal(mass) && isfinite(mass) && mass > zero(mass)) ||
        throw(ArgumentError("The species mass must be a positive, finite mass"))
end

function normalise_rate_table(entries, energies)
    result = Dict{Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},Any}()
    for (pair, rate) in entries
        lo, hi = convert.(NoHyperfineNumberSpec, Tuple(pair))
        for level in (lo, hi)
            haskey(energies, level) ||
                throw(ArgumentError("Unknown level '$level' in rate data"))
        end
        energies[lo] < energies[hi] ||
            throw(ArgumentError("Rate keys must be lower => upper"))
        if !(rate isa Union{ElectricDipole,MagneticDipole})
            (
                rate isa Unitful.Frequency &&
                isreal(rate) &&
                isfinite(rate) &&
                rate >= zero(rate)
            ) || throw(
                ArgumentError("An Einstein A must be a non-negative, finite rate"),
            )
            # A bare rate is taken as the lowest electric multipole allowed by
            # parity, i.e. E2 here, whereas fine-structure partners decay by M1.
            lo.l == hi.l && throw(
                ArgumentError(
                    "A plain rate between the fine-structure levels '$lo' and " *
                    "'$hi' would be taken as electric-quadrupole; give their " *
                    "magnetic-dipole decay as a MagneticDipole",
                ),
            )
        end
        haskey(result, lo => hi) &&
            throw(ArgumentError("Duplicate transition '$lo' => '$hi'"))
        result[lo=>hi] = rate
    end
    result
end

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

Explicit intermediate levels contribute their frequency-dependent E1 sums;
the ionic core and omitted channels contribute a static scalar/tensor
remainder. This approximation requires the driving field to be far from the
omitted resonances. Used by [`polarisability`](@ref) and [`light_shift`](@ref).
"""
struct LevelPolarisability
    """
    Reduced electric-dipole matrix elements ``|⟨k‖d‖a⟩|`` from this level ``a``
    to the intermediate levels ``k`` treated explicitly.

    Keys are fine-structure levels and must also appear in the species'
    `energies` table.
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

Creates light-shift data from `level => reduced matrix element` pairs and
static **remainders**, not total polarisabilities. Both remainders default
to zero; [`ImplicitPolarisability`](@ref) accepts total static anchors instead.

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
Returns the Einstein A coefficient of the electric-dipole (E1) transition with
(angular) frequency `ω`, upper-level angular momentum `j_upper` and absolute
reduced dipole matrix element `d`,
``A = ω^3 d^2 / (3 π ε_0 ħ c^3 (2 j_\\mathrm{upper} + 1))``.
"""
einstein_a_from_electric_dipole(d, ω, j_upper) =
    uconvert(u"µs^-1", ω^3 * d^2 / (3π * u"ε0" * u"ħ" * u"c"^3 * (2 * j_upper + 1)))

"""
Returns the absolute reduced dipole matrix element implied by the Einstein A
coefficient `a` (the inverse of [`Levels.einstein_a_from_electric_dipole`](@ref)).
"""
electric_dipole_from_einstein_a(a, ω, j_upper) =
    uconvert(u"C*m", sqrt(3π * u"ε0" * u"ħ" * u"c"^3 * (2 * j_upper + 1) * a / ω^3))

# Validates the argument of a multipole strength specification (`type_name`) and
# converts it to the stored units: µs⁻¹ for an Einstein A coefficient, or
# `moment_unit` for a reduced matrix element (described as `moment_description`).
function multipole_strength(type_name, strength, moment_unit, moment_description)
    if !(isreal(strength) && isfinite(strength) && strength >= zero(strength))
        throw(
            ArgumentError(
                "A $type_name strength must be non-negative and finite, got $strength",
            ),
        )
    end
    dimension(strength) == dimension(u"s^-1") &&
        return convert(typeof(1.0u"µs^-1"), strength)
    dimension(strength) == dimension(moment_unit) &&
        return convert(typeof(1.0 * moment_unit), strength)
    throw(
        ArgumentError(
            "A $type_name strength must be an Einstein A coefficient (an inverse " *
            "time) or a $moment_description, got $strength",
        ),
    )
end

"""
The strength of an electric-dipole (E1) transition, specified either as its
Einstein A coefficient or through the absolute reduced dipole matrix element
``|⟨\\mathrm{upper}‖d‖\\mathrm{lower}⟩|``.

Accepted as a value in the `einstein_as` table by the species keyword
constructors. A matrix element is converted to the equivalent rate
([`Levels.einstein_a_from_electric_dipole`](@ref), with the frequency from the
level energies): the form for pairs whose best-determined literature quantity is
the matrix element (typically from a many-body calculation) rather than a directly
measured rate. A rate is equivalent to a plain `einstein_as` entry, which parity
already marks as E1, but states the multipole explicitly. Downstream only the rate
exists; there is one strength per transition, however it was specified.
"""
struct ElectricDipole
    "The Einstein A coefficient, or the absolute reduced matrix element."
    strength::Union{typeof(1.0u"µs^-1"),typeof(1.0u"C*m")}

    ElectricDipole(strength::Quantity) = new(
        multipole_strength(
            "ElectricDipole",
            strength,
            u"C*m",
            "reduced dipole matrix element (e.g. in Levels.DIPOLE_AU)",
        ),
    )
end

@doc """
    ElectricDipole(strength)

Creates the E1 strength specification from an Einstein A coefficient (an inverse
time, e.g. `1.4e8u"s^-1"`) or a reduced matrix element in any dipole-moment units
(e.g. `4.187 * Levels.DIPOLE_AU`).
""" ElectricDipole(::Quantity)

"""
Returns the Einstein A coefficient of the magnetic-dipole (M1) transition with
(angular) frequency `ω`, upper-level angular momentum `j_upper` and absolute
reduced magnetic-moment matrix element `μ`,
``A = ω^3 μ^2 / (3 π ε_0 ħ c^5 (2 j_\\mathrm{upper} + 1))`` — the electric-dipole
expression of [`Levels.einstein_a_from_electric_dipole`](@ref) with ``d → μ/c``.
"""
einstein_a_from_magnetic_dipole(μ, ω, j_upper) =
    uconvert(u"µs^-1", ω^3 * μ^2 / (3π * u"ε0" * u"ħ" * u"c"^5 * (2 * j_upper + 1)))

"""
The strength of a magnetic-dipole (M1) transition, specified either as its
Einstein A coefficient or through the absolute reduced matrix element
``|⟨\\mathrm{upper}‖μ‖\\mathrm{lower}⟩|`` of the magnetic moment.

Accepted as a value in the `einstein_as` table by the species keyword
constructors, which record the pair as M1 (cf. [`Levels.multipole`](@ref)): a
same-parity pair is otherwise taken as an electric-quadrupole transition. The
typical use is the decay between fine-structure partners, such as
D``_{5/2}`` → D``_{3/2}``, for which a plain rate is refused as ambiguous. A matrix
element is converted to the equivalent rate
([`Levels.einstein_a_from_magnetic_dipole`](@ref), with the frequency from the
level energies). Downstream only that rate exists, and it fixes the decay
(lifetimes, jump operators) and the coherent coupling ([`rabi_frequency`](@ref))
of the pair alike.
"""
struct MagneticDipole
    "The Einstein A coefficient, or the absolute reduced matrix element."
    strength::Union{typeof(1.0u"µs^-1"),typeof(1.0u"J/T")}

    MagneticDipole(strength::Quantity) = new(
        multipole_strength(
            "MagneticDipole",
            strength,
            u"J/T",
            "reduced magnetic-moment matrix element (e.g. in μB)",
        ),
    )
end

@doc """
    MagneticDipole(strength)

Creates the M1 strength specification from an Einstein A coefficient (an inverse
time, e.g. `5.6e-3u"s^-1"`) or a reduced matrix element in magnetic-moment units
(e.g. `1.55u"μB"`).
""" MagneticDipole(::Quantity)

resolve_rate(rate, pair, energies) = rate

function resolve_rate(spec::ElectricDipole, (lower, upper)::Pair, energies)
    Δ = energies[upper] - energies[lower]
    if Δ <= zero(Δ)
        throw(
            ArgumentError(
                "'$lower' is not below '$upper'; einstein_as keys are " *
                "(lower, upper) pairs",
            ),
        )
    end
    if parity_multipole(lower => upper) != :E1 || abs(upper.j - lower.j) > 1
        throw(
            ArgumentError(
                "An ElectricDipole strength applies only to electric-dipole " *
                "pairs (opposite parity, |ΔJ| ≤ 1), which '$lower' → '$upper' is not",
            ),
        )
    end
    spec.strength isa Unitful.Frequency && return spec.strength
    einstein_a_from_electric_dipole(spec.strength, Δ / u"ħ", upper.j)
end

function resolve_rate(spec::MagneticDipole, (lower, upper)::Pair, energies)
    if parity_multipole(lower => upper) != :E2 || !is_triangle(lower.j, 1, upper.j)
        throw(
            ArgumentError(
                "A MagneticDipole strength applies only to magnetic-dipole " *
                "pairs (equal parity, |ΔJ| ≤ 1), which '$lower' → '$upper' is not",
            ),
        )
    end
    spec.strength isa Unitful.Frequency && return spec.strength
    Δ = energies[upper] - energies[lower]
    einstein_a_from_magnetic_dipole(spec.strength, Δ / u"ħ", upper.j)
end

"""
Resolves any [`ElectricDipole`](@ref) or [`MagneticDipole`](@ref) entries in an
`einstein_as` table into the equivalent Einstein A coefficients, passing rate
entries through unchanged.
"""
function resolve_einstein_as(einstein_as, energies)
    isempty(einstein_as) && return Dict{
        Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},
        typeof(1.0u"µs^-1"),
    }()
    Dict(
        pair => float(uconvert(u"µs^-1", resolve_rate(rate, pair, energies))) for
        (pair, rate) in einstein_as
    )
end

"""
Returns the multipolarity of each entry of a normalised `einstein_as` table: M1 for
a [`MagneticDipole`](@ref), otherwise the lowest electric multipole allowed by
parity ([`Levels.parity_multipole`](@ref)).
"""
tabulate_multipoles(einstein_as) =
    Dict{Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},Symbol}(
        pair => (rate isa MagneticDipole ? :M1 : parity_multipole(pair)) for
        (pair, rate) in einstein_as
    )

"""
E1 channels and **total** static polarisabilities for one fine-structure level.
Species constructors resolve these into [`LevelPolarisability`](@ref) data,
subtracting the explicit channels' static contributions from each total.

A bare channel level derives its dipole from the species' Einstein A
coefficient. An explicit `level => dipole` is accepted only if the pair is
absent from `einstein_as`, so each transition has one independent strength.
Rate tables should remain complete for each upper level because they also
determine [`lifetime`](@ref).

The remainder includes any mismatch between the chosen dipoles and the data
underlying the static anchors; it need not represent a physical tail alone.
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

Quantities accept any compatible units, e.g. `74.62 * Levels.POLARISABILITY_AU`.
An anchor of `nothing` leaves the corresponding remainder zero; an anchor of
zero makes the **total** zero. `static_scalar` is required, while
`static_tensor` defaults to `nothing`. For ``J ≤ 1/2``, the electronic tensor
remainder vanishes; resolved hyperfine channels can still give an ``F`` level
a tensor response.
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
        if parity_multipole(level => channel) != :E1 || abs(channel.j - j) > 1
            throw(
                ArgumentError(
                    "'$level' → '$channel' is not an electric-dipole pair, so " *
                    "it cannot be a polarisability channel",
                ),
            )
        end
        Δ = energies[channel] - e_level
        pair = Δ > zero(Δ) ? (level => channel) : (channel => level)
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
            electric_dipole_from_einstein_a(a, abs(Δ) / u"ħ", pair[2].j)
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
        # Shared evaluation helpers keep the static anchors consistent.
        channels_scalar += channel_scalar_polarisability(d, j, Δ, zero(Δ))
        if j > 1//2
            channels_tensor +=
                channel_tensor_polarisability(d, j, channel.j, Δ, zero(Δ))
        end
    end
    # An absent anchor leaves a channels-only model.
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

    The neutral-atom mass minus one electron mass, plus the ionisation energy
    divided by ``c^2``.
    """
    mass::M

    """
    Energies for different levels.

    The point of reference is chosen arbitrarily.
    """
    energies::Dict{NoHyperfineNumberSpec,E}

    """
    Einstein A coefficients for (lower, upper) pair of levels.

    Rates in inverse time, without a factor of 2π. Their sum for an upper
    level is its inverse lifetime and its natural angular linewidth.
    """
    einstein_as::Dict{Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},A}

    """
    Multipolarity (`:E1`, `:E2` or `:M1`) of each `einstein_as` entry, under the
    same keys (cf. [`Levels.multipole`](@ref)).
    """
    multipoles::Dict{Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},Symbol}

    """
    Measured electronic g-factors overriding the LS-coupling Landé formula in
    [`lande_g`](@ref), where available.
    """
    lande_g_overrides::Dict{NoHyperfineNumberSpec,Float64}

    """
    Light-shift data for the levels it is known for, if any.

    A background light shift requires these data; near-resonant channels
    remain available without them (cf. [`light_shift`](@ref)).
    """
    polarisabilities::Dict{NoHyperfineNumberSpec,LevelPolarisability}
end

"""
    NoHyperfineOneElectronSpecies(; mass, energies, einstein_as,
        lande_g_overrides = Dict(), polarisabilities = Dict())

Creates the species from its level data.

`einstein_as` entries may be plain Einstein A coefficients, taken as the lowest
electric multipole parity allows, or [`ElectricDipole`](@ref) and
[`MagneticDipole`](@ref) strengths, each a rate or a reduced matrix element; a
magnetic-dipole transition must be given as a `MagneticDipole`. `polarisabilities`
may mix [`LevelPolarisability`](@ref) values with
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
    validate_mass(mass)
    energies = normalise_level_table(energies, validate_energy)
    einstein_as = normalise_rate_table(einstein_as, energies)
    lande_g_overrides = Dict{NoHyperfineNumberSpec,Float64}(
        normalise_level_table(lande_g_overrides, Float64),
    )
    polarisabilities = normalise_level_table(polarisabilities)
    resolved_as = resolve_einstein_as(einstein_as, energies)
    NoHyperfineOneElectronSpecies(
        mass,
        energies,
        resolved_as,
        tabulate_multipoles(einstein_as),
        lande_g_overrides,
        resolve_polarisabilities(polarisabilities, energies, resolved_as),
    )
end

"""
    StateBasis(species::NoHyperfineOneElectronSpecies, levels...)

Creates a basis from the given fine-structure levels, each contributing all its
Zeeman sublevels in order of increasing ``m_J`` — the same as
[`StateBasis(levels...)`](@ref StateBasis), with the species argument letting
code construct bases for either species kind alike (for a
[`HyperfineOneElectronSpecies`](@ref), the corresponding method expands
fine-structure levels into their ``F`` levels).
"""
StateBasis(::NoHyperfineOneElectronSpecies, levels...) = StateBasis(levels...)

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
HyperfineConstants(; a, b=zero(a)) =
    HyperfineConstants(promote(validate_energy(a), validate_energy(b))...)

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
    The conventional signed nuclear g-factor ``g_I = μ_I/(I μ_N)``, in
    nuclear magnetons. The Zeeman convention is
    ``H_Z = B (g_J μ_B m_J - g_I μ_N m_I)``.

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

    Rates in inverse time, without a factor of 2π. Their sum for an upper
    level is its inverse lifetime and its natural angular linewidth.
    """
    einstein_as::Dict{Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},A}

    """
    Multipolarity (`:E1`, `:E2` or `:M1`) of each `einstein_as` entry, under the
    same keys (cf. [`Levels.multipole`](@ref)).
    """
    multipoles::Dict{Pair{NoHyperfineNumberSpec,NoHyperfineNumberSpec},Symbol}

    """
    Measured electronic g-factors overriding the LS-coupling Landé formula in
    [`lande_g`](@ref), where available.
    """
    lande_g_overrides::Dict{NoHyperfineNumberSpec,Float64}

    """
    Light-shift data for the levels it is known for, if any.

    A background light shift requires these data; near-resonant channels
    remain available without them (cf. [`light_shift`](@ref)).
    """
    polarisabilities::Dict{NoHyperfineNumberSpec,LevelPolarisability}
end

"""
    HyperfineOneElectronSpecies(; mass, nuclear_spin, nuclear_g, energies, hyperfine,
        einstein_as, lande_g_overrides = Dict(), polarisabilities = Dict())

Creates the species from its level data.

As for [`NoHyperfineOneElectronSpecies`](@ref), `einstein_as` entries may be
plain rates or [`ElectricDipole`](@ref) and [`MagneticDipole`](@ref) strengths,
and `polarisabilities` may mix [`LevelPolarisability`](@ref) values with
[`ImplicitPolarisability`](@ref) specifications; both are resolved here.
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
    validate_mass(mass)
    (nuclear_spin >= 0 && isinteger(2nuclear_spin)) || throw(
        ArgumentError("nuclear_spin must be a non-negative integer or half-integer"),
    )
    isfinite(nuclear_g) || throw(ArgumentError("nuclear_g must be finite"))
    energies = normalise_level_table(energies, validate_energy)
    einstein_as = normalise_rate_table(einstein_as, energies)
    hyperfine = normalise_level_table(hyperfine)
    lande_g_overrides = Dict{NoHyperfineNumberSpec,Float64}(
        normalise_level_table(lande_g_overrides, Float64),
    )
    polarisabilities = normalise_level_table(polarisabilities)
    resolved_as = resolve_einstein_as(einstein_as, energies)
    HyperfineOneElectronSpecies(
        mass,
        nuclear_spin,
        nuclear_g,
        energies,
        hyperfine,
        resolved_as,
        tabulate_multipoles(einstein_as),
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
The level kind of a species: [`NoHyperfineNumberSpec`](@ref) for a
fine-structure species, [`HyperfineNumberSpec`](@ref) for a hyperfine one.
"""
level_kind(::NoHyperfineOneElectronSpecies) = NoHyperfineNumberSpec
level_kind(::HyperfineOneElectronSpecies) = HyperfineNumberSpec

"""
Parses a level and checks it is of the species' kind — an ``F`` level for a
hyperfine species, a ``J`` level for a fine-structure one — so that `m` keeps
its meaning; returns the canonical spec.
"""
function validate_level_kind(species::OneElectronSpecies, level)
    spec = parse_level(level)
    spec isa level_kind(species) && return spec
    if species isa HyperfineOneElectronSpecies
        throw(
            ArgumentError(
                "Levels of a hyperfine species must specify hyperfine (F) " *
                "levels, e.g. 'S_1/2 F=4', not '$spec'",
            ),
        )
    end
    throw(
        ArgumentError(
            "Level '$spec' carries hyperfine structure, but the species has none",
        ),
    )
end

"""
Checks that a state's level is of the species' kind (cf.
[`Levels.validate_level_kind`](@ref)), returning the canonical state.
"""
validate_state_kind(species::OneElectronSpecies, state::StateSpec) =
    StateSpec(validate_level_kind(species, state.level), state.m)

"""
Checks that the levels of a basis are of the species' kind (cf.
[`Levels.validate_level_kind`](@ref)).
"""
function validate_basis_kind(species::OneElectronSpecies, basis::StateBasis)
    for level in basis.levels
        validate_level_kind(species, level)
    end
    basis
end

"""
Returns the energy of the given level relative to the species' reference point:
the stored (centroid) energy for a fine-structure level, plus the zero-field
hyperfine shift for a hyperfine level.
"""
level_energy(species::OneElectronSpecies, level::NoHyperfineNumberSpec) =
    get(species.energies, level) do
        throw(ArgumentError("Unknown level '$level'"))
    end
level_energy(species::HyperfineOneElectronSpecies, level::HyperfineNumberSpec) =
    level_energy(species, fine_structure(level)) +
    u"ħ" * hyperfine_shift(species, level)
level_energy(species::OneElectronSpecies, level::AbstractString) =
    level_energy(species, parse_level(level))

function level_energy(species::OneElectronSpecies, level::LevelSpec)
    throw(
        ArgumentError(
            "Level '$level' is not supported by a $(nameof(typeof(species)))",
        ),
    )
end

"""
    transition_frequency(species::OneElectronSpecies, lower => upper)

Returns the zero-field angular frequency of the given transition.

For hyperfine levels (of a [`HyperfineOneElectronSpecies`](@ref)), the
zero-field hyperfine shifts of any ``F``-resolved levels are included on top of
the centroid splitting. A pair of states gives the zero-field interval of
their levels (their Zeeman components are degenerate there); the interval at
a static field is [`transition_frequency(s::Eigenstructure, lower => upper)`](@ref).

Both levels must have tabulated energies; a decay channel between them is not
required. A pair of levels must be in increasing energy order. A pair of
states returns the signed upper-minus-lower interval, including within one
fine-structure manifold.
"""
function transition_frequency(species::OneElectronSpecies, transition::Pair)
    lower, upper = transition
    if lower isa StateSpec && upper isa StateSpec
        lo = validate_state_kind(species, lower)
        hi = validate_state_kind(species, upper)
        centroid =
            (
                level_energy(species, fine_structure(hi.level)) -
                level_energy(species, fine_structure(lo.level))
            ) / u"ħ"
        return uconvert(
            u"rad/µs",
            centroid + hyperfine_shift(species, hi.level) -
            hyperfine_shift(species, lo.level),
        )
    end
    e1 = level_energy(species, lower)
    e2 = level_energy(species, upper)
    e1 <= e2 || throw(ArgumentError("'$lower' is of higher energy than '$upper'"))
    uconvert(u"rad/µs", (e2 - e1) / u"ħ")
end

"""
    einstein_a(species::OneElectronSpecies, lower => upper)

Returns the Einstein A coefficient for the given two levels.

The coefficient is a decay rate in inverse time, without a factor of 2π.
Summing over all lower levels gives the inverse upper-level lifetime, equal
in value to its natural angular linewidth.
Hyperfine levels are resolved to their fine-structure part, as A coefficients
are fine-structure quantities (the ``F``-resolved rates are branching fractions
of them; cf. [`Levels.hyperfine_reduction`](@ref)).

Returns `nothing` if there is no decay from `upper` to `lower`.
"""
function einstein_a(species::OneElectronSpecies, transition::Pair)
    lower, upper = transition
    get(species.einstein_as, fine_structure(lower) => fine_structure(upper), nothing)
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
    saturation_intensity(species::OneElectronSpecies, lower => upper)

Returns the saturation intensity in the Oxford ion-trap group's convention
(`[SATSUM]`):

```math
I_0 = \\frac{ħ ω^3}{6π c^2 τ} = \\frac{4π^2 ħ c}{3λ^3 τ},
```

where ``ω`` is the pair's zero-field angular frequency and ``τ`` the total
upper-level [`lifetime`](@ref). A tabulated Einstein A coefficient for the
pair is not required; unknown levels, reversed energy order or missing
upper-level decay data raise an error.

For a fine-structure E1 pair at unit geometric amplitude,
``|Ω|^2 = C^2 (I/I_0) A/τ``, with ``C`` the Clebsch–Gordan coefficient.
On a closed stretched transition, ``I = I_0`` gives ``|Ω| = Γ = 1/τ``.
Thus ``I_0 = 2 I_\\mathrm{sat}`` for the two-level convention
``s = I/I_\\mathrm{sat} = 2|Ω|^2/Γ^2`` of `[Foot2005]`.

# References

- `[SATSUM]`: A. M. Steane, "A summary of saturation intensities, or: How to find
  scattering rates and Rabi frequencies quickly", internal Oxford ion-trap group note
  (2001-12-03).
- `[Foot2005]`: C. J. Foot, "Atomic Physics", Oxford University Press (2005),
  §7.6.1; also Metcalf & van der Straten, "Laser Cooling and Trapping"
  (1999), §2.4.
"""
function saturation_intensity(species::OneElectronSpecies, transition::Pair)
    lower, upper = transition
    ω = transition_frequency(species, lower => upper)
    τ = lifetime(species, upper)
    isnothing(τ) && throw(ArgumentError("No decay data for level '$upper'"))
    uconvert(u"W/m^2", u"ħ" * ω^3 / (6 * π * τ * u"c"^2))
end

"Returns the known fine-structure levels in increasing energy order."
levels(species::OneElectronSpecies) =
    sort!(collect(keys(species.energies)); by=l->species.energies[l])

"Returns the known decay transitions as lower => upper pairs, in energy order."
transitions(species::OneElectronSpecies) = sort!(
    [lo => hi for (lo, hi) in keys(species.einstein_as)];
    by=p->(species.energies[p.second], species.energies[p.first]),
)

function Base.show(io::IO, species::OneElectronSpecies)
    print(
        io,
        nameof(typeof(species)),
        "(",
        length(species.energies),
        " levels, mass=",
        species.mass,
    )
    species isa HyperfineOneElectronSpecies && print(io, ", I=", species.nuclear_spin)
    print(io, ")")
end

export Species,
    OneElectronSpecies,
    NoHyperfineOneElectronSpecies,
    HyperfineOneElectronSpecies,
    HyperfineConstants,
    LevelPolarisability,
    ImplicitPolarisability,
    ElectricDipole,
    MagneticDipole,
    einstein_a,
    lifetime,
    level_polarisability,
    saturation_intensity,
    transition_frequency,
    level_energy,
    levels,
    transitions
public BOHR_RADIUS, DIPOLE_AU, POLARISABILITY_AU
