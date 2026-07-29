"""
Specifies a particular electronic level within a species.

That is, all the quantum numbers to specify a state, except for the projection
of the total angular momentum on the quantization axis.
"""
abstract type LevelSpec end

"""
Specifies a particular electronic state using the relevant low-field
quantum numbers.

The total spin quantum number is implicitly assumed to be ``1/2`` (cf.
[`NoHyperfineOneElectronSpecies`](@ref)).
"""
struct NoHyperfineNumberSpec <: LevelSpec
    "Orbital angular momentum quantum number."
    l::Int

    "Total angular momentum quantum number."
    j::Rational{Int}
end

NoHyperfineNumberSpec(s::AbstractString) = convert(NoHyperfineNumberSpec, s)

"""
Specifies a particular electronic state using spectroscopic notation.
"""
struct SpectroscopicSpec <: LevelSpec
    val::String
end

"""
Orbital angular momentum symbols in spectroscopic notation, in order of
increasing ``l`` (which is thus the zero-based index into this string).
"""
const ORBITAL_SYMBOLS = "SPDFG"

Base.convert(::Type{T}, s::AbstractString) where {T<:LevelSpec} =
    convert(T, SpectroscopicSpec(s))

"""
    convert(::Type{NoHyperfineNumberSpec}, s::SpectroscopicSpec)

Parses a level in spectroscopic notation into its quantum numbers.

The orbital angular momentum is given by one of [`ORBITAL_SYMBOLS`](@ref), followed
by the total angular momentum as a fraction. The underscore separating the two is
optional, as are braces around the fraction and the second slash of Julia's `//`
rational literal syntax; `"S_{1/2}"`, `"D_5/2"` and `"D5//2"` are all accepted.
"""
function Base.convert(::Type{NoHyperfineNumberSpec}, s::SpectroscopicSpec)
    val = s.val
    if isempty(val)
        throw(ArgumentError("Empty level specification"))
    end

    l_sym = first(val)
    l_idx = findfirst(l_sym, ORBITAL_SYMBOLS)
    if isnothing(l_idx)
        throw(
            ArgumentError(
                "Unknown symbol '$l_sym' for orbital angular momentum quantum number: $(s.val)",
            ),
        )
    end
    l = l_idx - 1

    val = lstrip(val[2:end], '_')

    if startswith(val, '{')
        if !endswith(val, '}')
            throw(
                ArgumentError(
                    "Invalid total angular momentum quantum number specification '$val': $(s.val)",
                ),
            )
        end
        val = chopsuffix(chopprefix(val, "{"), "}")
    end

    fraction_idx = findfirst('/', val)
    if isnothing(fraction_idx) || fraction_idx < 2
        throw(
            ArgumentError(
                "Expected fractional total angular momentum quantum number, not '$val': $(s.val)",
            ),
        )
    end

    j_num_str = val[1:(fraction_idx-1)]
    j_num = tryparse(Int, j_num_str)
    if isnothing(j_num)
        throw(
            ArgumentError(
                "Invalid numerator '$j_num_str' for total angular momentum quantum number: $(s.val)",
            ),
        )
    end

    # Skip the fraction slash, plus an optional second one (Julia's `//` notation).
    j_den_str = lstrip(val[(fraction_idx+1):end], '/')
    if isempty(j_den_str)
        throw(
            ArgumentError(
                "Expected denominator for total angular momentum quantum number: $(s.val)",
            ),
        )
    end

    j_den = tryparse(Int, j_den_str)
    if isnothing(j_den)
        throw(
            ArgumentError(
                "Invalid denominator '$j_den_str' for total angular momentum quantum number: $(s.val)",
            ),
        )
    end

    NoHyperfineNumberSpec(l, j_num // j_den)
end

"""
Specifies a particular electronic state, i.e. a level plus the projection of the
total angular momentum on the quantisation axis.
"""
struct StateSpec{L<:LevelSpec}
    level::L

    """
    Total momentum projection on the quantisation (z) axis (usually denoted
    ``m_J`` or ``m_F`` depending on the species).
    """
    m::Rational{Int}
end

"""
    StateSpec(level, m)

Creates a state specification from any level specification and the total angular
momentum projection `m`.

Strings are parsed as spectroscopic notation into the canonical
[`NoHyperfineNumberSpec`](@ref) form.
"""
StateSpec(level::LevelSpec, m::Union{Integer,Rational}) =
    StateSpec{typeof(level)}(level, m)
StateSpec(level::AbstractString, m::Union{Integer,Rational}) =
    StateSpec(convert(NoHyperfineNumberSpec, level), m)

"""
    convert(::Type{StateSpec{L}}, s::StateSpec)

Converts the level part of a state specification, e.g. into its canonical
`StateSpec{NoHyperfineNumberSpec}` form from one given in spectroscopic notation.
"""
Base.convert(::Type{StateSpec{L}}, s::StateSpec) where {L<:LevelSpec} =
    StateSpec(convert(L, s.level), s.m)

"""
    state_pairs(lower_level, upper_level; Δm)

Enumerates the transitions between the Zeeman states of the two given levels, as
`lower => upper` pairs of [`StateSpec`](@ref)s.

`Δm` is a collection giving the allowed values of the signed projection difference
`upper.m - lower.m`. There is no default, as the appropriate choice depends on the
multipole order of the transition and the field geometry; e.g. `Δm = -1:1` covers
electric-dipole transitions, and `Δm = [-2, -1, 1, 2]` the eight ⁸⁸Sr⁺
S``_{1/2}`` ↔ D``_{5/2}`` components observed with the ``Δm = 0`` quadrupole
amplitude suppressed by geometry.

The pairs are ordered by increasing lower-state `m`, then upper-state `m`; re-sort
as needed, e.g. by magnetic-field sensitivity using [`zeeman_sensitivity`](@ref).
"""
function state_pairs(lower_level, upper_level; Δm)
    lower = convert(NoHyperfineNumberSpec, lower_level)
    upper = convert(NoHyperfineNumberSpec, upper_level)
    [
        StateSpec(lower, m_lower) => StateSpec(upper, m_upper) for
        m_lower in (-lower.j):(lower.j) for
        m_upper in (-upper.j):(upper.j) if m_upper - m_lower in Δm
    ]
end

export LevelSpec, NoHyperfineNumberSpec, StateSpec, SpectroscopicSpec, state_pairs
public ORBITAL_SYMBOLS
