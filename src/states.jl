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
Specifies a particular hyperfine level of a species with nuclear spin using the
relevant low-field quantum numbers.

The nuclear spin ``I`` itself is a property of the species (cf.
[`HyperfineOneElectronSpecies`](@ref)), not the level, so it is not part of the
specification; consistency of `f` with ``|I - j| ≤ F ≤ I + j`` is only checked
by the functions that take the species. The total spin quantum number is
implicitly assumed to be ``1/2``, as for [`NoHyperfineNumberSpec`](@ref).
"""
struct HyperfineNumberSpec <: LevelSpec
    "Orbital angular momentum quantum number."
    l::Int

    "Electronic total angular momentum quantum number."
    j::Rational{Int}

    "Total (electronic plus nuclear) angular momentum quantum number."
    f::Rational{Int}
end

HyperfineNumberSpec(s::AbstractString) = convert(HyperfineNumberSpec, s)

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
Splits an optional trailing hyperfine `F=<value>` suffix (separated by
whitespace) off a level string, returning the fine-structure part and the `F`
value string (or `nothing` if there is no suffix).
"""
function split_hyperfine_suffix(val::AbstractString)
    m = match(r"^(.*?)\s+F\s*=\s*(\S+)\s*$", val)
    isnothing(m) ? (val, nothing) : (m.captures[1], m.captures[2])
end

"""
Parses an angular momentum quantum number given as a fraction (with optional
braces and Julia's `//` tolerated), or, if `allow_integer` is set, as a plain
integer. `orig` is the full level string for error messages.
"""
function parse_momentum(
    val::AbstractString,
    orig::AbstractString;
    allow_integer::Bool=false,
)
    if startswith(val, '{')
        if !endswith(val, '}')
            throw(
                ArgumentError(
                    "Invalid total angular momentum quantum number specification '$val': $orig",
                ),
            )
        end
        val = chopsuffix(chopprefix(val, "{"), "}")
    end

    fraction_idx = findfirst('/', val)
    if isnothing(fraction_idx)
        if allow_integer
            int = tryparse(Int, val)
            if !isnothing(int)
                return int // 1
            end
        end
        throw(
            ArgumentError(
                "Expected fractional total angular momentum quantum number, not '$val': $orig",
            ),
        )
    end
    if fraction_idx < 2
        throw(
            ArgumentError(
                "Expected fractional total angular momentum quantum number, not '$val': $orig",
            ),
        )
    end

    num_str = val[1:(fraction_idx-1)]
    num = tryparse(Int, num_str)
    if isnothing(num)
        throw(
            ArgumentError(
                "Invalid numerator '$num_str' for total angular momentum quantum number: $orig",
            ),
        )
    end

    # Skip the fraction slash, plus an optional second one (Julia's `//` notation).
    den_str = lstrip(val[(fraction_idx+1):end], '/')
    if isempty(den_str)
        throw(
            ArgumentError(
                "Expected denominator for total angular momentum quantum number: $orig",
            ),
        )
    end

    den = tryparse(Int, den_str)
    if isnothing(den)
        throw(
            ArgumentError(
                "Invalid denominator '$den_str' for total angular momentum quantum number: $orig",
            ),
        )
    end

    num // den
end

"""
Parses the fine-structure part of a level in spectroscopic notation (an
[`ORBITAL_SYMBOLS`](@ref) letter plus the total angular momentum fraction).
`orig` is the full level string for error messages.
"""
function parse_fine_structure(val::AbstractString, orig::AbstractString)
    if isempty(val)
        throw(ArgumentError("Empty level specification"))
    end

    l_sym = first(val)
    l_idx = findfirst(l_sym, ORBITAL_SYMBOLS)
    if isnothing(l_idx)
        throw(
            ArgumentError(
                "Unknown symbol '$l_sym' for orbital angular momentum quantum number: $orig",
            ),
        )
    end
    l = l_idx - 1

    NoHyperfineNumberSpec(l, parse_momentum(lstrip(val[2:end], '_'), orig))
end

"""
    convert(::Type{NoHyperfineNumberSpec}, s::SpectroscopicSpec)

Parses a level in spectroscopic notation into its quantum numbers.

The orbital angular momentum is given by one of [`ORBITAL_SYMBOLS`](@ref), followed
by the total angular momentum as a fraction. The underscore separating the two is
optional, as are braces around the fraction and the second slash of Julia's `//`
rational literal syntax; `"S_{1/2}"`, `"D_5/2"` and `"D5//2"` are all accepted.

An error is raised if the string carries a hyperfine `F=` suffix; such levels
parse to [`HyperfineNumberSpec`](@ref) instead.
"""
function Base.convert(::Type{NoHyperfineNumberSpec}, s::SpectroscopicSpec)
    val, f_str = split_hyperfine_suffix(s.val)
    if !isnothing(f_str)
        throw(
            ArgumentError(
                "Level '$(s.val)' specifies a hyperfine F quantum number; " *
                "use HyperfineNumberSpec (or drop the suffix)",
            ),
        )
    end
    parse_fine_structure(val, s.val)
end

"""
    convert(::Type{HyperfineNumberSpec}, s::SpectroscopicSpec)

Parses a hyperfine level in spectroscopic notation into its quantum numbers.

The fine-structure part follows the
[`NoHyperfineNumberSpec`](@ref) grammar, followed by whitespace and the total
angular momentum as `F=` with an integer or fractional value; `"S_1/2 F=4"`,
`"S_{1/2} F = 4"` and `"D_5/2 F=7//2"` are all accepted.
"""
function Base.convert(::Type{HyperfineNumberSpec}, s::SpectroscopicSpec)
    val, f_str = split_hyperfine_suffix(s.val)
    if isnothing(f_str)
        throw(
            ArgumentError(
                "Level '$(s.val)' is missing the hyperfine F quantum number " *
                "(e.g. \"S_1/2 F=4\")",
            ),
        )
    end
    fs = parse_fine_structure(val, s.val)
    HyperfineNumberSpec(fs.l, fs.j, parse_momentum(f_str, s.val; allow_integer=true))
end

"""
    parse_level(level)

Parses a level given in spectroscopic notation into its canonical number-spec
form — [`HyperfineNumberSpec`](@ref) if a hyperfine `F=` suffix is present,
[`NoHyperfineNumberSpec`](@ref) otherwise. Number specs pass through unchanged.
"""
parse_level(level::AbstractString) = parse_level(SpectroscopicSpec(level))

function parse_level(level::SpectroscopicSpec)
    _, f_str = split_hyperfine_suffix(level.val)
    convert(isnothing(f_str) ? NoHyperfineNumberSpec : HyperfineNumberSpec, level)
end

parse_level(level::LevelSpec) = level

"""
    fine_structure(level)

Returns the fine-structure [`NoHyperfineNumberSpec`](@ref) part of the given
level, projecting out the hyperfine `F` quantum number if present.

This projection is deliberately not a `convert` method: a conversion would let a
hyperfine state specification silently reinterpret its ``m_F`` as an ``m_J``
where a no-hyperfine one is expected, whereas an explicit call site documents
that only the fine-structure identity is meant.
"""
fine_structure(level) = fine_structure(parse_level(level))
fine_structure(spec::NoHyperfineNumberSpec) = spec
fine_structure(spec::HyperfineNumberSpec) = NoHyperfineNumberSpec(spec.l, spec.j)

"""
    momentum(level)

Returns the total angular momentum quantum number of the given level — the one
whose projection is the `m` of a [`StateSpec`](@ref): ``j`` for a
[`NoHyperfineNumberSpec`](@ref), ``F`` for a [`HyperfineNumberSpec`](@ref).
"""
momentum(level) = momentum(parse_level(level))
momentum(spec::NoHyperfineNumberSpec) = spec.j
momentum(spec::HyperfineNumberSpec) = spec.f

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
[`NoHyperfineNumberSpec`](@ref) form, or [`HyperfineNumberSpec`](@ref) when a
hyperfine `F=` suffix is present (in which case `m` is the projection ``m_F``).
"""
StateSpec(level::LevelSpec, m::Union{Integer,Rational}) =
    StateSpec{typeof(level)}(level, m)
StateSpec(level::AbstractString, m::Union{Integer,Rational}) =
    StateSpec(parse_level(level), m)

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

Both levels must be of the same kind — two fine-structure levels, or two
hyperfine `F` levels (for which the projections are ``m_F``).

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
    lower = parse_level(lower_level)
    upper = parse_level(upper_level)
    if typeof(lower) != typeof(upper)
        throw(
            ArgumentError(
                "Levels '$lower_level' and '$upper_level' must both be " *
                "fine-structure or both hyperfine levels",
            ),
        )
    end
    [
        StateSpec(lower, m_lower) => StateSpec(upper, m_upper) for
        m_lower in (-momentum(lower)):momentum(lower) for
        m_upper in (-momentum(upper)):momentum(upper) if m_upper - m_lower in Δm
    ]
end

export LevelSpec,
    NoHyperfineNumberSpec,
    HyperfineNumberSpec,
    StateSpec,
    SpectroscopicSpec,
    fine_structure,
    state_pairs
public ORBITAL_SYMBOLS, momentum, parse_level
