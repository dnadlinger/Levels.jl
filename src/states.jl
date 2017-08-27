using Match

import Base: convert

"""
Specifies a particular electronic level within a species.

That is, all the quantum numbers to specify a state, except for the projection
of the total angular momentum on the quantization axis.
"""
abstract type LevelSpec end

"""
Specifies a particular electronic state using the relevant low-field
quantum numbers.

The total spin quantum number is implicitly assumed to be :math:`1 / 2` (cf.
`NoHyperfineOneElectronSpecies`).
"""
immutable NoHyperfineNumberSpec <: LevelSpec
    "Orbital angular momentum quantum number."
    l::Integer

    "Total angular momentum quantum number."
    j::Rational
end

"""
Specifies a particular electronic state using spectroscopic notation.
"""
immutable SpectroscopicSpec <: LevelSpec
    val::String
end

convert{T<:LevelSpec}(::Type{T}, s::String) = convert(T, SpectroscopicSpec(s))

function convert(::Type{NoHyperfineNumberSpec}, s::SpectroscopicSpec)
    val = s.val

    l_sym = val[1]
    l = findfirst(['S', 'P', 'D', 'F', 'G'], l_sym) - 1
    if l == -1
        throw(ArgumentError("Unknown symbol '$(s[1])' for orbital angular momentum quantum number: $(s.val)"))
    end

    val = lstrip(val[2 : end], '_')

    if val[1] == '{'
        if val[end] != '}'
            throw(ArgumentError("Invalid total angular momentum quantum number specification '$val': $(s.val)"))
        end
        val = val[2 : end - 1]
    end

    # Parse numerator.
    fraction_idx = findfirst(val, '/')
    if fraction_idx < 2
        throw(ArgumentError("Expected fractional total angular momentum quantum number, not '$val': $(s.val)"))
    end

    j_num_str = val[1 : fraction_idx - 1]
    j_num = try
        parse(Int64, j_num_str)
    catch
        throw(ArgumentError("Invalid numerator '$j_num_str' for total angular momentum quantum number: $(s.val)"))
    end

    # Skip fraction slash.
    val = val[fraction_idx + 1 : end]

    if length(val) == 0
        throw(ArgumentError("Expected denominator for total angular momentum quantum number: $(s.val)"))
    end

    # Skip optional second slash (Julia's fraction notation).
    val = lstrip(val, '/')

    # Parse denominator.
    j_den = try
        parse(Int64, val)
    catch
        throw(ArgumentError("Invalid denominator '$val' for total angular momentum quantum number: $(s.val)"))
    end

    NoHyperfineNumberSpec(l, j_num // j_den)
end

immutable StateSpec
    level::LevelSpec

    """
    Total momentum projection on the quantisation (z) axis (usually denoted
    :math:`m_J` or :math:`m_F` depending on the species).
    """
    m::Rational
end

export LevelSpec, NoHyperfineNumberSpec, convert, StateSpec, SpectroscopicSpec