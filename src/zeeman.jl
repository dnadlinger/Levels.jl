using Unitful

"""
The Bohr magneton ``μ_B = e ħ / (2 m_e)``, from the physical constants built into
Unitful.
"""
const BOHR_MAGNETON = Unitful.q * u"ħ" / (2 * Unitful.me)

"The nuclear magneton ``μ_N = e ħ / (2 m_p)``."
const NUCLEAR_MAGNETON = Unitful.q * u"ħ" / (2 * Unitful.mp)

# Magnitude of the free-electron g-factor, 2.00231930436092(36) [CODATA2022]
# (E. Tiesinga et al., Rev. Mod. Phys. 97, 025002 (2025), doi:10.1103/RevModPhys.97.025002).
const FREE_ELECTRON_G = 2.00231930436092

"""
Returns the LS-coupling Landé g-factor for a single electron (``s = 1/2``),
using the free-electron spin g-factor (CODATA 2022) and the leading-order
reduced-mass correction to the orbital g-factor, ``g_L = 1 - m_e/M`` with ``M``
the ion mass (any distinction from the nuclear mass is far beyond the accuracy
of this leading-order form).
"""
function ls_lande_g(mass, l, j)
    s = 1//2

    g_l = 1 - uconvert(NoUnits, Unitful.me / mass)
    c_l = (j * (j + 1) - s * (s + 1) + l * (l + 1)) // (2 * j * (j + 1))
    c_s = (j * (j + 1) + s * (s + 1) - l * (l + 1)) // (2 * j * (j + 1))
    g_l * c_l + FREE_ELECTRON_G * c_s
end

"""
    lande_g(species, level)

Returns the Landé g-factor of the given level.

For a fine-structure level this is ``g_J`` — a measured value where the species
provides one (its `lande_g_overrides`), the LS-coupling expression of
[`Levels.ls_lande_g`](@ref) otherwise. For a
hyperfine ``F`` level it is the **low-field** ``g_F``,

```math
g_F = g_J \\frac{F(F+1) + J(J+1) - I(I+1)}{2 F (F+1)}
    - g_I (μ_N/μ_B) \\frac{F(F+1) + I(I+1) - J(J+1)}{2 F (F+1)},
```

with ``g_I`` in the nuclear-magneton convention of the species; note that at
finite field the hyperfine Zeeman effect is nonlinear (cf.
[`zeeman_shift`](@ref) on an [`Eigenstructure`](@ref)), so ``g_F m_F`` is
only the ``B → 0`` slope.
"""
function lande_g(species::NoHyperfineOneElectronSpecies, level)
    spec = parse_level(level)
    if !(spec isa NoHyperfineNumberSpec)
        throw(
            ArgumentError(
                "Level '$level' carries hyperfine structure, but the species " *
                "has none",
            ),
        )
    end
    get(species.lande_g_overrides, spec) do
        ls_lande_g(species.mass, spec.l, spec.j)
    end
end

function lande_g(species::HyperfineOneElectronSpecies, level)
    spec = parse_level(level)
    fs = fine_structure(spec)
    g_j = get(species.lande_g_overrides, fs) do
        ls_lande_g(species.mass, fs.l, fs.j)
    end
    spec isa NoHyperfineNumberSpec && return g_j

    validate_hyperfine(species, spec)
    i = species.nuclear_spin
    j = spec.j
    f = spec.f
    iszero(f) && return zero(g_j)
    c_j = (f * (f + 1) + j * (j + 1) - i * (i + 1)) / (2 * f * (f + 1))
    c_i = (f * (f + 1) + i * (i + 1) - j * (j + 1)) / (2 * f * (f + 1))
    g_j * c_j -
    species.nuclear_g * uconvert(NoUnits, NUCLEAR_MAGNETON / BOHR_MAGNETON) * c_i
end

"""
    zeeman_sensitivity(species, lower::StateSpec => upper::StateSpec)

Returns the zero-field (first-order) magnetic-field sensitivity ``χ`` of the
transition between the two given states, in angular frequency units per
magnetic flux density (e.g. convertible to `rad/µs/mT`): the slope
``(g_u m_u - g_l m_l) μ_B/ħ`` with the Landé factors of [`lande_g`](@ref) —
``g_J`` for fine-structure states, the low-field ``g_F`` for hyperfine ones.
For a fine-structure species this is exact at any field; for a hyperfine
species the sensitivity at a working field is that of the
[`Eigenstructure`](@ref) form.

The pair may be one as produced by [`state_pairs`](@ref).
"""
function zeeman_sensitivity(
    species::OneElectronSpecies,
    transition::Pair{<:StateSpec,<:StateSpec},
)
    lower = validate_state_kind(species, transition.first)
    upper = validate_state_kind(species, transition.second)
    g_lower = lande_g(species, lower.level)
    g_upper = lande_g(species, upper.level)
    uconvert(
        u"rad/µs/mT",
        (g_upper * upper.m - g_lower * lower.m) * BOHR_MAGNETON / u"ħ",
    )
end

"""
    zeeman_hamiltonian(species, basis::StateBasis, b)

Returns the matrix of the Zeeman Hamiltonian of a magnetic field with the
Cartesian 3-vector `b` — ``\\sum_J g_J (μ_B / ħ) \\vec{b} ⋅ \\vec{J}`` for a
fine-structure species, ``\\vec{b} ⋅ (g_J μ_B \\vec{J} - g_I μ_N \\vec{I}) / ħ``
for a hyperfine one — over the given basis in the zero-field (coupled) basis
representation, in angular frequency units.

`b` is a field vector in the laboratory frame with the quantisation axis along
ẑ, and must be real for the result to be Hermitian. It is the building block
of perturbing fields (the ac drives of `Levels.PeriodicDynamics.zeeman_drives`),
not the static working field, which belongs to an [`Eigenstructure`](@ref);
[`zeeman_hamiltonian(s, b)`](@ref zeeman_hamiltonian(::Eigenstructure, ::Any))
gives the same operator in a structure's field eigenbasis. The result is
block-diagonal in the fine-structure levels, with matrix elements constructed
per state, so any basis ordering or subset of sublevels is supported; for a
hyperfine species it is not diagonal in ``F`` (the field mixes ``F`` within
each manifold while staying diagonal in ``m_F`` for a field along ẑ).
"""
function zeeman_hamiltonian(
    species::NoHyperfineOneElectronSpecies,
    basis::StateBasis{NoHyperfineNumberSpec},
    b,
)
    length(b) == 3 || throw(ArgumentError("The field must be a Cartesian 3-vector"))
    # b·J in terms of ladder operators: b_z J_z + (b_x - ib_y)/2 J₊ + (b_x + ib_y)/2 J₋.
    b_z = complex(float(b[3]))
    b_plus = (b[1] + im * b[2]) / 2
    b_minus = (b[1] - im * b[2]) / 2

    function element(bra, ket)
        j = ket.level.j
        m = ket.m
        coupling = if bra.level != ket.level
            zero(b_z)
        elseif bra.m == m
            float(m) * b_z
        elseif bra.m == m + 1
            sqrt(float(j * (j + 1) - m * (m + 1))) * b_minus
        elseif bra.m == m - 1
            sqrt(float(j * (j + 1) - m * (m - 1))) * b_plus
        else
            zero(b_z)
        end
        uconvert(
            u"rad/µs",
            lande_g(species, ket.level) * BOHR_MAGNETON * coupling / u"ħ",
        )
    end

    num_states = length(basis)
    [element(basis[i], basis[k]) for i in 1:num_states, k in 1:num_states]
end

export lande_g, zeeman_shift, zeeman_sensitivity, zeeman_hamiltonian
public BOHR_MAGNETON, NUCLEAR_MAGNETON, ls_lande_g
