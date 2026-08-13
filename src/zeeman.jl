using Unitful

"""
The Bohr magneton ``μ_B = e ħ / (2 m_e)``, from the physical constants built into
Unitful.
"""
const BOHR_MAGNETON = Unitful.q * u"ħ" / (2 * Unitful.me)

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
provides one (`lande_g_overrides` of a [`HyperfineOneElectronSpecies`](@ref)),
the LS-coupling expression of [`Levels.ls_lande_g`](@ref) otherwise. For a
hyperfine ``F`` level it is the **low-field** ``g_F``,

```math
g_F = g_J \\frac{F(F+1) + J(J+1) - I(I+1)}{2 F (F+1)}
    + g_I \\frac{F(F+1) + I(I+1) - J(J+1)}{2 F (F+1)},
```

with ``g_I`` in the Bohr-magneton convention of the species; note that at
finite field the hyperfine Zeeman effect is nonlinear (cf.
[`zeeman_shift`](@ref) and [`hyperfine_manifold`](@ref)), so ``g_F m_F`` is
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
    ls_lande_g(species.mass, spec.l, spec.j)
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
    g_j * c_j + species.nuclear_g * c_i
end

"""
    zeeman_shift(species, state::StateSpec, B)

Returns the first-order Zeeman shift of the given state in a static magnetic field
`B` along the quantisation axis, in angular frequency units.
"""
function zeeman_shift(species::NoHyperfineOneElectronSpecies, state::StateSpec, B)
    g = lande_g(species, state.level)
    uconvert(u"µs^-1", g * state.m * BOHR_MAGNETON * B / u"ħ")
end

"""
    zeeman_sensitivity(species, lower::StateSpec, upper::StateSpec)
    zeeman_sensitivity(species, transition::Pair)

Returns the first-order magnetic-field sensitivity ``χ`` of the transition
frequency between the two given states, in angular frequency units per magnetic
flux density (e.g. convertible to `µs^-1/mT`).

The `Pair` form accepts a `lower => upper` pair as produced by
[`state_pairs`](@ref).
"""
function zeeman_sensitivity(
    species::NoHyperfineOneElectronSpecies,
    lower::StateSpec,
    upper::StateSpec,
)
    g_lower = lande_g(species, lower.level)
    g_upper = lande_g(species, upper.level)
    uconvert(
        u"µs^-1/mT",
        (g_upper * upper.m - g_lower * lower.m) * BOHR_MAGNETON / u"ħ",
    )
end

zeeman_sensitivity(species::NoHyperfineOneElectronSpecies, transition::Pair) =
    zeeman_sensitivity(species, transition.first, transition.second)

"""
    zeeman_hamiltonian(species, basis::StateBasis, B)

Returns the matrix of the Zeeman Hamiltonian ``\\sum_J g_J (μ_B / ħ) \\vec{B} ⋅ \\vec{J}``
in the given basis, in angular frequency units.

`B` is the static Cartesian magnetic-field 3-vector, with the quantisation axis
along z; it must be real for the result to be Hermitian. The result is
block-diagonal in the levels of the basis, with matrix elements constructed
per-state, so any basis ordering (or subset of Zeeman sublevels) is supported.
"""
function zeeman_hamiltonian(
    species::NoHyperfineOneElectronSpecies,
    basis::StateBasis{NoHyperfineNumberSpec},
    B,
)
    # B·J in terms of ladder operators: B_z J_z + (B_x - iB_y)/2 J₊ + (B_x + iB_y)/2 J₋.
    b_z = complex(float(B[3]))
    b_plus = (B[1] + im * B[2]) / 2
    b_minus = (B[1] - im * B[2]) / 2

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
            u"µs^-1",
            lande_g(species, ket.level) * BOHR_MAGNETON * coupling / u"ħ",
        )
    end

    n = length(basis)
    [element(basis[i], basis[k]) for i in 1:n, k in 1:n]
end

export lande_g, zeeman_shift, zeeman_sensitivity, zeeman_hamiltonian
public BOHR_MAGNETON, ls_lande_g
