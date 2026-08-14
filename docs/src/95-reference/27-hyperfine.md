# [Hyperfine structure](@id reference-hyperfine)

Species with nuclear spin ([`HyperfineOneElectronSpecies`](@ref), e.g.
[`ca43`](@ref)) resolve each fine-structure level into its hyperfine ``F``
levels ([`HyperfineNumberSpec`](@ref), `"S_1/2 F=4"` in spectroscopic
notation), whose ``2F + 1`` sublevels carry the projection ``m_F``.

The canonical matrix/index convention is the coupled ``|F, m_F⟩`` basis (e.g.
`StateBasis(species, "S_1/2", "D_5/2")` expands the fine-structure manifolds
into all their ``F`` levels): the hyperfine interaction is diagonal there
(the Casimir shifts of [`hyperfine_shift`](@ref)), and every matrix is directly
indexable through the [`StateBasis`](@ref) machinery. The ``|m_J, m_I⟩``
product basis appears only as an internal construction device — operators that
are natural as Kronecker products (the magnetic-moment components of
[`Levels.moment_operators`](@ref), electronic slow / nuclear fast, each in
order of increasing projection) are conjugated with the Clebsch–Gordan unitary
of [`Levels.coupling_transform`](@ref) (``⟨J m_J; I m_I | F m_F⟩``,
Condon–Shortley phases, electron first). The kron storage order deliberately
mirrors the electron-first coupling order, but the two are independent
choices: the coupling order fixes the signs of the coupled states, while the
storage order is arbitrary bookkeeping with no sign consequences.

Three independent convention layers fix the signs of hyperfine amplitudes.
The Condon–Shortley phases of the spherical harmonics — and with them of the
Clebsch–Gordan coefficients, all real with the stretched components positive
— are shared by essentially all modern references and assumed throughout.
The form of the Wigner–Eckart theorem is a non-choice: for integer rank,
Racah's 3-j form defines exactly the same reduced matrix element as the
Clebsch–Gordan form used here. The one genuine choice is the
angular-momentum coupling order of the ``|F, m_F⟩`` states: the package
couples electron first, ``|(J I) F⟩``, as in most of the atomic-physics
literature. Sources coupling nuclear spin first, ``|(I J) F⟩``, have states
differing by ``(-1)^{I + J - F}``, so their multipole amplitudes between
``F`` levels differ by the per-multiplet sign ``(-1)^{(F' - F) - (J' - J)}``
(magnitudes and all physical predictions agree). See
[`Levels.hyperfine_reduction`](@ref) for details and the literature
cross-check.

The Zeeman convention is ``H_Z = μ_B B (g_J m_J + g_I m_I)``, i.e. the nuclear
g-factor is stated in Bohr magnetons with the electron-like sign absorbed
(``μ_I = -g_I I μ_B``), and it is the *effective* moment of the nucleus bound
in the ion, not corrected for diamagnetic shielding.

At working fields the Zeeman interaction competes with the hyperfine intervals
(the Breit–Rabi regime — for the ⁴³Ca⁺ D``_{5/2}`` manifold already at a few
gauss), so ``F`` is mixed and only ``m_F`` remains exact.
[`hyperfine_manifold`](@ref) diagonalises one manifold exactly and labels the
eigenstates **adiabatically**: within each exact-``m_F`` block, ``F`` labels
follow the energy order of the zero-field levels, which is rigorous because
levels of equal ``m_F`` do not cross as a function of the field. The labels
are therefore field-dependent bookkeeping for eigenstates, undefined at
``B = 0``, and [`zeeman_shift`](@ref) / [`zeeman_sensitivity`](@ref) for
hyperfine states are exact quantities derived from the eigen-solution (the
latter via the Hellmann–Feynman theorem), not first-order expressions.
[`insensitive_field`](@ref) finds the roots of the sensitivity, e.g. the
⁴³Ca⁺ clock-transition points at 3.38 G and 4.96 G.
[`eigenbasis_transform`](@ref) exposes the rotation from the coupled basis
into the labelled field eigenbasis as a matrix, so any coupled-basis operator
— a [`quadrupole_couplings`](@ref) matrix, a
[`zeeman_hamiltonian`](@ref) of an rf field — follows the eigenstates as
``V^† X V``; per-component transition amplitudes are available directly via
[`transition_amplitude`](@ref).

```@autodocs
Modules = [Levels]
Pages = ["hyperfine.jl"]
```
