# [Zeeman structure](@id reference-zeeman)

Landé g-factors (LS coupling with the reduced-mass correction to the orbital g-factor,
or measured values where a species provides them), the zero-field magnetic-field
sensitivities of transitions, and the coupled-basis Zeeman Hamiltonian of a field vector
in an arbitrary direction — the building block of perturbing (ac) fields, not of the
static working field.

The static field lives in the [`Eigenstructure`](@ref): `zeeman_shift(s, state)` and
`zeeman_sensitivity(s, pair)` evaluate the labelled eigenstates at its field, exactly
for both species kinds (for a hyperfine species via the Hellmann–Feynman theorem on the
manifold eigenstates; `lande_g` of an ``F`` level returns the low-field ``g_F`` that the
species form uses).

Frequencies are angular throughout.

```@autodocs
Modules = [Levels]
Pages = ["zeeman.jl"]
```
