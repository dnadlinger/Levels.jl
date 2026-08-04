# [Zeeman structure](@id reference-zeeman)

Landé g-factors (LS coupling with the reduced-mass correction to the orbital
g-factor, or measured values where a species provides them), first-order Zeeman
shifts of individual states, magnetic-field sensitivities of transitions, and
the Zeeman Hamiltonian for an arbitrary field direction.

For hyperfine species the same functions are exact rather than first-order
quantities, evaluated from the manifold eigen-solutions of the
[hyperfine-structure machinery](@ref reference-hyperfine); `lande_g` of an
``F`` level returns the low-field ``g_F``.

Frequencies are angular throughout.

```@autodocs
Modules = [Levels]
Pages = ["zeeman.jl"]
```
