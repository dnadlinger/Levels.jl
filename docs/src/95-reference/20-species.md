# [Species and transitions](@id reference-species)

An atomic species is described by its level energies and the Einstein A
coefficients connecting them; transition frequencies, level lifetimes and
saturation intensities are derived from those two tables. A species with
nuclear spin additionally carries its hyperfine coupling constants, nuclear
g-factor and measured electronic g-factors, while the energy and rate tables
stay keyed on the fine-structure levels — hyperfine level arguments resolve to
those, with the zero-field hyperfine shifts added where the quantity is
frequency-like (cf. the [hyperfine-structure reference](@ref
reference-hyperfine)).

Frequencies and rates are in angular units throughout.

```@autodocs
Modules = [Levels]
Pages = ["species.jl"]
```
