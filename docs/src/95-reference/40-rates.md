# [Rabi frequencies and coupling strengths](@id reference-rates)

Angular momentum coupling factors between individual Zeeman states, and the
resulting Rabi frequencies for a given drive intensity.

For hyperfine states, the amplitudes relative to the fine-structure reduced
matrix element factorise at zero field into an ``F``-basis Clebsch–Gordan
coefficient and the 6-j reduction factor of
[`Levels.hyperfine_reduction`](@ref); the Einstein A coefficients remain
fine-structure quantities, of which the squared reduction factors are the
branching fractions. At a working field the Zeeman interaction mixes ``F``
within each manifold, so [`transition_amplitude`](@ref) (and
[`rabi_frequency`](@ref)) also take the static flux density as a final
argument, evaluating the amplitude exactly between the adiabatically-labelled
eigenstates of [`hyperfine_manifold`](@ref); whole coupling matrices rotate
into the field eigenbasis with [`eigenbasis_transform`](@ref).

```@autodocs
Modules = [Levels]
Pages = ["rates.jl"]
```
