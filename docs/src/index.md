```@meta
CurrentModule = Levels
```

# Levels

Levels computes atomic level structures, transition strengths, light shifts and
laser-driven dynamics for trapped-ion spectroscopy. It includes atomic data for
⁸⁸Sr⁺, ⁴⁰Ca⁺ and ⁴³Ca⁺, and accepts user-defined alkali-like species.

```@example introduction
using Levels, Unitful
s = Eigenstructure(ca40, "S_1/2", "D_5/2"; B=0.5u"mT")
clock = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -1//2)
cyclic(transition_frequency(s, clock)) |> x -> uconvert(u"THz", x)
```

The [`Eigenstructure`](@ref) is the atom: a species, a basis of states and the
static magnetic field along the quantisation axis, solved once. It is the only
place the field enters, and every model — light shifts, periodically driven
transitions, optical Bloch equations — is built on it. The API accepts
spectroscopic strings for levels and `StateSpec(level, m)` for magnetic
substates. A transition is a `lower => upper` pair. Physical quantities use
[Unitful.jl](https://juliaphysics.github.io/Unitful.jl/) units.

Start with [Conventions and examples](05-conventions.md), especially the distinction
between angular frequencies and decay rates, then
[Bases and representations](06-bases.md) for the matrices behind the structure
and the construction of operators by hand. The API reference documents the core
calculations, `Levels.PeriodicDynamics` and `Levels.OpticalBloch`.

## Examples

- [Field-insensitive 729 nm transitions of ⁴³Ca⁺](@ref example-ca43-729): field
  sensitivity of the S₁/₂ ↔ D₅/₂ components, the first-order insensitive points of
  one of them and the second-order shift that remains there.

The pages are generated from the [Literate.jl](https://fredrikekre.github.io/Literate.jl)
scripts in `examples/`, which you can also run yourself in the `docs` environment, e.g.
`julia --project=docs examples/ca43_729_field_sensitivity.jl`. You can also convert them
to Jupyter notebooks if you prefer (see Literate.jl documentation).
