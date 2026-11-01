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
[Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) units.

Start with [Conventions and examples](05-conventions.md), especially the distinction
between angular frequencies and decay rates, then
[Bases and representations](06-bases.md) for the matrices behind the structure
and the construction of operators by hand. The API reference documents the core
calculations, `Levels.PeriodicDynamics` and `Levels.OpticalBloch`.

Julia 1.12 or newer is required. Install with `using Pkg; Pkg.add("Levels")` once
the package is registered; a development checkout can be loaded with
`Pkg.develop(path="/path/to/Levels.jl")`.
