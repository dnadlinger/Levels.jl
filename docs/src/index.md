```@meta
CurrentModule = Levels
```

# Levels

Levels computes atomic level structures, transition strengths, light shifts and
laser-driven dynamics for trapped-ion spectroscopy. It includes atomic data for ⁴⁰Ca⁺,
⁴³Ca⁺, ⁸⁸Sr⁺ and ¹³⁸Ba⁺, and accepts user-defined alkali-like species.

As a quick example, consider the ⁴³Ca⁺ 729 nm transition
``|S_{1/2}, F = 4, m_F = 4⟩ ↔ |D_{5/2}, F = 4, m_F = 3⟩``, which is first-order
insensitive to magnetic-field noise near 0.5 mT:

```@example introduction
using Levels, Unitful
clock = StateSpec("S_1/2 F=4", 4) => StateSpec("D_5/2 F=4", 3)
B₀ = insensitive_field(ca43, clock, (0.45u"mT", 0.55u"mT"))
```

At that field, let's calculate the Rabi frequency with which a 10 mW Gaussian beam
focused to a 20 µm waist drives it there, propagating perpendicular to the field and
polarised along it (for which hyperfine mixing is taken into account):

```@example introduction
s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=B₀)
direction, ε = beam_vectors(π / 2, 0)
Ω = rabi_frequency(s, clock, peak_intensity(10u"mW", 20u"µm"), ε, direction)
uconvert(u"kHz", cyclic(abs(Ω)))
```

The [`Eigenstructure`](@ref) is the atom: a species, a basis of states and the static
magnetic field along the quantisation axis, solved once. It is the only place the field
enters, and every model — light shifts, periodically driven transitions, optical Bloch
equations — is built on it. The API accepts spectroscopic strings for levels and
`StateSpec(level, m)` for magnetic substates. A transition is a `lower => upper` pair.
Physical quantities use [Unitful.jl](https://juliaphysics.github.io/Unitful.jl/) units.

Start with [Conventions](05-conventions.md) for units, states, basis transformations and
the construction of operators by hand. The API reference documents the core
calculations, `Levels.PeriodicDynamics` and `Levels.OpticalBloch`.

## Examples

- [Field-insensitive 729 nm transitions of ⁴³Ca⁺](@ref example-ca43-729): field
  sensitivity of the S₁/₂ ↔ D₅/₂ components, the first-order insensitive points of one
  of them and the second-order shift that remains there.

The pages are generated from the
[Literate.jl](https://fredrikekre.github.io/Literate.jl) scripts in `examples/`, which
you can also run yourself in the `docs` environment, e.g.
`julia --project=docs examples/ca43_729_field_sensitivity.jl`. You can also convert them
to Jupyter notebooks if you prefer (see Literate.jl documentation).
