# [Optical Bloch equations](@id reference-optical-bloch)

```@docs
Levels.OpticalBloch
```

## Specification

The problem is stated without reference to atomic data: which states, which
beams, which field, and — for the motional layer — which modes.

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["scheme.jl"]
```

## Rotating frame

Each fine-structure level carries one frame frequency, fixed along a spanning
tree of the beam graph so that every tree beam couples time-independently;
beams outside the tree leave beat notes.

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["frame.jl"]
```

## Lindblad model of the internal states

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["lindblad.jl"]
```

## Motional layer

The coupling of one motional mode to the internal dynamics, to first order in
the Lamb–Dicke parameters. Note the distinction between the *projected*
Lamb–Dicke factor of a beam (with the beam direction projected onto the mode)
and the *unprojected* Lamb–Dicke parameter of a decay transition (whose
emission direction is averaged into the recoil moment ``α``).

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["motion.jl"]
```

## Solvers

Steady states, populations, adiabatic-elimination cooling rates and the full
internal ⊗ motional model are provided through
[QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/) by the
`LevelsQuantumToolboxExt` package extension (see
[QuantumToolbox integration](@ref reference-quantum-toolbox)); the functions
below are defined there.

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["quantum_toolbox.jl"]
```
