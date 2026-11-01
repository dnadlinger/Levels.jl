# [Optical Bloch equations](@id reference-optical-bloch)

```@docs
Levels.OpticalBloch
```

## Specification

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["scheme.jl"]
```

## Rotating frame

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

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["motion.jl"]
```

## Numerical methods

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["liouvillian.jl"]
```

## Solvers

Steady states, populations, adiabatic-elimination cooling rates, the full
internal ⊗ motional model to any order in the Lamb–Dicke parameters and the
common [`cooling_metrics`](@ref Levels.OpticalBloch.cooling_metrics) front-end
of the cooling estimates are provided through
[QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/) by the
`LevelsQuantumToolboxExt` package extension (see
[QuantumToolbox integration](@ref reference-quantum-toolbox)); the functions
below are defined there, the types here.

```@autodocs
Modules = [Levels.OpticalBloch]
Pages = ["quantum_toolbox.jl"]
```
