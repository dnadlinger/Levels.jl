# [Periodic dynamics](@id reference-periodic-dynamics)

```@docs
Levels.PeriodicDynamics
```

## Driven-transition model

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["driving.jl"]
```

## Method selection

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["methods.jl"]
```

## Dressed-state engine

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["floquet.jl"]
```

## Quasienergy and ac Zeeman shifts

The time-averaged shifts an ac magnetic drive imparts on the states themselves,
via second-order Floquet perturbation theory (fast, with per-intermediate-state
diagnostics and near-resonance warnings) or the nonperturbative Floquet
dressing as its cross-check. For hyperfine species, these operate on the exact
[`Levels.hyperfine_manifold`](@ref) eigenstates.

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["quasienergy.jl"]
```

## Exact engine

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["tdse.jl"]
```
