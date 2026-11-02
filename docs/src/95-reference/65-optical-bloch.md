# [Optical Bloch equations](@id reference-optical-bloch)

```@meta
CurrentModule = Levels.OpticalBloch
```

```@docs
Levels.OpticalBloch
```

Loading [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/) activates
`LevelsQuantumToolboxExt`, which provides `steadystate(model)` and
`liouvillian(model)` for internal and motional models, as well as
[`motional_model`](@ref), [`populations`](@ref), [`mean_phonon_number`](@ref)
and [`cooling_metrics`](@ref). See also
[QuantumToolbox integration](@ref reference-quantum-toolbox).

## Internal-state model

### Lasers and rotating frame

```@docs
LaserBeam
LaserScheme
rotating_frame
RotatingFrame
frame_offset
```

### Master equation and decay

```@docs
lindblad_model
LindbladModel
decay_operators
dephasing_operators
DecayLabel
DephasingLabel
```

## Coupling to motion

### Mode geometry and recoil

```@docs
MotionalMode
motional_coupling
MotionalCoupling
lamb_dicke
recoil_moment
```

### Joint internal and motional model

The joint model and all three cooling methods require an internal model
without beat harmonics.

```@docs
motional_model
MotionalModel
RecoilLabel
HeatingLabel
```

## States and observables

`steadystate(model)` returns the internal or joint density matrix (the
period average for an internal model with a single beat frequency).
[`populations`](@ref) takes an internal state, obtained from a joint state
with `QuantumToolbox.ptrace(ρ, 1)`. [`mean_phonon_number`](@ref) and
[`fock_populations`](@ref) take the joint state directly.

```@docs
populations
level_projector
mean_phonon_number
fock_populations
phonon_number_operator
```

### Thermal initial states

The thermal helpers truncate and renormalise a geometric distribution.
Their `nbar` parameter is the mean of the untruncated distribution; its
truncated direct expectation can be smaller. [`fock_truncation`](@ref)
chooses a ladder size from the discarded thermal tail.

```@docs
thermal_populations
thermal_state
fock_truncation
```

## Cooling estimates and dynamics

[`cooling_metrics`](@ref) is the common entry point for the steady occupation
and a relaxation time per mode. The methods attach different meanings to
that time:

| Method | Motional treatment | Meaning of `τ_c` |
| :----- | :----------------- | :--------------- |
| [`AdiabaticElimination`](@ref) | Independent rate equations | Exponential occupation relaxation |
| [`IntegratedTransient`](@ref) | Joint internal ⊗ Fock model, with optional coherence restriction | Normalised transient area from a specified thermal initial state |
| [`LiouvillianSpectrum`](@ref) | Joint internal ⊗ Fock model, with optional coherence restriction | Time of a selected eigenmode with appreciable phonon-number weight |

The times agree for a single exponential in the rate-equation regime.
Outside that regime, the integral time depends on the initial state and
the spectral estimate depends on the mode selection. [`cooling_curve`](@ref)
resolves the transient directly.

```@docs
cooling_metrics
CoolingMetrics
CoolingMethod
```

### Adiabatic elimination

```@docs
AdiabaticElimination
```

### Integrated transients

```@docs
IntegratedTransient
```

### Liouvillian spectrum

```@docs
LiouvillianSpectrum
```

### Time evolution

```@docs
cooling_curve
```

## [Numerical methods](@id optical-bloch-numerical-methods)

The plain-matrix forms work without QuantumToolbox. Hamiltonians use inverse
`time_unit`, jumps inverse square-root `time_unit`. The tensor order is
internal states first, last Fock mode fastest. Model-based constructors
require the extension and supply the layout and units.

### Displacements and emission quadrature

```@docs
displacement_elements
joint_displacement
emission_rule
```

### Liouvillian construction and representation

```@docs
BandLiouvillian
band_vector
full_matrix
expectation_row
```

### Steady states and transient integrals

[`IntegratedTransientSolver`](@ref) supports arbitrary observables and
initial states, including internal-state problems without Fock factors.

```@docs
IntegratedTransientSolver
steady_density_matrix
steady_value
integral_relaxation_time
```

### Relaxation eigenmodes

```@docs
relaxation_modes
```

### Internal construction helpers

```@docs
beam_levels
fine_structure_levels
beam_centroid_offset
validate_scheme
split_beam_coupling
merge_harmonics
validate_lamb_dicke_order
full_dimension
```
