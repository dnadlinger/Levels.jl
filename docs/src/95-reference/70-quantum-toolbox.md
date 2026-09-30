# [QuantumToolbox integration](@id reference-quantum-toolbox)

Loading [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/) activates
the `LevelsQuantumToolboxExt` package extension, which maps Levels objects
onto `QuantumToolbox.QuantumObject`s so that its solvers (`mesolve`,
`steadystate_fourier`, …) apply directly:

- `QuantumObject(A, basis::StateBasis)` wraps an operator matrix (or a ket
  coefficient vector) over a [`StateBasis`](@ref). Unitful elements of pure
  time dimension ``𝐓^p`` are stripped in units of ``µs^p`` — µs⁻¹ for
  Hamiltonians in the angular convention, µs^(-1/2) for Lindblad jump
  operators — so times passed to the solvers are in µs; the reference time
  scale is configurable via the `time_unit` keyword.
- `basis(b::StateBasis, state)` and `projection(b::StateBasis, ket[, bra])`
  address kets and projectors through [`Levels.stateindex`](@ref).
- [`fourier_hamiltonians`](@ref Levels.PeriodicDynamics.fourier_hamiltonians)
  decomposes a [`DrivenTransition`](@ref Levels.PeriodicDynamics.DrivenTransition)
  into the Hamiltonian Fourier components consumed by
  `QuantumToolbox.steadystate_fourier`.

`examples/sr88-ac-zeeman-state-prep` is a worked optical-pumping steady-state
calculation on top of these conversions.

For the [`Levels.OpticalBloch`](@ref) master-equation models the extension
adds solver methods:

- `liouvillian(model::LindbladModel; time_unit)` and
  `steadystate(model::LindbladModel; …)` — the internal-state Liouvillian and
  steady state (a model with a single beat-note harmonic is solved for its
  periodic steady state with `steadystate_fourier`, of which the
  period-averaged component is returned);
- [`populations`](@ref Levels.OpticalBloch.populations),
  [`cooling_rates`](@ref Levels.OpticalBloch.cooling_rates) (adiabatic
  elimination of the internal dynamics for one or several modes),
  [`motional_model`](@ref Levels.OpticalBloch.motional_model) (the full
  internal ⊗ motional operators),
  [`mean_phonon_number`](@ref Levels.OpticalBloch.mean_phonon_number) and
  [`cooling_time`](@ref Levels.OpticalBloch.cooling_time), documented on the
  [Optical Bloch equations](@ref reference-optical-bloch) page.

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["quantum_toolbox.jl"]
```
