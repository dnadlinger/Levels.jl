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

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["quantum_toolbox.jl"]
```
