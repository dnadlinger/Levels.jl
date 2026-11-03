# [QuantumToolbox integration](@id reference-quantum-toolbox)

Loading [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/) activates the
`LevelsQuantumToolboxExt` package extension, which maps Levels objects onto
`QuantumToolbox.QuantumObject`s so that its solvers (`mesolve`, `steadystate_fourier`,
…) apply directly:

- `QuantumObject(A, basis::StateBasis)` and `QuantumObject(A, s::Eigenstructure)` wrap
  an operator matrix (or a ket coefficient vector) over a [`StateBasis`](@ref) or the
  basis of an [`Eigenstructure`](@ref). Unitful elements of pure time dimension ``𝐓^p``
  are stripped in units of ``µs^p`` — µs⁻¹ for Hamiltonians in the angular convention,
  µs^(-1/2) for Lindblad jump operators — so times passed to the solvers are in µs; the
  reference time scale is configurable via the `time_unit` keyword.
- `basis(b, state)` and `projection(b, ket[, bra])`, for a `StateBasis` or an
  `Eigenstructure`, address kets and projectors through [`Levels.stateindex`](@ref).
- [`fourier_hamiltonians`](@ref Levels.PeriodicDynamics.fourier_hamiltonians) decomposes
  a [`DrivenTransition`](@ref Levels.PeriodicDynamics.DrivenTransition) into the
  Hamiltonian Fourier components consumed by `QuantumToolbox.steadystate_fourier`.

`examples/sr88-ac-zeeman-state-prep` is a worked optical-pumping steady-state
calculation on top of these conversions.

For internal and motional master equations, the extension adds `liouvillian(model)` and
`steadystate(model)`. A single-beat internal model uses `steadystate_fourier` and
returns its period-averaged density matrix. Population, cooling and time-evolution
methods are documented with their models on the
[Optical Bloch equations](@ref reference-optical-bloch) page.

```@autodocs
Modules = [Levels.PeriodicDynamics]
Pages = ["quantum_toolbox.jl"]
```
