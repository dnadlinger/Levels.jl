# QuantumToolbox.jl integration: the function stub (and user-facing docstring)
# for the Fourier decomposition of a DrivenTransition implemented in the
# LevelsQuantumToolboxExt package extension, which also provides the
# `QuantumObject`/`basis`/`projection` conversions over a `StateBasis` (those
# add methods to QuantumToolbox functions, so they need no stubs here).

"""
    fourier_hamiltonians(dt::DrivenTransition; δ = 0, time_unit = u"µs")

Returns the Fourier components of the rotating-frame Hamiltonian of the model
at laser detuning `δ`, as
[QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/) operators: a
`NamedTuple` `(; H_0, H_p, H_m, ωd)` with

```math
H(t) = H_0 + H_p e^{i ω_d t} + H_m e^{-i ω_d t},
```

where `H_0` is the static part (frame diagonal, detuning, and the unmodulated
laser coupling), `H_p = Σ_i X_i e^{i φ_i} / 2` collects the ``e^{i Ω t}``
halves of the harmonic drives, and `H_m = H_p^†`. The `QuantumObject` matrices
are dimensionless in the angular `time_unit`⁻¹ convention (µs⁻¹ by default,
matching the normalisation of the dynamics layer), with the drive frequency
`ωd` a plain `Float64` in the same units (so times passed to the
QuantumToolbox solvers are in `time_unit`). The field order is chosen such
that

```julia
steadystate_fourier(fourier_hamiltonians(dt)..., c_ops)
```

directly computes the periodic steady state of the driven master equation with
jump operators `c_ops`.

Laser phase modulation is not representable in a single-harmonic
decomposition, so the components correspond to
`modulation = HarmonicPhaseModulation()` (none) in
[`full_hamiltonian`](@ref).

!!! note
    This method is defined in the `LevelsQuantumToolboxExt` package extension,
    i.e. only once QuantumToolbox.jl is loaded.
"""
function fourier_hamiltonians end

export fourier_hamiltonians
