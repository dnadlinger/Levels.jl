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

### Estimator comparison

Adiabatic elimination replaces the joint dynamics by a motional birth–death
process with cooling and heating rates ``nA_-`` and ``(n+1)A_+``. Without
background heating, the mean occupation obeys

```math
\dot{\bar n} = -W(\bar n - \bar n_∞), \qquad
W = A_- - A_+, \qquad
\bar n_∞ = \frac{A_+}{A_- - A_+},
```

so all three methods give the time constant ``1/W`` in this regime
([Eschner et al. (2003)](https://doi.org/10.1364/JOSAB.20.001003),
Eqs. (6)–(8)).

The reduction requires Lamb–Dicke confinement and internal relaxation much
faster than motional evolution. The sideband coupling must remain weak enough
to avoid saturating the internal resonance it addresses. In EIT cooling, the
narrow dressed-state linewidth can set this scale, so ``ηΩ/ω_m`` alone is not a
general validity test ([Morigi (2003)](https://doi.org/10.1103/PhysRevA.67.033402),
discussion after Eq. (33) and comparison with quantum Monte Carlo in Fig. 3).
Sideband matrix elements grow as ``\sqrt{n}``, making the initial occupation
relevant to validity. The Lamb–Dicke expansion and adiabatic elimination are
separate approximations: a first-order joint model (`lamb_dicke_order = 1`)
still retains coherent, saturated sideband dynamics.

Outside the rate-equation regime, cooling can be non-exponential and the
motional state need not remain thermal.
[Scharnhorst et al. (2018)](https://doi.org/10.1103/PhysRevA.98.023424)
demonstrate this for EIT cooling: strong spin–motion coupling disrupts the
dark state and slows the initial cooling (Sec. IV.B.2, Fig. 7). Their comparison
of analytical rates, master-equation simulations and experiment illustrates
the breakdown of the reduced theory (Sec. V.4, Fig. 13).

The two joint-model estimates describe different aspects of this relaxation.
For a diagonalisable generator, a decaying occupation transient has the form

```math
\bar n(t) - \bar n_∞ = ∑_j a_j e^{λ_j t}, \qquad
τ_\mathrm{int}
    = \frac{∫_0^∞ [\bar n(t) - \bar n_∞]\,\mathrm{d}t}{\bar n(0) - \bar n_∞}
    = -\frac{∑_j a_j/λ_j}{∑_j a_j},
```

where the sum runs over decaying modes and the coefficients ``a_j`` depend on
both the initial state and the phonon-number operator.
[`IntegratedTransient`](@ref) returns ``τ_\mathrm{int}``, an effective cooling
time for that preparation. This avoids choosing a fit window for an
exponential, but retains the dependence on initial occupation
([Kulosa et al. (2023)](https://doi.org/10.1088/1367-2630/acd13b),
Eqs. (10)–(14), Fig. 2(c) and Fig. C1). [`LiouvillianSpectrum`](@ref) returns
``-1/\mathrm{Re}\,λ`` for a selected mode with appreciable phonon-number
weight. The selection is independent of the preparation; the mode's
contribution to an actual transient is not. A weak tail can set the late-time
exponent while contributing little to the integral. Since the ``a_j`` can be
negative or complex, the integral is not generally an average of modal times
with positive weights.

#### Choosing a method

- Use [`AdiabaticElimination`](@ref) for parameter scans, then check promising
  points with a joint model. Compare both occupation and time; agreement at
  weak drive does not establish accuracy near a saturated optimum.
- Use [`LiouvillianSpectrum`](@ref) for late-time occupation relaxation.
  Check the selected mode as `num_fock`, `num_eigenpairs` and
  `weight_threshold` change. If the time keeps increasing with Fock size, it
  has not converged; evaluate cooling from the relevant preparation with
  [`IntegratedTransient`](@ref) or [`cooling_curve`](@ref).
- Use [`IntegratedTransient`](@ref) for an effective cooling time from a
  specified initial occupation, especially a hot state. Set `initial_nbar`
  to that preparation and allow [`fock_truncation`](@ref) to size the ladder.
  Use [`cooling_curve`](@ref) to verify the time to reach a target occupation
  number or the shape of the transient.

For either joint method, check convergence in Fock size, coherence bandwidth
and Lamb–Dicke order. A narrow `bandwidth` can greatly reduce memory and
factorisation cost; increase it until the observables stabilise. For a hot
state, use [`cooling_curve`](@ref) with `fock = true` to inspect the upper tail
during cooling. A first-order model can make a large ladder affordable, but
its accuracy must be checked over the occupied states.
When comparing the spectral and integral times, use the same model settings:
their defaults differ in both `lamb_dicke_order` and `bandwidth`.

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
