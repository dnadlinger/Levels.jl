# [Conventions and examples](@id conventions)

## [Frequencies, rates and energies](@id units)

| Quantity | Convention |
| --- | --- |
| Optical and trap frequencies, Rabi frequencies and detunings | Angular frequency, returned in `rad µs⁻¹` (`u"rad/µs"`) |
| Einstein A, population decay, cooling and heating rates | Events or quanta per time; no factor of 2π |
| Laser linewidth | Lorentzian FWHM in angular-frequency units |
| Hyperfine constants `a`, `b` | Energies: multiply a quoted cyclic frequency by `h` |
| Hamiltonian matrices | `H/ħ`, in angular-frequency units |
| Jump operators | Square root of a rate |
| Electric-dipole components | Coulomb metres (`dipole_components`) |

Angular frequencies carry Unitful's `rad` as a marker. The marker is
dimensionless and therefore soft: Unitful identifies Hz, s⁻¹ and rad s⁻¹ with
each other and cannot detect a missing factor of 2π, and a sum of a `rad µs⁻¹`
and a plain `µs⁻¹` quantity silently drops the marker. Use `angular(50u"MHz")`
for a stated cyclic frequency, and `cyclic(ω)` to display an angular result in
Hz. Do not apply either conversion to an Einstein A or a heating rate, which
stay plain inverse times. These conversions change the physical convention;
`uconvert` only changes the units.

`level_energy(species, level)` and `state_energy(s, state)` return energies in
joules relative to the species' energy reference. `state_offset(s, state)` and
`offsets(s)` return angular-frequency offsets from each fine-structure
centroid. Use `transition_frequency(s, lower => upper)` for an interval: it
keeps the centroid and offset terms separate, preserving small splittings
within a manifold.

`saturation_intensity(species, lower => upper)` returns `I_0 = ħω³/(6πc²τ)`,
with `τ` the total lifetime of the upper level. For a fine-structure E1 pair,
a component driven at unit geometric amplitude then has `|Ω|² = C² A/τ`, with
`C` its Clebsch–Gordan coefficient and `A` the pair's Einstein coefficient. On a closed
stretched transition, `|Ω| = Γ = 1/τ`. This is twice the two-level `I_sat` of the
`s = I/I_sat = 2|Ω|²/Γ²` saturation-parameter convention. The saturation
parameter of an individual component also depends on its angular and geometric
strength.

## States, fields and the structure

Levels accept strings such as `"D_5/2"`; hyperfine levels include an F suffix,
for example `"D_5/2 F=4"`. `StateSpec(level, m)` uses m_J for a fine-structure
level and m_F for a hyperfine level. Invalid magnetic projections are rejected.
`StateBasis(species, levels...)` expands each level into increasing-m substates;
for a hyperfine species a fine-structure shorthand expands all its F levels.
Explicit state lists fix any desired ordering.

Transitions use `lower => upper` pairs. The static magnetic field is specified
in exactly one place: the constructor of an [`Eigenstructure`](@ref),
`Eigenstructure(species, basis; B)`, with `B` a signed scalar flux density
along the quantisation axis ẑ (zero allowed, no default). The structure solves
every fine-structure manifold the basis touches and is the atom all further
calculations take:

```@example conventions
using Levels, Unitful
B = 0.5u"mT"
s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B)
```

**Species forms** use zero-field states, for example
`transition_amplitude(ca43, pair)` or `polarisability(ca43, state, ω)`.
**Structure forms** evaluate the labelled eigenstates at `field(s)`, for
example `transition_amplitude(s, pair)` or `polarisability(s, state, ω)`.
Model constructors take a structure: `LightShiftCoefficients(s, laser)`,
`DrivenTransition(s, probe; …)` and `lindblad_model(s, scheme)`.

Operators from a species form use the supplied coupled basis; operators from a
structure form use its field eigenbasis, indexed by `Levels.basis(s)`.
[Bases and representations](@ref bases) explains the eigenstate labels,
partial bases, phase conventions and construction of manual operators.

## Beam geometry and signs

`beam_vectors` returns `(direction, ε)`: a propagation direction and complex
polarisation. The field convention is `E(t) = Re[ε E₀ exp(-im*ω*t)]`; a
right-circular beam along ẑ drives Δm = +1. The complex phase of the Rabi
frequency is retained. The coupling Hamiltonian is `(C + C')/2`, with C
containing angular Rabi frequencies in the **upper-row, lower-column** block.
The physical dipole operator behind the Rabi frequency, and the sign convention
of `−d·E`, are documented at [`dipole_components`](@ref).

Detuning is laser minus atomic frequency. `RelativeFrequency(level_pair, offset)`
uses the zero-field level interval as its marker. The component form
`RelativeFrequency(s, state_pair, detuning)` folds the exact component shift
at the structure's field into that marker; it stores the frequency determined
then and does not follow the component into a structure at another field. A
light-shift sweep's `δ` overrides the marker offset; it is not necessarily
detuning from a particular Zeeman component.

Spherical channels use ordinary Julia indexing: index `q + rank + 1`, hence
`q + 2` for E1 and `q + 3` for E2. Fock occupation n starts at zero;
`num_fock` counts retained states, so the maximum occupation is `num_fock - 1`.

## Periodically driven transitions

The structure supplies the static energies; beam and field geometry supply
the laser coupling and ac drive matrices:

```@example conventions
using Levels.PeriodicDynamics
probe = StateSpec("S_1/2 F=4", 0) => StateSpec("D_5/2 F=4", 1)
direction, ε = beam_vectors(π/2, π/4)

# Structure-form matrices are in the field eigenbasis of s already.
C = coupling_matrix(s, "S_1/2" => "D_5/2", 1e6u"W/m^2", ε, direction)
drives = zeeman_drives(s, [10.0, 0.0, 0.0]u"µT")
dt = DrivenTransition(s, probe; drive_frequency=angular(50u"MHz"), coupling=C, drives)
nothing # hide
```

For a measured carrier Rabi frequency, use
[`rabi_normalised`](@ref) to scale the structure's relative coupling matrix.
For an operator written in the coupled basis, use [`to_eigenbasis`](@ref)
before passing it to a model; see [Manual operators](@ref manual-operators).

Treat species, bases, structures, schemes and models as read-only after
construction, including their arrays: accessors return the stored data directly.

## Numerical methods and physical results

```@example conventions
perturbative = sideband_rabi(dt, DressedStates(
    harmonic_cutoff=8, num_time_samples=256); sideband=1)
```

```@example conventions
reference = sideband_rabi(dt, Monodromy(reltol=1e-10); sideband=1)
```

Increase the harmonic cutoff and time sampling to check convergence. A cutoff K
retains `-K:K`, or 2K+1 harmonic blocks. `Monodromy` includes the laser
nonperturbatively within the model; its numerical accuracy is controlled by the
ODE and avoided-crossing search. A complex perturbative coupling amplitude is
available only from the dressed-state calculation.

Quasienergy shifts under harmonic driving use the same explicit method choice,
on the same structure. Ac Zeeman shifts are one application:

```@example conventions
drives = zeeman_drives(s, [10.0, 0.0, 0.0]u"µT")
shift = quasienergy_shift(s, probe, angular(50u"MHz"),
                         drives, PerturbationTheory())
```

```@example conventions
check = quasienergy_shift(s, probe, angular(50u"MHz"), drives, Floquet())
```

The pair form returns `(; shift, lower, upper)`; the state form returns a scalar
shift. [`quasienergy_shift_contributions`](@ref Levels.PeriodicDynamics.quasienergy_shift_contributions)
breaks down the perturbative sum by intermediate state.

Custom `HarmonicDrive`s can represent other interactions. These structure forms
require drives at one frequency, block-diagonal in the fine-structure manifolds,
with each probed manifold complete and contiguous in the basis. Couplings
between manifolds raise an error; a general model may instead use the
`quasienergy_shift(H0, H1, ω, i)` matrix form with full static energies.

For magnetic fields, [`ac_zeeman_shift`](@ref Levels.PeriodicDynamics.ac_zeeman_shift)
accepts the field phasor directly and constructs the drives:

```@example conventions
ac_zeeman_shift(s, probe, angular(50u"MHz"), [10.0, 0.0, 0.0]u"µT") == shift
```

For sideband `+1`, the dressed-state resonance minus the drive frequency is
the difference of the two dressed levels' ac shifts. `Floquet()` evaluates
those same shifts; `PerturbationTheory()` approximates them to second order.
The `Monodromy` resonance also includes probe-induced light shifts:

```@example conventions
perturbative.δ_res - dt.drive_frequency
```

## Cooling calculations

[`cooling_metrics`](@ref Levels.OpticalBloch.cooling_metrics) offers several
relaxation observables. Adiabatic elimination gives a rate-equation exponential
time constant; `SaturatedRates` and `IntegratedTransient` give an integral
relaxation time for a specified initial state, from a rate equation and a joint
model respectively; `LiouvillianSpectrum` gives an asymptotic relaxation
time; `SemiclassicalRates` gives the local rate of classical motion at a given
energy. They coincide for a single exponential and can differ for a hot ion or a
transient with several timescales. The [Optical Bloch reference](@ref
reference-optical-bloch) describes their assumptions and numerical controls.

Solver wrappers return unitful rates and times. Plain-matrix
`BandLiouvillian` constructors use a declared `time_unit` (µs by default);
model-based constructors inherit the model's scale. Numerical engines use
Float64/ComplexF64.
