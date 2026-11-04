# [Conventions](@id conventions)

## [Frequencies, rates and energies](@id units)

All frequencies, including Rabi frequencies, detunings and linewidths, are angular
frequencies, returned in `u"rad/µs"`; Hamiltonian matrices represent ``H/ħ`` in the same
units. For an excited state of lifetime ``τ`` decaying to a stable lower state, the
natural Lorentzian linewidth (FWHM) in angular units is numerically equal to the decay
rate ``τ^{-1}``.

Angular frequencies carry Unitful's `rad` as a marker. The marker is dimensionless and
therefore soft: Unitful identifies Hz, s⁻¹ and rad s⁻¹ with each other and cannot detect
a missing factor of 2π, and a sum of a `rad µs⁻¹` and a plain `µs⁻¹` quantity silently
drops the marker. Use `angular(50u"MHz")` for a stated cyclic frequency, and `cyclic(ω)`
to display an angular result in Hz. These conversions change the physical convention;
`uconvert` only changes the units.

`level_energy(species, level)` and `state_energy(s, state)` return energies in joules
relative to the species' energy reference. `state_offset(s, state)` and `offsets(s)`
return angular-frequency offsets from each fine-structure centroid. Use
`transition_frequency(s, lower => upper)` for an interval: it keeps the centroid and
offset terms separate, preserving small splittings within a manifold.

`saturation_intensity(species, lower => upper)` returns ``I_0 = ħω³/(6πc²τ)``, with
``τ`` the total upper-level lifetime. At this intensity a closed stretched E1 transition
has ``|Ω| = Γ = 1/τ``. This is twice the two-level ``I_{\mathrm{sat}}`` defined by
``s = I/I_{\mathrm{sat}} = 2|Ω|²/Γ²``.

## [States, fields and the structure](@id bases)

Levels accept strings such as `"D_5/2"`; hyperfine levels include an ``F`` suffix, for
example `"D_5/2 F=4"`. `StateSpec(level, m)` uses ``m_J`` for a fine-structure level and
``m_F`` for a hyperfine level. Invalid magnetic projections are rejected.
`StateBasis(species, levels...)` expands each level into increasing-``m`` substates; for
a hyperfine species a fine-structure shorthand expands all its ``F`` levels. Explicit
state lists fix any desired ordering.

Transitions use `lower => upper` pairs. The static magnetic field is specified in
exactly one place: the constructor of an [`Eigenstructure`](@ref),
`Eigenstructure(species, basis; B)`, with `B` a signed scalar flux density along the
quantisation axis ẑ (zero allowed, no default). The structure solves every
fine-structure manifold the basis touches:

```@example conventions
using Levels, LinearAlgebra, Unitful
B = 0.5u"mT"
s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B)
```

All Cartesian vectors use fixed laboratory axes with the quantisation axis along ẑ.
Reversing `B` leaves these axes, the ``m`` labels and polarisation handedness unchanged;
the same axes apply at zero field. Oscillating magnetic fields use phasors in this
frame, with ``\vec B_{\mathrm{ac}}(t) = \mathrm{Re}[\vec b\, e^{-iΩt}]``.

**Species forms**, such as `transition_amplitude(ca43, pair)`, use zero-field states.
**Structure forms**, such as `transition_amplitude(s, pair)`, evaluate the labelled
eigenstates at `field(s)`. Operator factories and model constructors taking `s` use its
eigenbasis, indexed by `Levels.basis(s)`; species-form operators use the supplied
coupled basis.

Vectors and operators are plain Julia arrays without basis metadata. Users must
therefore keep their basis, field, ordering and eigenvector phases consistent: equal
matrix sizes cannot detect a mismatch. Using one `Eigenstructure` instance throughout a
calculation keeps the built-in factories and model inputs consistent. The accessors and
transformations below support manually constructed operators. Treat species, bases,
structures and models, including their arrays, as read-only; accessors return stored
data directly.

### Product, coupled and field eigenbases

Hyperfine coupling and the static field relate three bases, each suited to a different
part of the calculation:

- The **product basis** ``|m_J, m_I⟩`` is used internally to construct magnetic moments
  and polarisability operators. Electronic projection varies slowest, nuclear projection
  fastest, both increasing.
- The **coupled basis** ``|F, m_F⟩`` uses zero-field quantum numbers; `coupled_basis(s)`
  contains every state of each manifold, with increasing ``F`` and then increasing
  ``m_F``. For a species without hyperfine structure, these are the ``|J, m_J⟩`` states.
- The **field eigenbasis**, indexed by `Levels.basis(s)`, contains the selected
  eigenstates at `field(s)`. For a hyperfine species, ``m_F`` remains exact; ``F``
  labels the branch connected to that zero-field level, rather than a conserved angular
  momentum.

### Transformations between bases

Two matrices relate these bases and determine how manually constructed operators enter a
structure:

- `U = Levels.coupling_transform(species, level)` contains the Clebsch–Gordan
  coefficients ``⟨J m_J; I m_I | F m_F⟩`` for one manifold. A product-basis operator
  becomes ``X_c = U^† X_p U`` in the coupled basis.
- `W = eigenvectors(s)` has rows indexed by `coupled_basis(s)`; column ``k`` contains
  the components of the eigenstate `Levels.basis(s)[k]`. Thus ``v_c = W v_e`` for state
  vectors and ``X_e = W^† X_c W`` for operators. The latter transformation is
  [`to_eigenbasis(s, X_c)`](@ref to_eigenbasis).

For a product-basis operator, embed ``U^† X_p U`` into the corresponding manifold block
of `coupled_basis(s)` before applying `to_eigenbasis`. Without hyperfine structure, the
field does not mix states: ``W`` only selects and orders the ``|J, m_J⟩`` states. It is
the identity matrix when `Levels.basis(s)` contains every state in exactly the order of
`coupled_basis(s)`.

### Partial and reordered bases

The eigenproblem is always solved on complete manifolds; `Levels.basis(s)` selects and
orders the eigenstates that operator matrices address. For a partial basis, ``W`` is
rectangular and `to_eigenbasis` projects a complete coupled-basis operator onto the
retained eigenstates. Scalar queries can still address any state of `coupled_basis(s)`.

These projected matrix elements are exact, but truncated dynamics need not be:
projection generally does not preserve products, ``(W^† X W)(W^† Y W) ≠ W^† X Y W``. In
particular, discarded states no longer contribute as intermediate states.
[`iscomplete`](@ref) identifies complete manifolds and documents which models require
them.

### Checking representations

The columns of ``W`` are orthonormal, ``W^† W = I``, and diagonalise the
centroid-relative Hamiltonian returned by `coupled_hamiltonian(s)`:

```math
W^† H_c W = \operatorname{diag}(\texttt{offsets(s)}).
```

This provides a check on the transformation and units:

```@example conventions
H_c = coupled_hamiltonian(s)
maximum(abs, ustrip.(u"rad/µs", to_eigenbasis(s, H_c) - Diagonal(offsets(s))))
```

For a whole fine-structure pair, transforming a species-form coupling gives the
structure-form coupling. Both use the same centroid frequency and Einstein A coefficient
to set the physical scale.

### Selecting ``F`` labels

An ``F`` restriction in a structure-form coupling selects the labelled eigenstates
**after** rotating the whole fine-structure pair. Restricting the zero-field operator
first discards components on other ``F`` levels that contribute through field-induced
mixing. These operations differ by tens of percent for the ⁴³Ca⁺ ``F = 4 → F = 4`` E2
block at 0.5 mT.

### [Manual operators](@id manual-operators)

An operator defined through zero-field quantum numbers must be constructed over
`coupled_basis(s)`, then transformed with `to_eigenbasis`. For example:

```@example conventions
cb = coupled_basis(s)
X = zeros(ComplexF64, length(cb), length(cb)) .* u"rad/µs"
X[stateindex(cb, "D_5/2 F=4", 1), stateindex(cb, "S_1/2 F=4", 0)] = angular(100u"kHz")
C = to_eigenbasis(s, X)
abs(C[stateindex(Levels.basis(s), "D_5/2 F=4", 1), stateindex(Levels.basis(s), "S_1/2 F=4", 0)])
```

A measured coupling between labelled field eigenstates is instead entered directly over
`Levels.basis(s)`. To preserve the relative phases from a beam geometry while setting
one measured Rabi frequency, scale the generated matrix with [`rabi_normalised`](@ref).

### Eigenvector phases

Hyperfine amplitudes use Condon–Shortley phases and electron-first coupling;
nuclear-first states differ by ``(-1)^{I+J-F}`` per state. Each field eigenvector is
chosen with a positive component on its own zero-field label; if that component
vanishes, its largest-magnitude component is made positive. This local choice allows
sign jumps in a field sweep, for example near 1.29 mT for the ⁴³Ca⁺ ``D_{5/2}`` branch
labelled ``(F, m_F) = (4, -1)``.

For a diagonal phase matrix ``D``, rephasing eigenvectors as ``W → W D`` transforms
operators as ``X_e → D^† X_e D``. Observables are unchanged when states and operators
use the same phases. Rabi-flopping measurements determine magnitudes; coherent relative
phases between manually entered couplings require a consistent choice of eigenvector
phases.

## Beam geometry and detunings

`beam_vectors` returns `(direction, ε)`: a propagation direction and complex
polarisation. The electric field is ``\vec E(t) = \mathrm{Re}[\vec ε E_0 e^{-iωt}]``; a
right-circular beam along ẑ drives ``Δm = +1``. Rabi frequencies retain their complex
phases. The coupling Hamiltonian is `(C + C')/2`, with `C` containing angular Rabi
frequencies in the **upper-row, lower-column** block. [`dipole_components`](@ref)
documents the dipole and interaction sign conventions.

Detuning is laser minus atomic frequency. `RelativeFrequency(level_pair, offset)` uses
the zero-field level interval as its marker;
`RelativeFrequency(s, state_pair, detuning)` sets the frequency relative to that
component at `field(s)`. Reuse the latter across structures for a fixed laser frequency;
rebuild it at each field for fixed detuning from the moving resonance. The marker pair
does not restrict which components a beam drives.

Spherical channels use index `q + rank + 1`, hence `q + 2` for E1 and M1 and `q + 3`
for E2. The multipole of a transition is that of its Einstein A entry
([`Levels.multipole`](@ref)), never inferred from parity alone; an M1 transition couples
to the beam's magnetic polarisation ``\hat n × \vec ε``.
Fock occupation ``n`` starts at zero; `num_fock` counts retained states, so the maximum
occupation is `num_fock - 1`.

## Zero field and degenerate branches

At `B = 0`, the eigenstates are the coupled states. If `Levels.basis(s)` contains every
state in the order of `coupled_basis(s)`, ``W = I``: no rotation or reordering is
required. For distinct zero-field ``F`` energies within each ``m_F`` block, this is also
the limit as ``B → 0``.

Degenerate ``F`` levels within a block have no unique adiabatic continuation. Their
finite-field eigenvectors need not approach the coupled states; labels follow a
deterministic energy-ordering policy, with zero-field ties ordered by increasing ``F``.
Different ``m_F`` sectors remain unmixed because they are diagonalised separately.

## [Fine-structure mixing](@id fine-structure-mixing)

Levels.jl currently neglects fine-structure mixing; calculations are done manifold by
manifold, each spanning a single fine-structure level ``J``. A given static magnetic
field mixes the ``F`` levels within each manifold, but the magnetic dipole coupling
between fine-structure partners such as D₃/₂ and D₅/₂ is not considered, even though we
do model M1 decays (see [`MagneticDipole`](@ref)). The quadratic Zeeman shift that this
mixing causes is therefore absent. For the doublet of orbital angular momentum ``l``, it
is in LS coupling

```math
Δf = \frac{(g_s - 1)^2 (μ_B B/h)^2}{Δf_\mathrm{FS}}
    \, \frac{(l + 1/2)^2 - m^2}{(2l + 1)^2},
```

raising the upper level and lowering the lower one by the same amount, with stretched
states unshifted. For D₅/₂ and the ``m = ±1/2`` states, where it is largest:

| ion | ``Δf_\mathrm{FS}`` / THz | ``(Δf / B^2)`` / Hz mT⁻² |
|:----|----:|----:|
| ⁴⁰Ca⁺ | 1.82 | +26 |
| ⁸⁸Sr⁺ | 8.40 | +5.6 |
| ¹³⁸Ba⁺ | 24.0 | +2.0 |

The same terms are missing from [`zeeman_sensitivity`](@ref) and
[`insensitive_field`](@ref), and from the ac Zeeman shifts that rf magnetic fields cause
in `Levels.PeriodicDynamics`. For hyperfine isotopes, the curvature from ``F`` mixing is
larger by orders of magnitude, so field-insensitive points barely move, but transition
frequencies still acquire these Hz-level offsets. At rf-trap fields the mixing is not
usually of significance, of relative order ``μ_B B / hΔν_\mathrm{FS} ≲ 10^{-5}``; but in
Penning traps it often is (0.035 for the Ca⁺ 3D doublet at 4.5 T).
