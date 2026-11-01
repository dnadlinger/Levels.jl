# [Bases and representations](@id bases)

A matrix's basis determines what its entries mean. This page describes the
index spaces and phase conventions used by Levels, and how to construct
operators in them.

## Three bases

| Basis | Labels | Where it appears |
| --- | --- | --- |
| Product (``m_J, m_I`` kets) | electronic and nuclear projections; storage order electronic slow, nuclear fast, each increasing | construction device only: `Levels.moment_operators` (the magnetic moment is a Kronecker sum there) and the tensor remainder of the polarisability |
| Coupled (``F, m_F`` kets; ``J, m_J`` for a fine-structure species) | zero-field quantum numbers; electron-first coupling ``⟨J m_J; I m_I ∣ F m_F⟩``, Condon–Shortley phases | `coupled_basis(s)`; every zero-field species-form result; the natural basis for an operator written down from quantum numbers |
| Field eigenbasis | adiabatic ``(F, m_F)`` labels at `field(s)` where the zero-field branches are resolved (``m_F`` exact) | `basis(s)`; operator factories and model inputs taking a structure |

For a fine-structure species the three coincide (``U = I``; ``W = I`` for a
canonical complete basis), so its angular amplitudes need no rotation in the
fixed-``J`` model; its energies and energy-dependent responses still depend on
the field. Partial or reordered bases still select or permute, even at zero.

## Two matrices

- ``U`` = `Levels.coupling_transform(species, level)`:
  ``U_{(m_J, m_I),(F, m_F)} = ⟨J m_J; I m_I | F m_F⟩``, real, nonzero only where
  ``m_J + m_I = m_F``, for one manifold. Product → coupled: ``X_c = U^† X_p U``.
- ``W`` = `eigenvectors(s)`: column ``k`` holds the coupled-basis components of
  the eigenstate labelled `basis(s)[k]`; rows are indexed by `coupled_basis(s)`.
  Coupled → eigen: ``X_e = W^† X_c W``, which is `to_eigenbasis(s, X_c)`.

Within one manifold, an operator in the product basis transforms as
``X_e = W^† U^† X_p U W``, with ``W`` restricted to that manifold's rows and
selected columns. For a structure spanning several manifolds, embed the
coupled operator into `coupled_basis(s)` before calling `to_eigenbasis`.

```@example bases
using Levels, LinearAlgebra, Unitful
s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=0.5u"mT")
W = eigenvectors(s)
size(W), maximum(abs, W' * W - I)
```

## The representation contract

> Operator factories and model constructors taking a structure use its
> eigenbasis, indexed by `basis(s)` (row and column ``i`` ↔ `basis(s)[i]`).
> Species-form operator factories use the coupled basis in the supplied
> `StateBasis` order. Conversion and diagnostic functions explicitly name both
> of their index spaces: `to_eigenbasis(s, X)` accepts `X` over
> `coupled_basis(s)`; `eigenvectors(s)` has coupled rows and selected
> eigenstate columns; `coupled_hamiltonian(s)` returns a coupled-basis
> operator.

Model constructors take their input matrices in that convention and perform
no automatic rotation. A plain matrix carries no basis metadata: equal sizes
cannot distinguish matrices from different bases, fields, orders or gauges.
Use one structure's factories throughout a calculation, or supply manual
operators as described [below](@ref manual-operators).

All Cartesian vectors use the same fixed laboratory axes, with the quantisation
axis along ẑ. `B` is signed: reversing it does not reverse the axes, relabel
``m`` or reinterpret a Jones vector's handedness. At `B = 0` the chosen axis
remains part of the representation convention. A perturbing field `b` is
expressed in the same frame; complex phasors go through `zeeman_drives`, with
``\vec B(t) = \mathrm{Re}[\vec b\, e^{-iΩt}]``.

## Conventions that fix signs

Three independent layers fix the sign of every hyperfine amplitude:
Condon–Shortley Clebsch–Gordan phases; the electron-first coupling order
(nuclear-first sources differ by ``(-1)^{I + J - F}`` per state; see
`Levels.hyperfine_reduction`); and the eigenvector sign. Under a whole-column
rephasing ``W → W D``, operators transform as ``X_e → D^† X_e D``, and state
vectors and density matrices transform with them: individual entries and their
relative phases change, consistently calculated observables and interference
predictions do not. Replacing individual eigenvector components by their
absolute values is not a valid phase repair.

Each eigenvector's component on its own zero-field label is made positive;
if that component is exactly zero, the largest-magnitude component is made
positive. For a canonical complete basis with resolved zero-field levels,
this gives ``W → I`` as ``B → 0``. A sign jump can occur when the chosen
component crosses zero, so signed amplitudes in a field sweep need not be
continuous. For example, the ⁴³Ca⁺ D``_{5/2}`` branch labelled
``(F, m_F) = (4, -1)`` has such a crossing near 1.29 mT.

Rabi-flopping measurements determine magnitudes, not coherent relative phases:
relative phases between several hand-entered connections must be supplied in
one consistent gauge.

## Zero field and degenerate branches

At `B = 0` the eigenstates are the coupled ``|F, m_F⟩`` states and the offsets
their Casimir shifts; for a complete canonical basis ``W = I``. Where the
zero-field ``F`` energies within each ``m_F`` block are distinct, this is also
the continuous limit of the field solutions:

```@example bases
deviation(b) = maximum(abs, eigenvectors(Eigenstructure(ca43, "S_1/2", "D_5/2"; B=b)) - I)
deviation(1e-3u"mT"), deviation(1e-6u"mT")
```

Degenerate zero-field ``F`` levels within a block (a species with vanishing
hyperfine constants but a nonzero nuclear moment, say) have no unique adiabatic
continuation: the finite-field eigenvectors are then fixed mixtures that do not
approach the identity, the labels follow a deterministic policy (ascending
``F`` among the degenerate levels, by energy order at the field), and the
structure makes no continuity claim for them. Diagonalising each exact-``m_F``
block separately keeps different ``m`` sectors from mixing at any degeneracy;
it does not remove degeneracy within a block.

## Partial and reordered bases

The eigenproblem is always solved on complete manifolds. `basis(s)` selects
which labelled eigenstates matrices address, and in which order: a matrix over
a partial basis holds the exact matrix elements between the retained
eigenstates, and omitting states truncates the *model* a consumer builds,
never the *solution*. ``W`` then has one orthonormal column per basis state
over the complete coupled rows, and `to_eigenbasis(s, X)` — which requires `X`
over `coupled_basis(s)`; an `n × n` matrix over a partial `basis(s)` is
refused — is the exact projection of a complete-manifold operator onto the
retained eigenstates. Scalar queries (`state_offset`, `state_energy`,
`zeeman_shift`, `transition_amplitude`, `rabi_frequency`, `intensity_for_rabi`,
`transition_frequency`, `zeeman_sensitivity`, `polarisability`) accept any
state of `coupled_basis(s)`, not only basis states; matrix queries are indexed
by `basis(s)` alone. `iscomplete(s, level)` tells whether a manifold is
present in full — `DrivenTransition` and `quasienergy_shift` require it, as
their within-manifold sums would silently change with a dropped state;
`lindblad_model` treats a dropped state as a modelling choice under its
`open_channels` policy.

Projection does not preserve products: for rectangular ``W``,
``(W^† X W)(W^† Y W) ≠ W^† X Y W`` in general, so squaring projected
angular-momentum operators or forming effective interactions through discarded
virtual states changes the physics. Exact projected matrix elements do not
imply exact truncated dynamics.

```@example bases
partial = Eigenstructure(ca43, [StateSpec("D_5/2 F=4", 3), StateSpec("S_1/2 F=4", 4)]; B=0.5u"mT")
size(eigenvectors(partial)), iscomplete(partial, "S_1/2"), state_offset(partial, StateSpec("S_1/2 F=3", 0))
```

## Checking a representation

The eigenvectors have orthonormal columns, ``W^† W = I``, and diagonalise the
centroid-relative Hamiltonian:

```math
W^† H_c W = \operatorname{diag}(\texttt{offsets(s)}).
```

For a whole fine-structure pair, a structure-form coupling equals the
transformed species-form coupling. The same identity holds for
`relative_couplings`, each `dipole_components` matrix, and
`zeeman_hamiltonian`. Coupling scales use the pair's centroid frequency and
Einstein A coefficient in both forms.

```@example bases
H = coupled_hamiltonian(s)
residual = maximum(abs, ustrip.(u"rad/µs", to_eigenbasis(s, H) - Diagonal(offsets(s))))
n, ε = beam_vectors(π / 2, π / 4)
C_e = coupling_matrix(s, "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n)
C_c = coupling_matrix(ca43, coupled_basis(s), "S_1/2" => "D_5/2", 1e4u"W/m^2", ε, n)
rotation = maximum(abs, ustrip.(u"rad/µs", C_e - to_eigenbasis(s, C_c)))
(; residual, rotation)
```

## Selecting ``F`` labels: rotate, then select

An ``F`` level named in a structure-form coupling factory selects the
eigenstates carrying that label **after** the whole fine-structure pair has
been rotated: ``P^e_{up} (W^† C W) P^e_{lo}``. Rotating a zero-field
``F``-restricted operator, ``W^† (P_{up} C P_{lo}) W``, is a different
operation — the eigenstates have components on other ``F`` levels — and the
two differ at the tens-of-percent level for the ⁴³Ca⁺ ``F = 4 → F = 4`` E2
block at 0.5 mT:

```@example bases
basis = Levels.basis(s)
r4s, r4d = staterange(basis, "S_1/2 F=4"), staterange(basis, "D_5/2 F=4")
selected = relative_couplings(s, "S_1/2 F=4" => "D_5/2 F=4", ε, n)[r4d, r4s]
whole = relative_couplings(s, "S_1/2" => "D_5/2", ε, n)[r4d, r4s]
rotated_restricted = to_eigenbasis(s, relative_couplings(ca43, coupled_basis(s), "S_1/2 F=4" => "D_5/2 F=4", ε, n))[r4d, r4s]
norm(selected - whole) / norm(whole), norm(rotated_restricted - whole) / norm(whole)
```

A `RelativeFrequency` pair is independently a frequency marker, not a
restriction on which components a beam drives.

## [Manual operators](@id manual-operators)

An operator written from zero-field quantum numbers uses `coupled_basis(s)`;
`to_eigenbasis` rotates and, for a partial basis, projects it:

```@example bases
cb = coupled_basis(s)
X = zeros(ComplexF64, length(cb), length(cb)) .* u"rad/µs"
X[stateindex(cb, "D_5/2 F=4", 1), stateindex(cb, "S_1/2 F=4", 0)] = angular(100u"kHz")
C = to_eigenbasis(s, X)
abs(C[stateindex(basis, "D_5/2 F=4", 1), stateindex(basis, "S_1/2 F=4", 0)])
```

A measured coupling between labelled field eigenstates is entered directly
over `Levels.basis(s)`. To retain the relative phases set by a beam geometry,
scale its generated matrix instead:
`rabi_normalised(relative_couplings(s, pair, ε, n), Levels.basis(s), probe, Ω0)`.

For a product-basis operator in one manifold, first form `U' * X_p * U` with
`U = Levels.coupling_transform(species(s), level)`. Its index space is
`StateBasis(species(s), level)`; embed that block into `coupled_basis(s)`,
then call `to_eigenbasis`.

A manual `HarmonicDrive(X, phase)` takes a Hermitian angular-frequency matrix
over the structure's eigenbasis and represents `X*cos(ω*t + phase)`.
`zeeman_drives(s, b)` builds these matrices from a laboratory field phasor.
`quasienergy_shift(s, state_or_pair, ω, drives)` requires these drives to be
block-diagonal in the fine-structure manifolds, and each probed manifold to be
complete and contiguous in the basis. Couplings between manifolds raise an error.
`ac_zeeman_shift(s, state_or_pair, ω, b)` constructs magnetic drives from `b`.
The keyword constructor of `DrivenTransition`, `BandLiouvillian(matrix, …)`
and `quasienergy_shift(H0, H1, …)` also support plain-matrix models; the caller
supplies a common basis and ordering for their inputs.

## Inspecting eigenstates and sweeping fields

`eigenvectors(s)[:, stateindex(Levels.basis(s), state)]` gives an eigenstate's
coupled-basis composition. Summing squared components within an ``F`` level
gives its population in that level:

```@example bases
v = eigenvectors(s)[:, stateindex(basis, "D_5/2 F=4", 3)]
[(level, sum(abs2, v[staterange(cb, level)])) for level in cb.levels if sum(abs2, v[staterange(cb, level)]) > 1e-6]
```

A `RelativeFrequency(s, probe, δ)` stores the laser frequency determined at
construction. Reuse it across structures for a fixed laser; rebuild it at
each field for a fixed detuning from the moving resonance.

## Physical dipole components

[`dipole_components`](@ref) returns the three E1 operator components in
coulomb metres. Contracting them with the beam geometry recovers
`rabi_frequency`; its docstring defines the reduced-element and interaction
sign conventions. Rate data fixes the magnitude of each reduced element,
but supplies no relative sign between different fine-structure pairs.

```@example bases
sp = Eigenstructure(ca43, "S_1/2", "P_1/2"; B=0.5u"mT")
D = dipole_components(sp, "S_1/2" => "P_1/2")
lo, up = StateSpec("S_1/2 F=4", 4), StateSpec("P_1/2 F=4", 3)
intensity = 1.0u"W/m^2"
E0 = sqrt(2 * intensity / (u"ε0" * u"c"))
εσ = [1.0, -im, 0.0] / sqrt(2)
g = dipole_geometry(εσ)
Ω_contracted = uconvert(u"rad/µs", E0 / u"ħ" * sum(g[q+2] * D[q+2][stateindex(Levels.basis(sp), up), stateindex(Levels.basis(sp), lo)] for q in -1:1))
Ω_contracted ≈ rabi_frequency(sp, lo => up, intensity, εσ, [0, 0, 1.0])
```

The components are the extension interface for custom second-order sums. The
S``_{1/2}`` polarisability at 674 nm from the explicit S → P channels plus the
lumped static remainder of the data, with the structure's energies in the
denominators:

```@example bases
s3 = Eigenstructure(sr88, "S_1/2", "P_1/2", "P_3/2"; B=0.5u"mT")
b3 = Levels.basis(s3)
state = StateSpec("S_1/2", 1//2)
λ = 674.025u"nm"
ħω = photon_energy(λ)
_, ε_lin = beam_vectors(0.4, 0.6)
w = Levels.polarisation_weights(ε_lin)
function custom_sum()
    α = level_polarisability(sr88, "S_1/2").static_scalar
    for upper in ("P_1/2", "P_3/2"), (k, t) in enumerate(b3)
        t.level == NoHyperfineNumberSpec(upper) || continue
        Dq = dipole_components(s3, "S_1/2" => upper)
        Δ = state_energy(s3, t) - state_energy(s3, state)
        for q in -1:1
            d = Dq[q+2][k, stateindex(b3, state)]
            α += uconvert(u"C*m^2/V", abs2(d) * (w[q+2] / (Δ - ħω) + w[2-q] / (Δ + ħω)))
        end
    end
    α
end
custom_sum() / polarisability(s3, state, λ, ε_lin)
```
