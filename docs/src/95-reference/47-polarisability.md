# [Light shifts and polarisabilities](@id reference-polarisability)

The light-shift model combines an electric-dipole (E1) **background** with
**near-resonant channels** resolved into individual Zeeman components. The
background uses explicit dipole strengths plus a static scalar/tensor
remainder; the channels use Einstein A coefficients and at-field couplings.
E1 shifts depend on polarisation. Electric-quadrupole (E2) channels also
depend on the beam direction.

| Quantity | API |
| --- | --- |
| E1 response of one state | [`polarisability`](@ref) |
| Shift of one state under a parked beam | [`light_shift`](@ref) |
| Displacement of a resonantly driven transition | [`driven_light_shift`](@ref) |
| Data cached for repeated shift evaluations | [`LightShiftCoefficients`](@ref) |
| Zero-field level decomposition | [`scalar_polarisability`](@ref), [`vector_polarisability`](@ref), [`tensor_polarisability`](@ref) |

## State polarisability

[`polarisability`](@ref) returns the E1 response of a state in SI units
(``C m^2/V``). With no polarisation argument, it returns three entries for
pure ``σ^-``, ``π`` and ``σ^+`` fields, indexed by ``q + 2``. A polarisation
argument contracts these with the normalised weights ``w_q = |ε_{-q}|^2``.
The ``π`` entry describes a field along the quantisation axis; the mean of
the two ``σ`` entries describes a linear field perpendicular to it.

A wavelength or angular frequency selects the dynamic response;
`0u"s^-1"` selects the static response. For example, the static differential
polarisability sets the electric Stark shift from micromotion in a Paul trap:

```@example polarisability
using Levels, Unitful

s = StateSpec("S_1/2 F=4", 4)
d = StateSpec("D_5/2 F=4", 3)
B0 = insensitive_field(ca43, s => d, (0.45u"mT", 0.55u"mT"))
Δα = polarisability(ca43, d, 0u"s^-1"; B=B0) -
     polarisability(ca43, s, 0u"s^-1"; B=B0)
Δα_parallel = Δα[2]
Δα_perpendicular = (Δα[1] + Δα[3]) / 2
θ = π / 4
Δα_tilted = cos(θ)^2 * Δα_parallel + sin(θ)^2 * Δα_perpendicular
(; Δα_parallel, Δα_perpendicular, Δα_tilted)
```

For a sinusoidal electric field with unit polarisation and peak amplitude
``E_0``, the cycle-averaged energy shift is ``-α E_0^2/4``. For a static field
it is ``-α E^2/2``. For a light beam of intensity ``I``, the angular-frequency
shift is ``-I α/(2c ε_0 ħ)``.

Without `B`, states are zero-field ``|J,m_J⟩`` or ``|F,m_F⟩`` states. With
`B`, hyperfine states denote the adiabatically labelled field eigenstates,
including ``F`` mixing. Fine-structure eigenvectors remain unchanged, but
their Zeeman energies enter the denominators. The scalar/vector/tensor level
accessors describe zero-field levels; their single-``F`` angular factors do
not apply to mixed field eigenstates.

## Parked and driven beams

For a parked beam, [`light_shift`](@ref) returns the angular-frequency shift
of one state. The shift of a transition in a Ramsey measurement is the upper
state's shift minus the lower state's shift.

A nearby laser is specified by [`RelativeFrequency`](@ref), an angular
offset from a named zero-field interval. The uncertain tabulated line centre
then cancels from the component detunings, leaving the offset and the
Zeeman/hyperfine splittings. See [Laser frequencies](@ref reference-laser).

```@example polarisability
laser = RelativeFrequency(s.level => d.level, 2π * 20u"MHz")
c = LightShiftCoefficients(ca43, [s, d], laser; B=B0)
intensity = 1u"W/m^2"
n, ε = beam_vectors(π / 2, π / 4, 0.3)
Δω_parked = light_shift(c, d, intensity, ε; n) -
            light_shift(c, s, intensity, ε; n)
```

For a resonantly driven component, [`driven_light_shift`](@ref) returns the
displacement of its resonance. It evaluates spectator channels with the laser
tuned to the probed component and excludes that component as the coherent
drive. A positive result raises the resonance frequency.

```@example polarisability
Δω_driven = driven_light_shift(c, s => d, intensity, ε; n)
```

Both APIs accept `parts = :total | :background | :resonant` (default
`:total`). The resonant part needs `B` and the Einstein A coefficient, but
no polarisability data. E2 channels are built for pairs within the basis;
E1 channels require that pair to be the `RelativeFrequency` reference.

Coefficients cache the atomic structure for repeated intensity, polarisation
and detuning evaluations. A parked evaluation's `δ` keyword changes only the
near-resonant detunings; the smooth background stays at the construction
frequency. Coefficients used in driven mode must have a construction
frequency within relative tolerance ``10^{-3}`` of the probed line.

## Background and channel split

For an E1 reference pair, the near-resonant term moves from the background
into the channel sum. The counter-rotating term stays in the background,
so the split counts every term once at the construction frequency. The
reference channel must be explicit in the [`LevelPolarisability`](@ref)
data; a lumped remainder cannot be split.

Bare frequencies within ~100 GHz of an explicit E1 line are refused: precise
component detunings require a reference that cancels the uncertain line
centre. The state [`polarisability`](@ref) accessor has the same guard and
rejects `RelativeFrequency`; the level scalar/vector/tensor accessors remain
unguarded curves.

For hyperfine states, explicit E1 channels resolve the known hyperfine
splittings of both levels. This can produce a tensor response even for an
``F`` level of a ``J = 1/2`` parent. With `B`, amplitudes rotate into both
field eigenbases before squaring, and denominators use their eigen-energies.
The static remainder is an electronic operator evaluated in the source
eigenstate, retaining coherences between its ``F`` components.

An intermediate level without hyperfine constants is approximated as
degenerate at its centroid, with no Zeeman shift. Field-resolved source
levels and near-resonant partners require known constants.

## Validity and symmetry

The background requires polarisability data: S``_{1/2}`` and D``_{5/2}`` for
[`sr88`](@ref), and also D``_{3/2}`` for [`ca40`](@ref) and [`ca43`](@ref).
Evaluating a background for a level without data raises an error; its
near-resonant channels can still be used with `parts = :resonant`.
The static remainder approximates unresolved E1 channels and is suitable
only far from their resonances. E2 contributions outside the selected
near-resonant pair and counter-rotating E2 terms are omitted.

The channel sum uses second-order perturbation theory with unregularised
``1/detuning`` poles. Evaluation warns if a channel's Rabi frequency exceeds
one tenth of its detuning or the laser lies within five upper-level
linewidths of a pole. Linewidths enter these warnings, not the denominators;
the model provides no scattering rate.

For fine-structure E2 transitions, spectator detunings scale with `B`, so
the resonant shift scales as ``I/B``. With linear polarisation, reversing
both magnetic projections changes this part's sign; it cancels in a
``±m`` pair average (`[Lindvall2025]`, Sec. III F 2). Ellipticity breaks that
cancellation. The E1 background generally remains in the pair average.

Hyperfine spectators in other ``F`` levels generally prevent the resonant
pair cancellation and the ``I/B`` scaling. For linear polarisation, the
field-reversal identity is ``shift(-m,B) = shift(+m,-B)``, with both
projections reversed; for elliptical polarisation the field reversal also
requires conjugating `ε`.

```@autodocs
Modules = [Levels]
Pages = ["polarisability.jl"]
```
