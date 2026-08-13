# [Light shifts and polarisabilities](@id reference-polarisability)

The ac Stark (light) shift of the Zeeman states of a laser-illuminated ion,
composed of two parts: the far-detuned electric-dipole **background** — the
per-state polarisability at the laser frequency, summed over the explicit
dipole channels of the level — and the **near-resonant channels** of one level
pair the laser sits close to, resolved down to the individual Zeeman
components at the static field.

The atomic structure is fixed once by [`LightShiftCoefficients`](@ref);
evaluating a shift for a particular intensity and polarisation afterwards is
cheap enough to sit inside the inner loop of a laser-parameter fit. Which
levels the background is available for depends on the
[`LevelPolarisability`](@ref) data of the species (for [`sr88`](@ref), the
S``_{1/2}`` ↔ D``_{5/2}`` clock levels); the channels need only the Einstein A
coefficient.

Note that, unlike the electric-quadrupole coupling, an E1 light shift does not
depend on the beam direction — only on the polarisation. The two therefore
constrain different combinations of the beam parameters.

## Two observables

Driving one Zeeman component of a transition also couples the two states
involved, off resonantly, to every other component sharing one of them, and
the detunings involved are mere Zeeman splittings.
[`driven_light_shift`](@ref) returns the displacement of the *driven*
resonance: the laser is pinned to the probed component — so the laser
frequency drops out entirely — and that component itself is excluded as the
coherent drive. This is the shift ⁸⁸Sr⁺ clock evaluations quote as the "E2 ac
Stark shift", and — needing no polarisability data — it covers the 687-nm
S``_{1/2}`` → D``_{3/2}`` line just as well.

A beam *parked* at a stated frequency, resonantly driving nothing, instead
shifts every state it illuminates: [`light_shift`](@ref) with the laser given
as a [`RelativeFrequency`](@ref) — an offset from the zero-field interval of a
named level pair — adds the near-resonant channel sum at the exact at-field
detunings to the background. The parametrisation is what preserves precision:
the tabulated line centre cancels between the laser and the components,
leaving only the stated offset and exactly-known Zeeman and hyperfine
splittings (see [Laser frequencies](@ref reference-laser)). This is the
regime of a laser placed *within* the
manifold of Zeeman components, poles and all; second-order perturbation
theory holds while every channel Rabi frequency stays small against its
detuning, and the evaluation warns otherwise. The Ramsey-type shift of a pair
from a parked beam is simply the difference of the two states' shifts.

Being second order in the coupling, the driven shift scales as the intensity
over the field, so its size relative to the background is set by the field
alone: for ⁸⁸Sr⁺ at 674 nm it is of the same order at 0.5 mT, but up to two
orders of magnitude larger at the 4.8 µT of an optical-clock field. What saves
the clock is a symmetry rather than a small size — with a linear polarisation
the two components of a ``±m`` Zeeman pair shift equally and oppositely, so
the shift cancels in the pair average the servo is steered to, and only the
ellipticity leaks through (`[Lindvall2025]`, Sec. III F 2).

## Electric-dipole lines

The channel machinery is rank-agnostic: naming an electric-dipole pair in the
`RelativeFrequency` (e.g. 422 nm S``_{1/2}`` → P``_{1/2}``) moves the
co-rotating part of that one polarisability channel out of the background and
into the Zeeman-resolved sum — an exact split, with the smooth
counter-rotating term staying behind. This is the regime of a beam parked
some tens of MHz to a few GHz from a broad dipole line, where tens-of-MHz
Zeeman splittings modify the result at a level a Zeeman-blind polarisability
cannot see. A *bare* wavelength that close (~100 GHz) to an explicit dipole
channel is refused outright for the same reason.

## Hyperfine species

For hyperfine states, the electric-dipole background channels are summed over
the hyperfine levels of each intermediate fine-structure level with the
``F``-basis angular factors, and the intermediate splittings are resolved in
the detunings — which is what gives e.g. the ⁴³Ca⁺ S``_{1/2}`` ``F`` levels a
small tensor polarisability their ``J = 1/2`` parent lacks. The near-resonant
channels are built from the eigenbasis-rotated amplitudes and exact
intermediate-field eigen-energy differences of the manifolds at the
construction field, so the ``F`` mixing is accounted for exactly; the
reference pair of a `RelativeFrequency` names two ``F`` levels, whose
zero-field interval — where the Zeeman components are degenerate — anchors
the offset. Note that the ``±m`` pair cancellation quoted above then survives
only approximately: spectators in *other* ``F`` levels sit at
hyperfine-interval detunings that are even under ``m → -m``, so their
contribution — often the dominant one — does not average out.

```@autodocs
Modules = [Levels]
Pages = ["polarisability.jl"]
```
