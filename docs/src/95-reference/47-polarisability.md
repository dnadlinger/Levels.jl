# [Light shifts and polarisabilities](@id reference-polarisability)

The ac Stark (light) shift a probe laser induces on a transition — the offset
between the frequency at which resonant Rabi flopping is observed and the
unperturbed transition frequency — evaluated as a sum over the electric-dipole
channels of the two levels involved.

The atomic structure is fixed once per laser frequency by
[`LightShiftCoefficients`](@ref); evaluating [`light_shift`](@ref) for a
particular intensity and polarisation afterwards is cheap enough to sit inside
the inner loop of a laser-parameter fit. Which levels this is available for
depends on the [`LevelPolarisability`](@ref) data of the species; for
[`sr88`](@ref) it covers the S``_{1/2}`` ↔ D``_{5/2}`` clock transition.

Note that, unlike the electric-quadrupole Rabi frequency, an E1 light shift does
not depend on the beam direction — only on the polarisation. The two therefore
constrain different combinations of the beam parameters.

## Near-resonant quadrupole shifts

Those far-detuned dipole channels are not the whole story for a laser sitting on
a narrow quadrupole line. Driving one Zeeman component also couples both of its
states, off resonantly, to every other component sharing one of them, and the
detunings involved are mere Zeeman splittings. Passing the beam direction and
the magnetic flux density to [`light_shift`](@ref) as the `n` and `B` keyword
arguments adds this contribution; [`quadrupole_light_shift`](@ref) returns it on
its own, and — needing only the Einstein A coefficient rather than the
[`LevelPolarisability`](@ref) data — also covers the S``_{1/2}`` → D``_{3/2}``
line.

Being second order in the coupling, the shift scales as the intensity over the
field, so its size relative to the dipole shift is set by the field alone: for
⁸⁸Sr⁺ at 674 nm it is of the same order at 0.5 mT, but up to two orders of
magnitude larger at the 4.8 µT of an optical-clock field. What saves the clock
is a symmetry rather than a small size — with a linear polarisation the two
components of a ``±m`` Zeeman pair shift equally and oppositely, so the shift
cancels in the pair average the servo is steered to, and only the ellipticity
leaks through (`[Lindvall2025]`, Sec. III F 2).

The perturbative treatment requires every Rabi frequency to be small against the
Zeeman splittings, which is what breaks down as the field goes to zero.

## Hyperfine species

For hyperfine states, the electric-dipole channels are summed over the
hyperfine levels of each intermediate fine-structure level with the ``F``-basis
angular factors, and the intermediate splittings are resolved in the detunings
— which is what gives e.g. the ⁴³Ca⁺ S``_{1/2}`` ``F`` levels a small tensor
polarisability their ``J = 1/2`` parent lacks. The near-resonant quadrupole
shift becomes field-resolved: in the Breit–Rabi regime the spectator detunings
are exact eigen-energy differences rather than a common ``1/B`` scaling, so the
static field enters at [`LightShiftCoefficients`](@ref) construction (`B`
keyword). Note that the ``±m`` pair cancellation quoted above then survives
only approximately: spectators in *other* ``F`` levels sit at
hyperfine-interval detunings that are even under ``m → -m``, so their
contribution — often the dominant one — does not average out.

```@autodocs
Modules = [Levels]
Pages = ["polarisability.jl"]
```
