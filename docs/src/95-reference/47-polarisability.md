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

```@autodocs
Modules = [Levels]
Pages = ["polarisability.jl"]
```
