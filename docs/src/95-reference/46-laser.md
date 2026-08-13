# [Laser frequencies](@id reference-laser)

A laser frequency enters the package in one of three forms: a wavelength, an
absolute angular frequency, or a [`RelativeFrequency`](@ref) — an angular
offset from the zero-field interval of a named level pair.

The relative form exists for precision. Tabulated transition energies are
uncertain at the tens-of-MHz level, far coarser than the Zeeman and hyperfine
splittings the package computes exactly. Anchoring the laser to the line
centre lets that coarse tabulated value cancel out of every frequency
*difference* against components of the reference pair — as in the
near-resonant channels of [`light_shift`](@ref) — leaving only the stated
offset and exactly-known splittings. Only resolving the laser to an absolute
frequency via [`Levels.photon_energy`](@ref) reintroduces the tabulated value,
and with it its uncertainty.

The reference pair is a frequency marker, not a transition to be driven; for a
hyperfine species it names two ``F`` levels, whose Zeeman components are
degenerate at zero field.

```@autodocs
Modules = [Levels]
Pages = ["laser.jl"]
```
