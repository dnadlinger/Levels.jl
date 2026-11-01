# [Laser frequencies](@id reference-laser)

A laser frequency may be a wavelength, an absolute angular frequency, or a
[`RelativeFrequency`](@ref). The relative form anchors the laser to a named
zero-field level interval. When forming detunings from components of that
pair, the tabulated centroid interval cancels, leaving the stated offset and
hyperfine/Zeeman splittings. Their accuracy is set by the species data and
the level-structure model.

The reference pair is a frequency marker. Beam couplings are evaluated over
the whole fine-structure pair, including all its ``F`` levels for a hyperfine
species.

```@autodocs
Modules = [Levels]
Pages = ["laser.jl"]
```
