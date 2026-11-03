# Levels

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://dnadlinger.github.io/Levels.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://dnadlinger.github.io/Levels.jl/dev)
[![Test workflow status](https://github.com/dnadlinger/Levels.jl/actions/workflows/Test.yml/badge.svg?branch=master)](https://github.com/dnadlinger/Levels.jl/actions/workflows/Test.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/dnadlinger/Levels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dnadlinger/Levels.jl)
[![Lint workflow Status](https://github.com/dnadlinger/Levels.jl/actions/workflows/Lint.yml/badge.svg?branch=master)](https://github.com/dnadlinger/Levels.jl/actions/workflows/Lint.yml?query=branch%3Amaster)
[![Docs workflow Status](https://github.com/dnadlinger/Levels.jl/actions/workflows/Docs.yml/badge.svg?branch=master)](https://github.com/dnadlinger/Levels.jl/actions/workflows/Docs.yml?query=branch%3Amaster)

Levels.jl provides unitful atomic-structure and laser-dynamics calculations for
trapped-ion spectroscopy, with data for ⁴⁰Ca⁺, ⁴³Ca⁺, ⁸⁸Sr⁺, and ¹³⁸Ba⁺. It requires
Julia 1.12 or newer.

```julia
using Levels, Unitful
clock = StateSpec("S_1/2 F=4", 4) => StateSpec("D_5/2 F=4", 3)
B₀ = insensitive_field(ca43, clock, (0.45u"mT", 0.55u"mT"))  # 0.4958 mT
s = Eigenstructure(ca43, "S_1/2", "D_5/2"; B=B₀)
direction, ε = beam_vectors(π / 2, 0)  # k ⟂ B, ε ∥ B
Ω = rabi_frequency(s, clock, peak_intensity(10u"mW", 20u"µm"), ε, direction)
cyclic(abs(Ω))  # 197.3 kHz
```

This brief example finds the field at which a 729 nm component of ⁴³Ca⁺ is first-order
insensitive to field noise, and the Rabi frequency with which a Gaussian beam drives it
there (including the hyperfine mixing effects).

An `Eigenstructure` is the atom — a species, a basis of states and the static magnetic
field, solved once — on which every calculation is built. See the
[conventions](docs/src/05-conventions.md) for units, states, basis transformations and
manual operator construction.

## How to Cite

If you use Levels.jl in your work, please cite using the reference given in
[CITATION.cff](https://github.com/dnadlinger/Levels.jl/blob/master/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first take a look into our
[contributing guide directly on GitHub](docs/src/90-contributing.md) or the
[contributing page on the website](https://dnadlinger.github.io/Levels.jl/dev/90-contributing/)
