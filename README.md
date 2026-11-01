# Levels

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://dnadlinger.github.io/Levels.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://dnadlinger.github.io/Levels.jl/dev)
[![Test workflow status](https://github.com/dnadlinger/Levels.jl/actions/workflows/Test.yml/badge.svg?branch=master)](https://github.com/dnadlinger/Levels.jl/actions/workflows/Test.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/dnadlinger/Levels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dnadlinger/Levels.jl)
[![Lint workflow Status](https://github.com/dnadlinger/Levels.jl/actions/workflows/Lint.yml/badge.svg?branch=master)](https://github.com/dnadlinger/Levels.jl/actions/workflows/Lint.yml?query=branch%3Amaster)
[![Docs workflow Status](https://github.com/dnadlinger/Levels.jl/actions/workflows/Docs.yml/badge.svg?branch=master)](https://github.com/dnadlinger/Levels.jl/actions/workflows/Docs.yml?query=branch%3Amaster)

Levels.jl provides unitful atomic-structure and laser-dynamics calculations for
trapped-ion spectroscopy, with data for ⁸⁸Sr⁺, ⁴⁰Ca⁺ and ⁴³Ca⁺. It requires Julia
1.12 or newer.

```julia
using Levels, Unitful
s = Eigenstructure(ca40, "S_1/2", "D_5/2"; B=0.5u"mT")
clock = StateSpec("S_1/2", -1//2) => StateSpec("D_5/2", -1//2)
cyclic(transition_frequency(s, clock))
```

An `Eigenstructure` is the atom — a species, a basis of states and the static
magnetic field, solved once — on which every calculation is built. See the
[conventions and examples](docs/src/05-conventions.md) for angular
frequencies, decay rates and transition pairs, and
[bases and representations](docs/src/06-bases.md) for the matrices behind the
structure and manual operator construction.

## How to Cite

If you use Levels.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/dnadlinger/Levels.jl/blob/master/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first take a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://dnadlinger.github.io/Levels.jl/dev/90-contributing/)
