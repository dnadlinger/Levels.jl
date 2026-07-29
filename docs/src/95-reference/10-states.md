# [Levels and states](@id reference-states)

Types identifying an electronic level within a species, and the parser for
spectroscopic notation.

`NoHyperfineNumberSpec` is the canonical form; `SpectroscopicSpec` wraps a notation
string such as `"D_5/2"`. Since strings convert to level specifications implicitly,
the APIs in the rest of the package accept level strings directly.

```@autodocs
Modules = [Levels]
Pages = ["states.jl"]
```
