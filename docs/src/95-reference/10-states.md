# [Levels and states](@id reference-states)

Types identifying an electronic level within a species, and the parser for
spectroscopic notation.

`NoHyperfineNumberSpec` is the canonical fine-structure form, and
`HyperfineNumberSpec` its hyperfine counterpart carrying an additional ``F``
quantum number; `SpectroscopicSpec` wraps a notation string such as `"D_5/2"`
or `"S_1/2 F=4"` (the `F=` suffix selecting the hyperfine form). Since strings
convert to level specifications implicitly, the APIs in the rest of the package
accept level strings directly.

```@autodocs
Modules = [Levels]
Pages = ["states.jl"]
```
