# [Hyperfine structure](@id reference-hyperfine)

For a species with nuclear spin ``I``, each fine-structure level of angular
momentum ``J`` splits into levels ``F = |I-J|, …, I+J``. The
[`HyperfineConstants`](@ref) determine their zero-field energies through
[`hyperfine_shift`](@ref). `StateBasis(species, "S_1/2", "D_5/2")` expands
both manifolds into all their ``|F, m_F⟩`` states.

The nuclear g-factor is the signed value ``g_I = μ_I/(I μ_N)``, so the Zeeman
Hamiltonian is ``H_Z = B (g_J μ_B J_z - g_I μ_N I_z)``, with ``J_z`` and
``I_z`` in units of ``ħ``. The species data uses the effective nuclear moment
in the ion, without a diamagnetic-shielding correction.

At finite field, ``F`` levels mix and ``m_F`` remains exact. An
[`Eigenstructure`](@ref) solves the hyperfine + Zeeman Hamiltonian;
[`zeeman_shift`](@ref) and [`zeeman_sensitivity`](@ref) then evaluate its
labelled eigenstates. [`insensitive_field`](@ref) finds a sensitivity zero
within a supplied bracket. [Bases and representations](@ref bases) defines
the adiabatic labels, electron-first coupling order and eigenvector phases.

```@autodocs
Modules = [Levels]
Pages = ["hyperfine.jl"]
```
