using WignerSymbols

"""
    beam_vectors(φ_k, γ_pol, η_pol = 0.0) -> (n, ε)

Returns the unit propagation direction `n` and (complex) polarisation vector `ε`
for a beam at angle `φ_k` to the quantisation axis ẑ, with k in the x–z plane.

`γ_pol` rotates the polarisation from the in-(k, z)-plane vector towards ŷ;
`η_pol` is the relative phase of the two components (`0`: linear polarisation,
`±π/2` at `γ_pol = π/4`: circular). The physical field is
``E(t) = \\mathrm{Re}(ε E_0 e^{-i ω t})``, so `beam_vectors(0.0, π/4, π/2)`
gives ``ε = (x̂ + i ŷ)/\\sqrt{2}``, a σ⁺ beam along the quantisation axis
(driving ``Δm = +1``).
"""
function beam_vectors(φ_k, γ_pol, η_pol=0.0)
    n = [sin(φ_k), 0.0, cos(φ_k)]
    ε =
        cos(γ_pol) .* [complex(cos(φ_k)), 0.0, complex(-sin(φ_k))] .+
        (cis(η_pol) * sin(γ_pol)) .* [0.0, 1.0, 0.0]
    n, ε
end

"""
    dipole_geometry(ε)

Returns the geometric amplitudes ``d_q`` of the electric-dipole coupling for the
``Δm = q`` transition channels for polarisation `ε` — the coefficients of the
atomic operators ``r C^{(1)}_q`` in ``ε · r``, i.e. ``d_q = (-1)^q ε_{-q}`` in
terms of the (unconjugated) spherical components of the polarisation — as a
length-3 complex vector indexed by `q + 2` for ``q = -1, 0, 1``.

A σ⁺-polarised beam, ``ε = (x̂ + i ŷ)/\\sqrt{2}``, thus carries its full weight
in the ``Δm = +1`` channel. The components are kept complex so that the relative
phases between the ``Δm`` channels are preserved.
"""
dipole_geometry(ε) = [(-1.0)^q * spherical_component(ε, -q) for q in -1:1]

"""
The Clebsch–Gordan coefficients ``⟨1 μ; 1 ν | 2 (μ + ν)⟩`` coupling two rank-1
tensors to rank two, indexed by `[μ + 2, ν + 2]`.

These nine fixed numbers are tabulated rather than requested from WignerSymbols
on each use, as [`quadrupole_geometry`](@ref) sits in the hot path of
Rabi-frequency and light-shift evaluations.
"""
const RANK2_CG = [
    1.0 sqrt(1 / 2) sqrt(1 / 6)
    sqrt(1 / 2) sqrt(2 / 3) sqrt(1 / 2)
    sqrt(1 / 6) sqrt(1 / 2) 1.0
]

"""
    quadrupole_geometry(ε, n)

Returns the geometric amplitudes ``Γ_q`` of the electric-quadrupole coupling for
the ``Δm = q`` transition channels for polarisation `ε` and beam direction `n` —
the coefficients of the atomic operators ``r^2 (r̂ ⊗ r̂)^{(2)}_q`` in
``(ε · r)(n · r)``, i.e. ``Γ_q = (-1)^q (ε ⊗ n)^{(2)}_{-q}`` in terms of the
irreducible tensor product of polarisation and beam direction — as a length-5
complex vector indexed by `q + 3` for ``q = -2, …, 2``.

A σ⁺-polarised beam along the quantisation axis thus carries its full weight in
the ``Δm = +1`` channel. The components are kept complex so that the relative
phases between the ``Δm`` channels are preserved.
"""
function quadrupole_geometry(ε, n)
    Γ = zeros(ComplexF64, 5)
    for q in -2:2, μ in -1:1
        ν = -q - μ
        abs(ν) <= 1 || continue
        Γ[q+3] +=
            (-1.0)^q *
            RANK2_CG[μ+2, ν+2] *
            spherical_component(ε, μ) *
            spherical_component(n, ν)
    end
    Γ
end

"""
    quadrupole_couplings(basis::StateBasis, lower_level, upper_level, ε, n)

Returns the relative electric-quadrupole coupling amplitudes
``c_{m→m'} = ⟨j_l m; 2 Δm | j_u m'⟩ Γ_{Δm}`` between the Zeeman states of the two
given levels for the given beam geometry (cf. [`quadrupole_geometry`](@ref)), as
a complex matrix over the given basis with entries in the upper⟨row|lower⟩⟨col|
block (zero elsewhere).

The overall scale is arbitrary; it is typically fixed by normalising one probed
transition to a known carrier Rabi frequency using [`rabi_normalised`](@ref)
(e.g. one computed from the beam intensity with [`rabi_frequency`](@ref)).
"""
function quadrupole_couplings(basis::StateBasis, lower_level, upper_level, ε, n)
    lower = convert(NoHyperfineNumberSpec, lower_level)
    upper = convert(NoHyperfineNumberSpec, upper_level)
    Γ = quadrupole_geometry(ε, n)
    C = zeros(ComplexF64, length(basis), length(basis))
    for (i, lower_state) in enumerate(basis)
        lower_state.level == lower || continue
        for (k, upper_state) in enumerate(basis)
            upper_state.level == upper || continue
            q = upper_state.m - lower_state.m
            abs(q) <= 2 || continue
            C[k, i] =
                clebschgordan(lower.j, lower_state.m, 2, q, upper.j, upper_state.m) *
                Γ[Int(q)+3]
        end
    end
    C
end

"""
    quadrupole_couplings(species, basis::StateBasis, lower, upper, ε, n)

Species-first form: for a [`HyperfineOneElectronSpecies`](@ref), returns the
relative ``F``-resolved coupling amplitudes
``c = ⟨F m; 2 Δm | F' m'⟩ β^{(2)}(F → F') Γ_{Δm}`` (cf.
[`Levels.hyperfine_reduction`](@ref)) over a hyperfine basis. `lower` and
`upper` may each be a single hyperfine ``F`` level or a fine-structure
manifold, in which case all its ``F`` levels present in the basis contribute.

These are zero-field amplitudes; at finite field the ``F`` mixing within the
manifolds modifies them. The exact at-field matrix over a complete-manifold
basis is `V' * C * V` with `V` from [`eigenbasis_transform`](@ref) (the
rotation the hyperfine `Levels.PeriodicDynamics.DrivenTransition` applies
internally — so pass *this* zero-field matrix there, never a pre-rotated one);
individual at-field components are available via
[`transition_amplitude`](@ref).

For a [`NoHyperfineOneElectronSpecies`](@ref) this simply forwards to the
species-less form.
"""
function quadrupole_couplings(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    lower_level,
    upper_level,
    ε,
    n,
)
    lower = hyperfine_level_list(species, lower_level)
    upper = hyperfine_level_list(species, upper_level)
    Γ = quadrupole_geometry(ε, n)
    C = zeros(ComplexF64, length(basis), length(basis))
    for lo in lower, hi in upper
        β = hyperfine_reduction(species.nuclear_spin, lo, hi; rank=2)
        iszero(β) && continue
        for (i, lower_state) in enumerate(basis)
            lower_state.level == lo || continue
            for (k, upper_state) in enumerate(basis)
                upper_state.level == hi || continue
                q = upper_state.m - lower_state.m
                abs(q) <= 2 || continue
                C[k, i] =
                    clebschgordan(lo.f, lower_state.m, 2, q, hi.f, upper_state.m) *
                    β *
                    Γ[Int(q)+3]
            end
        end
    end
    C
end

quadrupole_couplings(
    species::NoHyperfineOneElectronSpecies,
    basis::StateBasis{NoHyperfineNumberSpec},
    lower_level,
    upper_level,
    ε,
    n,
) = quadrupole_couplings(basis, lower_level, upper_level, ε, n)

export beam_vectors, dipole_geometry, quadrupole_geometry, quadrupole_couplings
