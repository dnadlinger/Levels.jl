using StaticArrays: SVector

"""
    beam_vectors(φ_k, γ_pol, η_pol = 0.0) -> (direction, ε)

Returns the unit propagation direction `direction` and (complex) polarisation vector `ε`
(as `SVector{3}`s) for a beam at angle `φ_k` to the quantisation axis ẑ, with k
in the x–z plane.

`γ_pol` rotates the polarisation from the in-(k, z)-plane vector towards ŷ;
`η_pol` is the relative phase of the two components (`0`: linear polarisation,
`±π/2` at `γ_pol = π/4`: circular). The physical field is
``E(t) = \\mathrm{Re}(ε E_0 e^{-i ω t})``, so `beam_vectors(0.0, π/4, π/2)`
gives ``ε = (x̂ + i ŷ)/\\sqrt{2}``, a σ⁺ beam along the quantisation axis
(driving ``Δm = +1``).
"""
function beam_vectors(φ_k, γ_pol, η_pol=0.0)
    direction = SVector{3,Float64}(sin(φ_k), 0.0, cos(φ_k))
    ε =
        cos(γ_pol) .* SVector{3,ComplexF64}(cos(φ_k), 0.0, -sin(φ_k)) .+
        (cis(η_pol) * sin(γ_pol)) .* SVector{3,ComplexF64}(0.0, 1.0, 0.0)
    direction, ε
end

"""
    dipole_geometry(ε)

Returns the geometric amplitudes ``d_q`` of the electric-dipole coupling for the
``Δm = q`` transition channels for polarisation `ε` — the coefficients of the
atomic operators ``r C^{(1)}_q`` in ``ε · r``, i.e. ``d_q = (-1)^q ε_{-q}`` in
terms of the (unconjugated) spherical components of the polarisation — as a
length-3 complex `SVector` indexed by `q + 2` for ``q = -1, 0, 1``.

A σ⁺-polarised beam, ``ε = (x̂ + i ŷ)/\\sqrt{2}``, thus carries its full weight
in the ``Δm = +1`` channel. The components are kept complex so that the relative
phases between the ``Δm`` channels are preserved.
"""
dipole_geometry(ε) = SVector{3,ComplexF64}(
    -spherical_component(ε, 1),
    spherical_component(ε, 0),
    -spherical_component(ε, -1),
)

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
    quadrupole_geometry(ε, direction)

Returns the geometric amplitudes ``Γ_q`` of the electric-quadrupole coupling for
the ``Δm = q`` transition channels for polarisation `ε` and beam direction `direction` —
the coefficients of the atomic operators ``r^2 (r̂ ⊗ r̂)^{(2)}_q`` in
``(ε · r)(direction · r)``, i.e. ``Γ_q = (-1)^q (ε ⊗ direction)^{(2)}_{-q}`` in terms of the
irreducible tensor product of polarisation and beam direction — as a length-5
complex `SVector` indexed by `q + 3` for ``q = -2, …, 2``.

A σ⁺-polarised beam along the quantisation axis thus carries its full weight in
the ``Δm = +1`` channel. The components are kept complex so that the relative
phases between the ``Δm`` channels are preserved.
"""
function quadrupole_geometry(ε, direction)
    Γ = zeros(ComplexF64, 5)
    for q in -2:2, μ in -1:1
        ν = -q - μ
        abs(ν) <= 1 || continue
        Γ[q+3] +=
            (-1.0)^q *
            RANK2_CG[μ+2, ν+2] *
            spherical_component(ε, μ) *
            spherical_component(direction, ν)
    end
    SVector{5,ComplexF64}(Γ)
end

export beam_vectors, dipole_geometry, quadrupole_geometry
