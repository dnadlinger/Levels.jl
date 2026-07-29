using WignerSymbols

"""
    beam_vectors(φ_k, γ_pol, η_pol = 0.0) -> (n, ε)

Returns the unit propagation direction `n` and (complex) polarisation vector `ε`
for a beam at angle `φ_k` to the quantisation axis ẑ, with k in the x–z plane.

`γ_pol` rotates the polarisation from the in-(k, z)-plane vector towards ŷ;
`η_pol` is the relative phase of the two components (`0`: linear polarisation,
`±π/2` at `γ_pol = π/4`: circular).
"""
function beam_vectors(φ_k, γ_pol, η_pol=0.0)
    n = [sin(φ_k), 0.0, cos(φ_k)]
    ε =
        cos(γ_pol) .* [complex(cos(φ_k)), 0.0, complex(-sin(φ_k))] .+
        (cis(η_pol) * sin(γ_pol)) .* [0.0, 1.0, 0.0]
    n, ε
end

"""
    quadrupole_geometry(ε, n)

Returns the rank-2 spherical tensor components ``Γ_q = (ε ⊗ n)^{(2)}_q`` of the
outer product of polarisation and beam direction — the geometric factor of the
electric-quadrupole coupling for ``Δm = q`` — as a length-5 complex vector
indexed by `q + 3` for ``q = -2, …, 2``.

The components are kept complex so that the relative phases between the ``Δm``
channels are preserved.
"""
function quadrupole_geometry(ε, n)
    Γ = zeros(ComplexF64, 5)
    for q in -2:2, μ in -1:1
        ν = q - μ
        abs(ν) <= 1 || continue
        Γ[q+3] +=
            clebschgordan(1, μ, 1, ν, 2, q) *
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
transition to a known carrier Rabi frequency using [`rabi_normalised`](@ref).
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
    rabi_normalised(couplings, basis::StateBasis, transition::Pair, Ω0)

Scales the given relative coupling matrix such that the entry for the
`lower => upper` `transition` has magnitude `Ω0` (its carrier Rabi frequency, in
angular units), fixing the physical scale of all the couplings.

An error is raised if the amplitude for the given transition (nearly) vanishes in
`couplings`, as the normalisation would then be ill-defined.
"""
function rabi_normalised(
    couplings::AbstractMatrix,
    basis::StateBasis,
    transition::Pair,
    Ω0,
)
    c = couplings[
        stateindex(basis, transition.second),
        stateindex(basis, transition.first),
    ]
    if abs(c) <= 1e-9 * maximum(abs, couplings)
        throw(
            ArgumentError(
                "Coupling for transition $(transition.first) => $(transition.second) " *
                "(nearly) vanishes for the given geometry",
            ),
        )
    end
    (Ω0 / abs(c)) .* couplings
end

export beam_vectors, quadrupole_geometry, quadrupole_couplings, rabi_normalised
