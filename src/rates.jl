using Unitful
using WignerSymbols

"""
    multipole_rank(lower, upper)

Returns the electric-multipole rank of the transition between the two given
levels — 1 (E1) if the level parities differ, 2 (E2) if they are equal — as
determined from the orbital angular momenta. Magnetic multipoles and electric
orders beyond E2 are not considered.
"""
function multipole_rank(lower, upper)
    lo = convert(NoHyperfineNumberSpec, lower)
    hi = convert(NoHyperfineNumberSpec, upper)
    isodd(lo.l + hi.l) ? 1 : 2
end

"""
    clebsch_gordan(lower::StateSpec, upper::StateSpec)

Returns the Clebsch–Gordan coefficient ``⟨j m; R Δm | j' m'⟩`` giving the
amplitude of the transition between the two given states relative to the reduced
matrix element of the transition, with the rank ``R`` the electric-multipole
order from [`multipole_rank`](@ref).
"""
function clebsch_gordan(lower::StateSpec, upper::StateSpec)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    rank = multipole_rank(lo, hi)
    Δm = upper.m - lower.m
    abs(Δm) <= rank || return 0.0
    Float64(clebschgordan(lo.j, lower.m, rank, Δm, hi.j, upper.m))
end

"""
    rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity, ε, n)

Returns the Rabi frequency (in angular units, such that the excitation
probability oscillates as ``\\sin^2(Ω t / 2)`` on resonance) for driving the
given transition with a running wave of the given intensity, polarisation `ε`
and propagation direction `n` (cf. [`beam_vectors`](@ref); `n` only enters for
quadrupole transitions).

The coupling strength is obtained from the Einstein A coefficient of the
transition following `[James1998]`,

``Ω = \\sqrt{6 π c^2 I A / (ħ ω^3)} \\, |⟨j m; 1 Δm | j' m'⟩ \\, d_{Δm}|``

for electric-dipole and

``Ω = \\sqrt{20 π c^2 I A / (ħ ω^3)} \\, |⟨j m; 2 Δm | j' m'⟩ \\, Γ_{Δm}|``

for electric-quadrupole transitions, with ``d_q`` and ``Γ_q`` the geometric
channel amplitudes from [`dipole_geometry`](@ref) and
[`quadrupole_geometry`](@ref).

# References

- `[James1998]`: D. F. V. James, "Quantum dynamics of cold trapped ions with
  application to quantum computation", Appl. Phys. B **66**, 181 (1998),
  [doi:10.1007/s003400050373](https://doi.org/10.1007/s003400050373); Eqs.
  (5.11)–(5.13) in the
  [arXiv:quant-ph/9702053](https://arxiv.org/abs/quant-ph/9702053) numbering.
"""
function rabi_frequency(species, lower::StateSpec, upper::StateSpec, intensity, ε, n)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    a = einstein_a(species, lo, hi)
    isnothing(a) && throw(ArgumentError("No known transition between '$lo' and '$hi'"))
    ω = transition_frequency(species, lo, hi)
    rank = multipole_rank(lo, hi)
    Δm = upper.m - lower.m
    prefactor, geometry = if abs(Δm) > rank
        0.0, 0.0im
    elseif rank == 1
        6.0, dipole_geometry(ε)[Int(Δm)+2]
    else
        20.0, quadrupole_geometry(ε, n)[Int(Δm)+3]
    end
    scale = prefactor * π * u"c"^2 * intensity * a / (u"ħ" * ω^3)
    uconvert(u"µs^-1", sqrt(scale) * abs(clebsch_gordan(lower, upper) * geometry))
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

export rabi_frequency, rabi_normalised
public clebsch_gordan, multipole_rank
