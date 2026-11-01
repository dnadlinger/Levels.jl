# Second-order (Floquet) perturbation-theory quasienergy shifts, with per-term
# diagnostics and near-resonance warnings; the Floquet method wraps the
# nonperturbative dress_manifold engine as the cross-check.

# The perturbation-sum terms for level α of a manifold with (stripped, rad µs⁻¹)
# energies e driven by h1 e^{iΩt} + h1† e^{−iΩt}: named tuples
# (k, sideband, coupling, detuning, shift), all in rad µs⁻¹. Couplings vanishing to
# rounding (selection rules) are skipped.
function pt_terms(e, h1, ω, α)
    wmax = maximum(abs, h1)
    terms = @NamedTuple{
        k::Int,
        sideband::Int,
        coupling::Float64,
        detuning::Float64,
        shift::Float64,
    }[]
    for k in eachindex(e)
        k == α && continue
        Δ = e[α] - e[k]
        for (coupling, sideband) in ((abs(h1[k, α]), -1), (abs(h1[α, k]), +1))
            coupling <= 1e-12 * wmax && continue
            detuning = Δ + sideband * ω
            push!(
                terms,
                (; k, sideband, coupling, detuning, shift=coupling^2 / detuning),
            )
        end
    end
    terms
end

"""
    quasienergy_shift(H0::Diagonal, H1::AbstractMatrix, drive_frequency, α::Integer,
                    method = PerturbationTheory())

Returns the quasienergy shift relative to the static energy of state `α`
driven by ``H(t) = H_0 + H_1 e^{iΩt} + H_1^† e^{−iΩt}`` (the convention of
[`dress_manifold`](@ref)), from second-order Floquet perturbation theory:

```math
δ_α = \\sum_{k ≠ α} \\frac{|H_{1,kα}|^2}{E_α - E_k - Ω}
    + \\frac{|H^†_{1,kα}|^2}{E_α - E_k + Ω}.
```

All quantities are in angular frequency units. A warning is emitted whenever a
denominator is smaller than `method.warn_ratio` times the corresponding coupling
matrix element (near-resonant drive, where the perturbative treatment degrades
— cross-check such cases with `quasienergy_shift(..., Floquet())`).
"""
function quasienergy_shift(
    H0::Diagonal,
    H1::AbstractMatrix,
    drive_frequency,
    α::Integer,
    method::PerturbationTheory=PerturbationTheory();
    warn_context=(;),
)
    e = ustrip.(u"rad/µs", H0.diag)
    h1 = ustrip.(u"rad/µs", H1)
    ω = ustrip(u"rad/µs", drive_frequency)

    total = 0.0
    for term in pt_terms(e, h1, ω, α)
        if abs(term.detuning) < method.warn_ratio * term.coupling
            @warn "Near-resonant term in quasienergy perturbation sum; " *
                  "second-order perturbation theory may be inaccurate " *
                  "(cross-check with quasienergy_shift(..., Floquet()))." context =
                warn_context k = term.k sideband = term.sideband coupling =
                term.coupling * u"rad/µs" detuning = term.detuning * u"rad/µs"
        end
        total += term.shift
    end
    total * u"rad/µs"
end

"""
    quasienergy_shift_contributions(H0::Diagonal, H1::AbstractMatrix, drive_frequency, α::Integer)

Returns the per-intermediate-state breakdown of [`quasienergy_shift`](@ref): a
vector of named tuples `(k, sideband, coupling, detuning, shift)`, sorted by
decreasing `|shift|`. `k` is a basis index, `sideband` is ±1, and `coupling`,
`detuning` and `shift` are angular frequencies.
"""
function quasienergy_shift_contributions(
    H0::Diagonal,
    H1::AbstractMatrix,
    drive_frequency,
    α::Integer,
)
    e = ustrip.(u"rad/µs", H0.diag)
    h1 = ustrip.(u"rad/µs", H1)
    ω = ustrip(u"rad/µs", drive_frequency)
    terms = sort(pt_terms(e, h1, ω, α); by=t -> abs(t.shift), rev=true)
    [
        (;
            t.k,
            t.sideband,
            coupling=t.coupling * u"rad/µs",
            detuning=t.detuning * u"rad/µs",
            shift=t.shift * u"rad/µs",
        ) for t in terms
    ]
end

# The manifold index range of a state of the structure, requiring the whole
# manifold to be present: the perturbation sum runs over all of its states.
function manifold_range(s::Eigenstructure, state::StateSpec)
    fs = fine_structure(validate_state_kind(s.species, state).level)
    if !iscomplete(s, fs)
        throw(
            ArgumentError(
                "The structure must contain the complete '$fs' manifold: " *
                "quasienergy shifts include couplings to all of its states",
            ),
        )
    end
    staterange(s.basis, fs)
end

manifold_range(s::Eigenstructure, state::Integer) = manifold_range(s, s.basis[state])

state_index(s::Eigenstructure, state::Integer) = Int(state)
state_index(s::Eigenstructure, state::StateSpec) = stateindex(s.basis, state)

# Validates drives before restricting them to individual manifolds, so that
# couplings between manifolds cannot be silently discarded.
function validate_shift_drives(
    s::Eigenstructure,
    drives::AbstractVector{<:HarmonicDrive},
)
    isempty(drives) && throw(ArgumentError("At least one drive is required"))
    num_states = length(s.basis)
    manifolds = [fine_structure(state.level) for state in s.basis]
    for drive in drives
        if size(drive.amplitude) != (num_states, num_states)
            throw(
                ArgumentError(
                    "Drive amplitude matrix does not match the structure's basis size",
                ),
            )
        end
        if !all(
            iszero(drive.amplitude[i, k]) for i in 1:num_states for
            k in 1:num_states if manifolds[i] != manifolds[k]
        )
            throw(
                ArgumentError(
                    "Harmonic drives must be block-diagonal in the fine-structure " *
                    "manifolds; use quasienergy_shift(H0, H1, ...) for a model " *
                    "with couplings between manifolds",
                ),
            )
        end
    end
    drives
end

# The (H0, H1) blocks of a driven manifold of the structure in its field
# eigenbasis: H0 the diagonal eigen-energies, H1 = Σᵢ ½ e^{iφᵢ} Xᵢ.
function eigenbasis_blocks(s::Eigenstructure, state, drives)
    r = manifold_range(s, state)
    α = state_index(s, state) - first(r) + 1
    Diagonal(offsets(s)[r]), harmonic_block(validate_shift_drives(s, drives), r), α
end

"""
    quasienergy_shift(s::Eigenstructure, state, drive_frequency, drives, method = PerturbationTheory())

Returns the quasienergy shift of the labelled eigenstate `state` (a `StateSpec`,
or an index into `basis(s)`) of the structure `s` — which fixes the species,
the static field and the basis — for driving at `drive_frequency` (angular
units), as an angular frequency. The default method uses second-order
perturbation theory; `Floquet()` performs nonperturbative dressing.

`drives` is a vector of [`HarmonicDrive`](@ref)s at the same `drive_frequency`,
over `basis(s)` in the structure's field eigenbasis. Every drive must be
block-diagonal in the fine-structure manifolds; couplings between different
``F`` levels of one manifold are allowed. The state's whole fine-structure
manifold must be present and contiguous in the basis
(cf. [`Levels.iscomplete`](@ref)). The plain-matrix form supports models with
couplings between manifolds when the caller supplies their full static
energies and drive matrix.

Magnetic drives may be constructed with [`zeeman_drives`](@ref)`(s, b)`;
[`ac_zeeman_shift`](@ref) accepts the field phasor directly. Works for both species kinds:
on the exact hyperfine + Zeeman eigenstates of a hyperfine species, on the
Zeeman-shifted ``|J, m_J⟩`` states of a fine-structure one.
"""
function quasienergy_shift(
    s::Eigenstructure,
    state::Union{StateSpec,Integer},
    drive_frequency,
    drives::AbstractVector{<:HarmonicDrive},
    method::QuasienergyShiftMethod=PerturbationTheory(),
)
    h0, h1, α = eigenbasis_blocks(s, state, drives)
    quasienergy_shift(
        h0,
        h1,
        drive_frequency,
        α,
        method;
        warn_context=(; state=s.basis[state_index(s, state)]),
    )
end

"""
    quasienergy_shift_contributions(s::Eigenstructure, state, drive_frequency, drives)

Returns the per-intermediate-state breakdown of the structure form of
[`quasienergy_shift`](@ref), with each entry carrying the adiabatic
[`StateSpec`](@ref) label of its intermediate state:
`(state, sideband, coupling, detuning, shift)`, sorted by decreasing `|shift|`.
The drive and manifold restrictions are the same as for [`quasienergy_shift`](@ref).
"""
function quasienergy_shift_contributions(
    s::Eigenstructure,
    state::Union{StateSpec,Integer},
    drive_frequency,
    drives::AbstractVector{<:HarmonicDrive},
)
    h0, h1, α = eigenbasis_blocks(s, state, drives)
    r = manifold_range(s, state)
    [
        (; state=s.basis[first(r)+t.k-1], t.sideband, t.coupling, t.detuning, t.shift)
        for t in quasienergy_shift_contributions(h0, h1, drive_frequency, α)
    ]
end

"""
    quasienergy_shift(s::Eigenstructure, lower::StateSpec => upper::StateSpec,
                    drive_frequency, drives, method = PerturbationTheory())

Returns the differential quasienergy shift of the transition between the two labelled
eigenstates of the structure `s` at its static field, for driving at
`drive_frequency` (angular units), as a named tuple `(; shift, lower, upper)`
with `shift = upper - lower` (positive: the transition frequency increases).

`drives` is a vector of [`HarmonicDrive`](@ref)s over `basis(s)`, at the same
`drive_frequency` and block-diagonal in the fine-structure manifolds.
The two probed states may lie in the same fine-structure manifold
(e.g. a microwave qubit) or in two different ones; their manifolds must be
complete and contiguous in `s`.
"""
function quasienergy_shift(
    s::Eigenstructure,
    probe::Pair{<:StateSpec,<:StateSpec},
    drive_frequency,
    drives::AbstractVector{<:HarmonicDrive},
    method::QuasienergyShiftMethod=PerturbationTheory(),
)
    δ_lower = quasienergy_shift(s, probe.first, drive_frequency, drives, method)
    δ_upper = quasienergy_shift(s, probe.second, drive_frequency, drives, method)
    (shift=δ_upper - δ_lower, lower=δ_lower, upper=δ_upper)
end

"""
    quasienergy_shift(H0::Diagonal, H1::AbstractMatrix, drive_frequency, α::Integer, method::Floquet)

Returns the quasienergy shift of state `α` from the nonperturbative Floquet dressing
([`dress_manifold`](@ref), truncated at `harmonic_cutoff` harmonics): the difference
between the dressed quasienergy and the static eigen-energy.

Exact in the drive amplitude up to the harmonic truncation, so it serves as the
cross-check for [`quasienergy_shift`](@ref) near resonances; a warning is emitted
when the dressed-state identification becomes ambiguous (harmonic-0 overlap
below `method.min_overlap`).
"""
function quasienergy_shift(
    H0::Diagonal,
    H1::AbstractMatrix,
    drive_frequency,
    α::Integer,
    method::Floquet;
    warn_context=(;),
)
    ε, _ = dress_manifold(
        H0,
        H1,
        drive_frequency,
        method.harmonic_cutoff,
        1;
        min_overlap=method.min_overlap,
    )
    uconvert(u"rad/µs", ε[α] - H0.diag[α])
end

"""
    quasienergy_shift(dt::DrivenTransition, i::Integer, method = PerturbationTheory())

Returns the quasienergy shift of basis state `i`
of the driven-transition model from its harmonic drives, over the manifold
containing `i` (cf. the manifold blocks used by [`dress`](@ref)).
"""
function quasienergy_shift(
    dt::DrivenTransition,
    i::Integer,
    method::QuasienergyShiftMethod=PerturbationTheory(),
)
    range = i in dt.lower_range ? dt.lower_range : dt.upper_range
    h0, h1 = manifold_blocks(dt, range)
    quasienergy_shift(
        h0,
        h1,
        dt.drive_frequency,
        i - first(range) + 1,
        method;
        warn_context=(; state=dt.basis[i]),
    )
end

"""
    ac_zeeman_shift(s::Eigenstructure, state, drive_frequency, b, method = PerturbationTheory())
    ac_zeeman_shift(s::Eigenstructure, lower::StateSpec => upper::StateSpec,
                    drive_frequency, b, method = PerturbationTheory())

Returns the ac Zeeman shift under a magnetic field with Cartesian amplitude
phasor `b`, ``\\vec{B}(t) = \\mathrm{Re}[\\vec{b} e^{-iΩt}]``, where
`drive_frequency` is the angular frequency ``Ω``. Real phasors describe
linearly polarised fields; complex phasors also describe elliptical fields.

Constructs the magnetic [`HarmonicDrive`](@ref)s with [`zeeman_drives`](@ref)`(s, b)`
and forwards to [`quasienergy_shift`](@ref), with the same method choice and
manifold requirements. `state` is a `StateSpec` or an index into `basis(s)`.
The state form returns an angular frequency; the pair form returns
`(; shift, lower, upper)` with `shift = upper - lower`.
"""
function ac_zeeman_shift(
    s::Eigenstructure,
    probe::Union{StateSpec,Integer,Pair{<:StateSpec,<:StateSpec}},
    drive_frequency,
    b::AbstractVector{<:Quantity},
    method::QuasienergyShiftMethod=PerturbationTheory(),
)
    quasienergy_shift(s, probe, drive_frequency, zeeman_drives(s, b), method)
end

export quasienergy_shift, quasienergy_shift_contributions, ac_zeeman_shift
