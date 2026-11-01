# Second-order (Floquet) perturbation-theory quasienergy shifts, with per-term
# diagnostics and near-resonance warnings; the Floquet method wraps the
# nonperturbative dress_manifold engine as the cross-check.

# The perturbation-sum terms for level α of a manifold with (stripped, µs⁻¹)
# energies e driven by h1 e^{iΩt} + h1† e^{−iΩt}: named tuples
# (k, sideband, coupling, detuning, shift), all in µs⁻¹. Couplings vanishing to
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
    e = ustrip.(u"µs^-1", H0.diag)
    h1 = ustrip.(u"µs^-1", H1)
    ω = ustrip(u"µs^-1", drive_frequency)

    total = 0.0
    for term in pt_terms(e, h1, ω, α)
        if abs(term.detuning) < method.warn_ratio * term.coupling
            @warn "Near-resonant term in quasienergy perturbation sum; " *
                  "second-order perturbation theory may be inaccurate " *
                  "(cross-check with quasienergy_shift(..., Floquet()))." context =
                warn_context k = term.k sideband = term.sideband coupling =
                term.coupling * u"µs^-1" detuning = term.detuning * u"µs^-1"
        end
        total += term.shift
    end
    total * u"µs^-1"
end

"""
    quasienergy_shift_contributions(H0::Diagonal, H1::AbstractMatrix, drive_frequency, α::Integer)

Returns the per-intermediate-state breakdown of [`quasienergy_shift`](@ref): a
vector of named tuples `(k, sideband, coupling, detuning, shift)` (angular
frequency units), sorted by decreasing `|shift|`.
"""
function quasienergy_shift_contributions(
    H0::Diagonal,
    H1::AbstractMatrix,
    drive_frequency,
    α::Integer,
)
    e = ustrip.(u"µs^-1", H0.diag)
    h1 = ustrip.(u"µs^-1", H1)
    ω = ustrip(u"µs^-1", drive_frequency)
    terms = sort(pt_terms(e, h1, ω, α); by=t -> abs(t.shift), rev=true)
    [
        (;
            t.k,
            t.sideband,
            coupling=t.coupling * u"µs^-1",
            detuning=t.detuning * u"µs^-1",
            shift=t.shift * u"µs^-1",
        ) for t in terms
    ]
end

# The (H0, H1) blocks of a driven hyperfine manifold in its field eigenbasis:
# H0 the diagonal eigen-energies, H1 = Σᵢ ½ e^{iφᵢ} Xᵢ with the drive
# amplitudes rotated from the coupled basis.
function eigenbasis_blocks(m::HyperfineManifold, drives)
    isempty(drives) && throw(ArgumentError("At least one drive is required"))
    num_states = length(m.basis)
    for drive in drives
        if size(drive.amplitude) != (num_states, num_states)
            throw(
                ArgumentError(
                    "Drive amplitude matrix does not match the manifold basis size",
                ),
            )
        end
    end
    h1 = harmonic_block(drives, 1:num_states)
    Diagonal(m.energies), m.states' * h1 * m.states
end

"""
    quasienergy_shift(m::HyperfineManifold, state, drive_frequency, drives, method = PerturbationTheory())

Returns the quasienergy shift of the
adiabatically-labelled eigenstate `state` (a `StateSpec` or basis index) of the
hyperfine manifold solution `m`, for driving at `drive_frequency` (angular units).
The default method uses second-order perturbation theory; `Floquet()` performs
nonperturbative dressing.

`drives` is a vector of [`HarmonicDrive`](@ref)s at `drive_frequency` over the
manifold's coupled basis (e.g. from [`zeeman_drives`](@ref)), rotated into the
field eigenbasis internally. [`ac_zeeman_shift`](@ref) accepts a magnetic field
phasor directly.
"""
function quasienergy_shift(
    m::HyperfineManifold,
    state,
    drive_frequency,
    drives::AbstractVector{<:HarmonicDrive},
    method::QuasienergyShiftMethod=PerturbationTheory(),
)
    h0, h1 = eigenbasis_blocks(m, drives)
    α = manifold_state_index(m, state)
    quasienergy_shift(
        h0,
        h1,
        drive_frequency,
        α,
        method;
        warn_context=(; manifold=m.level, state=m.basis[α]),
    )
end

"""
    quasienergy_shift_contributions(m::HyperfineManifold, state, drive_frequency, drives)

Returns the per-intermediate-state breakdown of the manifold form of
[`quasienergy_shift`](@ref), with each entry carrying the adiabatic
[`StateSpec`](@ref) label of its intermediate state:
`(state, sideband, coupling, detuning, shift)`, sorted by decreasing `|shift|`.
"""
function quasienergy_shift_contributions(
    m::HyperfineManifold,
    state,
    drive_frequency,
    drives::AbstractVector{<:HarmonicDrive},
)
    h0, h1 = eigenbasis_blocks(m, drives)
    α = manifold_state_index(m, state)
    [
        (; state=m.basis[t.k], t.sideband, t.coupling, t.detuning, t.shift) for
        t in quasienergy_shift_contributions(h0, h1, drive_frequency, α)
    ]
end

manifold_state_index(m::HyperfineManifold, state::Integer) = Int(state)
manifold_state_index(m::HyperfineManifold, state::StateSpec) =
    stateindex(m.basis, state)

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
    uconvert(u"µs^-1", ε[α] - H0.diag[α])
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
    ac_zeeman_shift(m::HyperfineManifold, state, drive_frequency, b, method = PerturbationTheory())

Returns the ac Zeeman shift of the adiabatically-labelled eigenstate `state` of
the hyperfine manifold solution `m` under a magnetic field with Cartesian
amplitude phasor `b`, ``\\vec{B}(t) = \\mathrm{Re}[\\vec{b} e^{-iΩt}]``, where
`drive_frequency` is the angular frequency ``Ω``. Real phasors describe linearly
polarised fields; complex phasors also describe elliptical fields.

Constructs the magnetic [`HarmonicDrive`](@ref)s with [`zeeman_drives`](@ref) and
forwards to [`quasienergy_shift`](@ref), with the same method choice.
"""
function ac_zeeman_shift(
    m::HyperfineManifold,
    state,
    drive_frequency,
    b::AbstractVector{<:Quantity},
    method::QuasienergyShiftMethod=PerturbationTheory(),
)
    quasienergy_shift(
        m,
        state,
        drive_frequency,
        zeeman_drives(m.species, m.basis, b),
        method,
    )
end

"""
    ac_zeeman_shift(species, lower::StateSpec => upper::StateSpec, drive_frequency, b,
                    method = PerturbationTheory(); B)

Returns the ac Zeeman shift of the transition between the two labelled states at
the static field `B` along ẑ, under a magnetic field with Cartesian amplitude
phasor `b`, as a named tuple `(; shift, lower, upper)` with `shift = upper - lower`
(positive: the transition frequency increases). The two states may lie in the
same fine-structure manifold (e.g. a microwave qubit) or in two different ones.
Hyperfine states denote the adiabatically-labelled eigenstates.
"""
function ac_zeeman_shift(
    species::HyperfineOneElectronSpecies,
    probe::Pair,
    drive_frequency,
    b::AbstractVector{<:Quantity},
    method::QuasienergyShiftMethod=PerturbationTheory();
    B,
)
    lower, upper = probe
    fs_lower = fine_structure(parse_level(lower.level))
    fs_upper = fine_structure(parse_level(upper.level))
    m_lower = hyperfine_manifold(species, fs_lower; B)
    m_upper = fs_upper == fs_lower ? m_lower : hyperfine_manifold(species, fs_upper; B)
    δ_lower = ac_zeeman_shift(m_lower, lower, drive_frequency, b, method)
    δ_upper = ac_zeeman_shift(m_upper, upper, drive_frequency, b, method)
    (shift=δ_upper - δ_lower, lower=δ_lower, upper=δ_upper)
end

function ac_zeeman_shift(
    species::NoHyperfineOneElectronSpecies,
    probe::Pair,
    drive_frequency,
    b::AbstractVector{<:Quantity},
    method::QuasienergyShiftMethod=PerturbationTheory();
    B,
)
    function shift(state)
        basis = StateBasis(species, state.level)
        H0 = Diagonal([zeeman_shift(species, s; B) for s in basis])
        H1 = harmonic_block(zeeman_drives(species, basis, b), 1:length(basis))
        quasienergy_shift(
            H0,
            H1,
            drive_frequency,
            stateindex(basis, state),
            method;
            warn_context=(; state),
        )
    end
    lower, upper = shift(probe.first), shift(probe.second)
    (; shift=upper - lower, lower, upper)
end

export quasienergy_shift, quasienergy_shift_contributions, ac_zeeman_shift
