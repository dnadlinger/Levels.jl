# Second-order (Floquet) perturbation-theory ac Zeeman shifts, with per-term
# diagnostics and near-resonance warnings; floquet_zeeman_shift wraps the
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
    ac_zeeman_shift(H0::Diagonal, H1::AbstractMatrix, Ω, α::Integer;
                    warn_ratio = 20)

Returns the time-averaged (quasienergy) shift of level `α` of one manifold
driven by ``H(t) = H_0 + H_1 e^{iΩt} + H_1^† e^{−iΩt}`` (the convention of
[`dress_manifold`](@ref)), from second-order Floquet perturbation theory:

```math
δ_α = \\sum_{k ≠ α} \\frac{|H_{1,kα}|^2}{E_α - E_k - Ω}
    + \\frac{|H^†_{1,kα}|^2}{E_α - E_k + Ω}.
```

All quantities are in angular frequency units. A warning is emitted whenever a
denominator is smaller than `warn_ratio` times the corresponding coupling
matrix element (near-resonant drive, where the perturbative treatment degrades
— cross-check such cases with [`floquet_zeeman_shift`](@ref)).
"""
function ac_zeeman_shift(
    H0::Diagonal,
    H1::AbstractMatrix,
    Ω,
    α::Integer;
    warn_ratio=20,
    warn_context=(;),
)
    e = ustrip.(u"µs^-1", H0.diag)
    h1 = ustrip.(u"µs^-1", H1)
    ω = ustrip(u"µs^-1", Ω)

    total = 0.0
    for term in pt_terms(e, h1, ω, α)
        if abs(term.detuning) < warn_ratio * term.coupling
            @warn "Near-resonant term in ac Zeeman perturbation sum; " *
                  "second-order perturbation theory may be inaccurate " *
                  "(cross-check with floquet_zeeman_shift)." context = warn_context k =
                term.k sideband = term.sideband coupling = term.coupling * u"µs^-1" detuning =
                term.detuning * u"µs^-1"
        end
        total += term.shift
    end
    total * u"µs^-1"
end

"""
    ac_zeeman_contributions(H0::Diagonal, H1::AbstractMatrix, Ω, α::Integer)

Returns the per-intermediate-state breakdown of [`ac_zeeman_shift`](@ref): a
vector of named tuples `(k, sideband, coupling, detuning, shift)` (angular
frequency units), sorted by decreasing `|shift|`.
"""
function ac_zeeman_contributions(H0::Diagonal, H1::AbstractMatrix, Ω, α::Integer)
    e = ustrip.(u"µs^-1", H0.diag)
    h1 = ustrip.(u"µs^-1", H1)
    ω = ustrip(u"µs^-1", Ω)
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
    n = length(m.basis)
    for drive in drives
        if size(drive.amplitude) != (n, n)
            throw(
                ArgumentError(
                    "Drive amplitude matrix does not match the manifold basis size",
                ),
            )
        end
    end
    h1 = sum(
        (0.5 * cis(drive.phase)) .* (m.states' * drive.amplitude * m.states) for
        drive in drives
    )
    Diagonal(m.energies), h1
end

# Resolves drives given either directly or as a complex Cartesian field phasor.
resolve_drives(m::HyperfineManifold, drives::AbstractVector{<:HarmonicDrive}) = drives
resolve_drives(m::HyperfineManifold, b) = zeeman_drives(m.species, m.basis, b)

"""
    ac_zeeman_shift(m::HyperfineManifold, state, Ω, drives; warn_ratio = 20)

Returns the second-order perturbation-theory ac Zeeman shift of the
adiabatically-labelled eigenstate `state` (a `StateSpec` or basis index) of the
hyperfine manifold solution `m`, for driving at the angular frequency `Ω`.

`drives` is either a vector of [`HarmonicDrive`](@ref)s over the manifold's
coupled basis (e.g. from [`zeeman_drives`](@ref)), or a complex Cartesian
field-amplitude phasor ``\\vec{b}`` with
``\\vec{B}(t) = \\mathrm{Re}[\\vec{b} e^{-iΩt}]``; either is rotated into the
field eigenbasis internally.
"""
function ac_zeeman_shift(m::HyperfineManifold, state, Ω, drives; warn_ratio=20)
    h0, h1 = eigenbasis_blocks(m, resolve_drives(m, drives))
    α = manifold_state_index(m, state)
    ac_zeeman_shift(
        h0,
        h1,
        Ω,
        α;
        warn_ratio,
        warn_context=(; manifold=m.level, state=m.basis[α]),
    )
end

"""
    ac_zeeman_contributions(m::HyperfineManifold, state, Ω, drives)

Returns the per-intermediate-state breakdown of the manifold form of
[`ac_zeeman_shift`](@ref), with each entry carrying the adiabatic
[`StateSpec`](@ref) label of its intermediate state:
`(state, sideband, coupling, detuning, shift)`, sorted by decreasing `|shift|`.
"""
function ac_zeeman_contributions(m::HyperfineManifold, state, Ω, drives)
    h0, h1 = eigenbasis_blocks(m, resolve_drives(m, drives))
    α = manifold_state_index(m, state)
    [
        (; state=m.basis[t.k], t.sideband, t.coupling, t.detuning, t.shift) for
        t in ac_zeeman_contributions(h0, h1, Ω, α)
    ]
end

manifold_state_index(m::HyperfineManifold, state::Integer) = Int(state)
manifold_state_index(m::HyperfineManifold, state::StateSpec) =
    stateindex(m.basis, state)

"""
    ac_zeeman_shift(species::HyperfineOneElectronSpecies, probe::Pair,
                    static_field, Ω, drives; warn_ratio = 20)

Returns the ac Zeeman shift of the transition between the two
adiabatically-labelled hyperfine eigenstates `probe` (`lower => upper`
[`StateSpec`](@ref)s) at the static field `static_field` along ẑ, for driving
at the angular frequency `Ω`, as a named tuple `(; shift, lower, upper)` with
`shift = upper - lower` (positive: the transition frequency increases).

`drives` is a complex Cartesian field phasor or a vector of
[`HarmonicDrive`](@ref)s over the coupled basis of *each* manifold (so the
phasor form is the natural one here). The two probed states may lie in the
same fine-structure manifold (e.g. a microwave qubit) or in two different
ones.
"""
function ac_zeeman_shift(
    species::HyperfineOneElectronSpecies,
    probe::Pair,
    static_field,
    Ω,
    drives;
    warn_ratio=20,
)
    lower, upper = probe
    fs_lower = fine_structure(parse_level(lower.level))
    fs_upper = fine_structure(parse_level(upper.level))
    m_lower = hyperfine_manifold(species, fs_lower, static_field)
    m_upper =
        fs_upper == fs_lower ? m_lower :
        hyperfine_manifold(species, fs_upper, static_field)
    δ_lower = ac_zeeman_shift(m_lower, lower, Ω, drives; warn_ratio)
    δ_upper = ac_zeeman_shift(m_upper, upper, Ω, drives; warn_ratio)
    (shift=δ_upper - δ_lower, lower=δ_lower, upper=δ_upper)
end

"""
    floquet_zeeman_shift(m::HyperfineManifold, state, Ω, drives; nharm = 8)

Returns the ac Zeeman shift of the adiabatically-labelled eigenstate `state` of
the hyperfine manifold solution `m` from the nonperturbative Floquet dressing
([`dress_manifold`](@ref), truncated at `nharm` harmonics): the difference
between the dressed quasienergy and the static eigen-energy.

Exact in the drive amplitude up to the harmonic truncation, so it serves as the
cross-check for [`ac_zeeman_shift`](@ref) near resonances; a warning is emitted
when the dressed-state identification becomes ambiguous (harmonic-0 overlap
below 0.5).
"""
function floquet_zeeman_shift(m::HyperfineManifold, state, Ω, drives; nharm=8)
    h0, h1 = eigenbasis_blocks(m, resolve_drives(m, drives))
    α = manifold_state_index(m, state)
    ε, _ = dress_manifold(h0, h1, Ω, nharm, 1; min_overlap=0.5)
    uconvert(u"µs^-1", ε[α] - m.energies[α])
end

"""
    ac_zeeman_shift(dt::DrivenTransition, i::Integer; warn_ratio = 20)

Returns the second-order perturbation-theory ac Zeeman shift of basis state `i`
of the driven-transition model from its harmonic drives, over the manifold
containing `i` (cf. the manifold blocks used by [`dress`](@ref)).
"""
function ac_zeeman_shift(dt::DrivenTransition, i::Integer; warn_ratio=20)
    range = i in dt.lower_range ? dt.lower_range : dt.upper_range
    h0, h1 = manifold_blocks(dt, range)
    ac_zeeman_shift(
        h0,
        h1,
        dt.drive_frequency,
        i - first(range) + 1;
        warn_ratio,
        warn_context=(; state=dt.basis[i]),
    )
end

export ac_zeeman_shift, ac_zeeman_contributions, floquet_zeeman_shift
