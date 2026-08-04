using Unitful
using WignerSymbols

"""
    photon_energy(laser)

Returns the photon energy ``ħ ω`` of a laser given either by its wavelength or
by its angular frequency.
"""
photon_energy(wavelength::Unitful.Length) =
    uconvert(u"J", 2π * u"ħ" * u"c" / wavelength)
photon_energy(frequency::Unitful.Frequency) = uconvert(u"J", u"ħ" * frequency)

"""
Returns ``⟨j m; 1 q | j' m + q⟩``, or zero where that is not a valid
electric-dipole channel.
"""
function dipole_cg(j, m, q, j_upper)
    abs(j - j_upper) <= 1 || return 0.0
    abs(m + q) <= j_upper || return 0.0
    Float64(clebschgordan(j, m, 1, q, j_upper, m + q))
end

"""
    polarisation_weights(ε)

Returns the fractions ``w_q = |ε_{-q}|^2`` of the beam intensity driving each of
the ``Δm = q`` electric-dipole channels, as a length-3 tuple indexed by
`q + 2`.

These are the squared magnitudes of the [`dipole_geometry`](@ref) amplitudes,
normalised to sum to one; only the direction of `ε` matters, as the overall
scale is set by the intensity.
"""
function polarisation_weights(ε)
    w = ntuple(i -> abs2(spherical_component(ε, 2 - i)), 3)
    total = sum(w)
    if total <= 0
        throw(ArgumentError("Polarisation vector must be non-zero"))
    end
    w ./ total
end

"""
    quadrupole_weights(ε, n)

Returns the weights ``|Γ_q|^2`` with which the beam intensity drives each of the
``Δm = q`` electric-quadrupole channels, as a length-5 tuple indexed by
`q + 3`.

These are the squared magnitudes of the [`quadrupole_geometry`](@ref) amplitudes
for a unit polarisation and beam direction, so only the directions of `ε` and
`n` matter; the overall scale is set by the intensity. Unlike the
electric-dipole [`polarisation_weights`](@ref) they do not sum to one, as only
part of ``ε ⊗ n`` is of rank two.

The polarisation must be transverse to the beam direction, ``ε ⋅ n = 0``, as it
is for any propagating field; pairs that do not belong to one physical beam are
rejected.
"""
function quadrupole_weights(ε, n)
    ε_scale = sqrt(sum(abs2, ε))
    n_scale = sqrt(sum(abs2, n))
    if iszero(ε_scale) || iszero(n_scale)
        throw(ArgumentError("Polarisation and beam direction must be non-zero"))
    end
    ε_unit = ε ./ ε_scale
    n_unit = n ./ n_scale
    if abs(sum(ε_unit .* n_unit)) > 1e-3
        throw(
            ArgumentError(
                "Polarisation must be transverse to the beam direction (ε ⋅ n = 0)",
            ),
        )
    end
    Γ = quadrupole_geometry(ε_unit, n_unit)
    ntuple(i -> abs2(Γ[i]), 5)
end

"""
Returns the polarisability of `state` for each of the three ``Δm`` channels
separately, as a length-3 vector indexed by `q + 2`.

For hyperfine states of a [`HyperfineOneElectronSpecies`](@ref), the explicit
channels are summed over the hyperfine levels of each intermediate
fine-structure level, with the ``F``-basis angular factors
(``\\mathrm{CG} × β``, cf. [`Levels.hyperfine_reduction`](@ref)) and the
intermediate-state hyperfine splittings resolved in the detunings (levels
without known hyperfine constants keep degenerate ``F`` levels at their
centroid). The hyperfine-resolved detunings are what gives e.g. the ⁴³Ca⁺
S``_{1/2}`` ``F`` levels their small tensor polarisability (``∝ A_P/Δ``),
which the fine-structure ``J = 1/2`` level lacks; with the splittings zeroed,
the result collapses exactly to the standard 6-j re-projection of the
fine-structure polarisability. The tensor part of the lumped static remainder
is re-projected with the expectation value of the ``J``-basis tensor operator
in the coupled state (diagonal-in-``F`` approximation; the ``F``-off-diagonal
elements are second order in the hyperfine mixing).
"""
function state_polarisabilities(species, state::StateSpec, ħω)
    level = convert(NoHyperfineNumberSpec, state.level)
    data = level_polarisability(species, level)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$(state.level)'"))
    end
    j, m = level.j, state.m
    if abs(m) > j || !isinteger(j - m)
        throw(ArgumentError("Invalid projection m = $m for level with j = $j"))
    end
    e_level = species.energies[level]

    α = fill(0.0u"C*m^2/V", 3)
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        r = d^2 / (2 * upper.j + 1)
        for q in -1:1
            rotating = dipole_cg(j, m, q, upper.j)^2
            counter = dipole_cg(j, m, -q, upper.j)^2
            α[q+2] +=
                uconvert(u"C*m^2/V", r * (rotating / (Δ - ħω) + counter / (Δ + ħω)))
        end
    end

    # The tensor part of the lumped remainder scales with the same angular factors
    # as an explicit channel would; the scalar part is channel-independent.
    tensor = j > 1//2 ? (3m^2 - j * (j + 1)) / (j * (2j - 1)) : 0//1
    for q in -1:1
        α[q+2] +=
            data.static_scalar + data.static_tensor * ((3 * (q == 0) - 1) / 2) * tensor
    end
    α
end

function state_polarisabilities(
    species::HyperfineOneElectronSpecies,
    state::StateSpec{HyperfineNumberSpec},
    ħω,
)
    spec = validate_hyperfine(species, state.level)
    fs = fine_structure(spec)
    data = level_polarisability(species, fs)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$(state.level)'"))
    end
    f, m = spec.f, state.m
    if abs(m) > f || !isinteger(f - m)
        throw(ArgumentError("Invalid projection m = $m for level with F = $f"))
    end
    e_level = level_energy(species, spec)

    # Zero-field hyperfine shift where constants are known; intermediate levels
    # without them keep their F levels degenerate at the centroid.
    shift_of(level) =
        haskey(species.hyperfine, fine_structure(level)) ?
        u"ħ" * hyperfine_shift(species, level) : zero(1.0u"J")

    α = fill(0.0u"C*m^2/V", 3)
    for (upper, d) in data.reduced_dipoles
        r = d^2 / (2 * upper.j + 1)
        for upper_f in hyperfine_levels(species, upper)
            β = hyperfine_reduction(species.nuclear_spin, spec, upper_f; rank=1)
            iszero(β) && continue
            Δ = species.energies[upper] + shift_of(upper_f) - e_level
            for q in -1:1
                rotating = dipole_cg(f, m, q, upper_f.f)^2
                counter = dipole_cg(f, m, -q, upper_f.f)^2
                α[q+2] += uconvert(
                    u"C*m^2/V",
                    r * β^2 * (rotating / (Δ - ħω) + counter / (Δ + ħω)),
                )
            end
        end
    end

    # The J-basis tensor operator of the lumped remainder, evaluated in the
    # coupled state (diagonal-in-F approximation): its angular factor is the
    # CG-weighted average of (3m_J² − J(J+1))/(J(2J−1)) over the state's
    # |m_I, m_J⟩ decomposition. Zero for J = 1/2, where the operator vanishes
    # identically — the genuinely-new hyperfine tensor term of such levels
    # comes from the resolved detunings above instead.
    j = fs.j
    i_nuc = species.nuclear_spin
    tensor = if j > 1//2
        sum(
            Float64(clebschgordan(i_nuc, m - m_j, j, m_j, f, m))^2 *
            (3m_j^2 - j * (j + 1)) / (j * (2j - 1)) for
            m_j in (-j):j if abs(m - m_j) <= i_nuc;
            init=0.0,
        )
    else
        0.0
    end
    for q in -1:1
        α[q+2] +=
            data.static_scalar + data.static_tensor * ((3 * (q == 0) - 1) / 2) * tensor
    end
    α
end

"""
Returns the near-resonant electric-quadrupole shift coefficients ``κ_q`` of the
given `lower` → `upper` transition, per ``Δm = q`` channel as a length-5 tuple
indexed by `q + 3`, or `nothing` if the two states are not connected by an
electric-quadrupole transition with a known Einstein A coefficient, or the
probed component itself is not a drivable one (``|Δm| ≤ 2``).

Contracted with the [`quadrupole_weights`](@ref) and scaled by the intensity
over the magnetic flux density, these give the shift; see
[`quadrupole_light_shift`](@ref) for the physics and the conventions.
"""
function quadrupole_shift_coefficients(species, lower::StateSpec, upper::StateSpec)
    lo = convert(NoHyperfineNumberSpec, lower.level)
    hi = convert(NoHyperfineNumberSpec, upper.level)
    multipole_rank(lo, hi) == 2 || return nothing
    abs(lo.j - hi.j) <= 2 <= lo.j + hi.j || return nothing
    Δm = upper.m - lower.m
    abs(Δm) <= 2 || return nothing
    a = einstein_a(species, lo, hi)
    isnothing(a) && return nothing
    ω = transition_frequency(species, lo, hi)

    # Ω² = rabi_scale × intensity × |⟨j m; 2 q|j' m'⟩ Γ_q|², cf. rabi_frequency().
    rabi_scale = 20π * u"c"^2 * a / (u"ħ" * ω^3)

    # With the laser on resonance with the probed transition, the line centre
    # drops out of every detuning: a channel `q` sharing the probed upper state
    # couples it to the lower sublevel at `upper.m - q`, and what is left of its
    # detuning is the Zeeman splitting of the two lower states involved — and
    # vice versa for the channels sharing the probed lower state. Both shift the
    # observed resonance in the same direction, hence the common sign below.
    w = fill(0.0u"µs*mT", 5)
    for q in -2:2
        q == Δm && continue
        m_lower = upper.m - q
        if abs(m_lower) <= lo.j
            c = clebschgordan(Float64, lo.j, m_lower, 2, q, hi.j, upper.m)
            other = StateSpec(lower.level, m_lower)
            w[q+3] += c^2 / zeeman_sensitivity(species, lower, other)
        end
        m_upper = lower.m + q
        if abs(m_upper) <= hi.j
            c = clebschgordan(Float64, lo.j, lower.m, 2, q, hi.j, m_upper)
            other = StateSpec(upper.level, m_upper)
            w[q+3] += c^2 / zeeman_sensitivity(species, other, upper)
        end
    end

    ntuple(i -> uconvert(u"m^2*T/J", -rabi_scale * w[i] / 4), 5)
end

"""
Raises an error unless `B` is a sensible signed scalar flux density for the
near-resonant quadrupole shift.
"""
function validate_shift_field(B)
    if !(B isa Unitful.BField)
        throw(
            ArgumentError(
                "B must be the signed scalar magnetic flux density along the " *
                "quantisation axis, not e.g. the Cartesian field vector " *
                "zeeman_hamiltonian() accepts",
            ),
        )
    end
    if iszero(B)
        throw(
            ArgumentError(
                "The near-resonant quadrupole shift is undefined at zero magnetic " *
                "field, where the Zeeman components are degenerate",
            ),
        )
    end
end

# Relative quadrupole amplitudes between the adiabatically-labelled eigenstates
# of two solved manifolds: the coupled-basis CG × β amplitudes conjugated with
# the eigenvector matrices ([upper eigenstate, lower eigenstate], real).
function rotated_quadrupole_amplitudes(
    species::HyperfineOneElectronSpecies,
    m_lower::HyperfineManifold,
    m_upper::HyperfineManifold,
)
    c = zeros(length(m_upper.basis), length(m_lower.basis))
    for (i, ls) in enumerate(m_lower.basis), (k, us) in enumerate(m_upper.basis)
        q = us.m - ls.m
        abs(q) <= 2 || continue
        β = hyperfine_reduction(species.nuclear_spin, ls.level, us.level; rank=2)
        iszero(β) && continue
        c[k, i] = Float64(clebschgordan(ls.level.f, ls.m, 2, q, us.level.f, us.m)) * β
    end
    m_upper.states' * c * m_lower.states
end

# The near-resonant shift coefficients of the probed eigen-pair (il, iu) from
# the rotated amplitudes and the exact eigen-energy detunings. As m_F is exact,
# every spectator pair belongs to a well-defined Δm channel; unlike the
# fine-structure case, spectators sharing the probed state's m_F (but a
# different F label) contribute too, at hyperfine-interval detunings.
function kappa_from_rotated(rotated, m_lower, m_upper, il, iu, rabi_scale)
    w = fill(0.0u"µs", 5)
    lower_m = m_lower.basis[il].m
    upper_m = m_upper.basis[iu].m
    for idx in eachindex(m_lower.energies)
        idx == il && continue
        q = upper_m - m_lower.basis[idx].m
        abs(q) <= 2 || continue
        w[Int(q)+3] +=
            abs2(rotated[iu, idx]) / (m_lower.energies[idx] - m_lower.energies[il])
    end
    for idx in eachindex(m_upper.energies)
        idx == iu && continue
        q = m_upper.basis[idx].m - lower_m
        abs(q) <= 2 || continue
        w[Int(q)+3] +=
            abs2(rotated[idx, il]) / (m_upper.energies[iu] - m_upper.energies[idx])
    end
    ntuple(i -> uconvert(u"µs^-1*m^2/W", -rabi_scale * w[i] / 4), 5)
end

"""
    quadrupole_shift_coefficients(species::HyperfineOneElectronSpecies,
                                  lower::StateSpec, upper::StateSpec, B)

Hyperfine form of the near-resonant shift coefficients: in the Breit–Rabi
regime the spectator detunings are the exact eigen-energy differences of the
manifolds at the static field `B` (neither linear in ``m`` nor sharing a
common ``1/B`` factor), so the field enters the *construction*, and the
coefficients directly scale the intensity (`shift = intensity × Σ w_q κ_q`),
with no ``1/B`` at evaluation time. The amplitudes are the eigenbasis-rotated
``F``-basis amplitudes, accounting exactly for the ``F`` mixing.
"""
function quadrupole_shift_coefficients(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec,
    upper::StateSpec,
    B,
)
    validate_shift_field(B)
    lo = parse_level(lower.level)
    hi = parse_level(upper.level)
    if !(lo isa HyperfineNumberSpec && hi isa HyperfineNumberSpec)
        throw(
            ArgumentError(
                "States must specify hyperfine (F) levels for a hyperfine species",
            ),
        )
    end
    fs_lo = fine_structure(lo)
    fs_hi = fine_structure(hi)
    multipole_rank(fs_lo, fs_hi) == 2 || return nothing
    abs(fs_lo.j - fs_hi.j) <= 2 <= fs_lo.j + fs_hi.j || return nothing
    abs(upper.m - lower.m) <= 2 || return nothing
    a = einstein_a(species, fs_lo, fs_hi)
    isnothing(a) && return nothing
    ω = transition_frequency(species, fs_lo, fs_hi)
    rabi_scale = 20π * u"c"^2 * a / (u"ħ" * ω^3)

    m_lower = hyperfine_manifold(species, fs_lo, B)
    m_upper = hyperfine_manifold(species, fs_hi, B)
    rotated = rotated_quadrupole_amplitudes(species, m_lower, m_upper)
    kappa_from_rotated(
        rotated,
        m_lower,
        m_upper,
        stateindex(m_lower.basis, StateSpec(lo, lower.m)),
        stateindex(m_upper.basis, StateSpec(hi, upper.m)),
        rabi_scale,
    )
end

"""
Entry of [`LightShiftCoefficients`](@ref)`.quadrupole_shifts`: the near-resonant
shift coefficients `κ` of one transition (cf.
[`quadrupole_shift_coefficients`](@ref)) and the photon energy `ħω` of its
resonance.
"""
const QuadrupoleShiftEntry =
    @NamedTuple{ħω::typeof(1.0u"J"), κ::NTuple{5,typeof(1.0u"m^2*T/J")}}

"""
Hyperfine counterpart of [`Levels.QuadrupoleShiftEntry`](@ref): the
coefficients are field-resolved (computed at the stored flux density `B`, which
evaluation checks) and scale the intensity directly, with no ``1/B``.
"""
const HyperfineQuadrupoleShiftEntry = @NamedTuple{
    ħω::typeof(1.0u"J"),
    B::typeof(1.0u"mT"),
    κ::NTuple{5,typeof(1.0u"µs^-1*m^2/W")},
}

"""
Light-shift data for a set of states, precomputed for one laser frequency.

Construct with [`LightShiftCoefficients`](@ref)`(species, basis, laser)` and
evaluate with [`light_shift`](@ref); see there for the sign and unit
conventions.
"""
struct LightShiftCoefficients{L<:LevelSpec,E<:Quantity,T<:Quantity,Q}
    "The states the coefficients refer to, fixing the row order."
    basis::StateBasis{L}

    "The photon energy the coefficients were computed for."
    photon_energy::E

    """
    Polarisability of each state for each ``Δm`` channel, as a
    `length(basis) × 3` matrix indexed by `[stateindex, q + 2]`.

    States of levels without [`LevelPolarisability`](@ref) data carry `NaN`
    entries; evaluating an electric-dipole shift involving them raises an
    error.
    """
    polarisabilities::Matrix{T}

    """
    Near-resonant electric-quadrupole shift coefficients of the transitions
    between the basis states, keyed by their `(lower, upper)` pair of basis
    indices (cf. [`quadrupole_light_shift`](@ref)).

    Only pairs that are actually connected by an electric-quadrupole transition
    with known data appear. Unlike the polarisabilities, the coefficients do
    not depend on the laser frequency, which the model instead fixes to
    resonance with the probed transition; the photon energy of that resonance
    is stored alongside them so that evaluation can check the premise against
    the frequency the coefficients were computed for. The entries are
    [`Levels.QuadrupoleShiftEntry`](@ref)s for a no-hyperfine basis, or
    field-resolved [`Levels.HyperfineQuadrupoleShiftEntry`](@ref)s for a
    hyperfine one.
    """
    quadrupole_shifts::Dict{Tuple{Int,Int},Q}
end

"""
    LightShiftCoefficients(species, basis::StateBasis, laser)
    LightShiftCoefficients(species, levels_or_states::AbstractVector, laser)

Precomputes the ac Stark shift of every state in the given basis for a laser of
the given wavelength or angular frequency.

All the atomic structure enters here, so that evaluating the shift for a
particular intensity and polarisation afterwards costs only a handful of
arithmetic operations — the intended use when fitting laser parameters against
many observed transition frequencies.

The electric-dipole part requires [`LevelPolarisability`](@ref) data. Levels
without it are still admitted, but evaluating any shift involving their
polarisability raises an error. The near-resonant electric-quadrupole
coefficients only need the Einstein A coefficient; they are computed for
whichever pairs of basis states support them, and simply left out for the rest.

For a hyperfine basis (of a [`HyperfineOneElectronSpecies`](@ref)), the
quadrupole coefficients are field-resolved (cf. the four-argument
[`Levels.quadrupole_shift_coefficients`](@ref)), so the static flux density
must be supplied via the `B` keyword if they are wanted; without it only the
electric-dipole part is precomputed. Evaluation then checks that the `B` it is
given matches the one the coefficients were built for.
"""
function LightShiftCoefficients(
    species,
    basis::StateBasis{NoHyperfineNumberSpec},
    laser,
)
    ħω = photon_energy(laser)
    rows = [
        isnothing(level_polarisability(species, state.level)) ?
        fill(NaN * u"C*m^2/V", 3) : state_polarisabilities(species, state, ħω) for
        state in basis
    ]
    quadrupole = Dict{Tuple{Int,Int},QuadrupoleShiftEntry}()
    for (i, lower) in enumerate(basis), (k, upper) in enumerate(basis)
        κ = quadrupole_shift_coefficients(species, lower, upper)
        isnothing(κ) && continue
        resonance = uconvert(
            u"J",
            u"ħ" * transition_frequency(species, lower.level, upper.level),
        )
        quadrupole[(i, k)] = (ħω=resonance, κ=κ)
    end
    LightShiftCoefficients(basis, ħω, permutedims(reduce(hcat, rows)), quadrupole)
end

function LightShiftCoefficients(
    species::HyperfineOneElectronSpecies,
    basis::StateBasis{HyperfineNumberSpec},
    laser;
    B=nothing,
)
    ħω = photon_energy(laser)
    rows = [
        isnothing(level_polarisability(species, state.level)) ?
        fill(NaN * u"C*m^2/V", 3) : state_polarisabilities(species, state, ħω) for
        state in basis
    ]
    quadrupole = Dict{Tuple{Int,Int},HyperfineQuadrupoleShiftEntry}()
    if !isnothing(B)
        validate_shift_field(B)
        # Hoist the manifold solutions and rotated amplitude matrices out of
        # the pair loop; per pair only the O(n) spectator sums remain.
        manifolds = Dict(
            fs => hyperfine_manifold(species, fs, B) for
            fs in unique!(fine_structure.(basis.levels))
        )
        for (fs_lo, m_lower) in manifolds, (fs_hi, m_upper) in manifolds
            fs_lo == fs_hi && continue
            multipole_rank(fs_lo, fs_hi) == 2 || continue
            abs(fs_lo.j - fs_hi.j) <= 2 <= fs_lo.j + fs_hi.j || continue
            a = einstein_a(species, fs_lo, fs_hi)
            isnothing(a) && continue
            ω = transition_frequency(species, fs_lo, fs_hi)
            rabi_scale = 20π * u"c"^2 * a / (u"ħ" * ω^3)
            rotated = rotated_quadrupole_amplitudes(species, m_lower, m_upper)
            for (i, lower) in enumerate(basis), (k, upper) in enumerate(basis)
                fine_structure(lower.level) == fs_lo || continue
                fine_structure(upper.level) == fs_hi || continue
                abs(upper.m - lower.m) <= 2 || continue
                κ = kappa_from_rotated(
                    rotated,
                    m_lower,
                    m_upper,
                    stateindex(m_lower.basis, lower),
                    stateindex(m_upper.basis, upper),
                    rabi_scale,
                )
                resonance = uconvert(
                    u"J",
                    u"ħ" * transition_frequency(species, lower.level, upper.level),
                )
                quadrupole[(i, k)] = (ħω=resonance, B=uconvert(u"mT", B), κ=κ)
            end
        end
    end
    LightShiftCoefficients(basis, ħω, permutedims(reduce(hcat, rows)), quadrupole)
end

LightShiftCoefficients(species, levels_or_states::AbstractVector, laser; kwargs...) =
    LightShiftCoefficients(species, StateBasis(levels_or_states), laser; kwargs...)

"""
Returns the polarisability of the basis state with the given index for the
[`polarisation_weights`](@ref) `w`.
"""
function state_polarisability(c::LightShiftCoefficients, index::Integer, w)
    if isnan(c.polarisabilities[index, 1])
        throw(
            ArgumentError(
                "No light-shift data known for level '$(c.basis[index].level)'",
            ),
        )
    end
    w[1] * c.polarisabilities[index, 1] +
    w[2] * c.polarisabilities[index, 2] +
    w[3] * c.polarisabilities[index, 3]
end

"""
Converts a polarisability into the corresponding shift at the given intensity.
"""
shift_at_intensity(α, intensity) =
    uconvert(u"µs^-1", -intensity * α / (2 * u"c" * u"ε0" * u"ħ"))

"""
Contracts near-resonant quadrupole shift coefficients with the
[`quadrupole_weights`](@ref) `w` for the given intensity and magnetic flux
density.
"""
function quadrupole_shift_at(κ, intensity, w, B)
    validate_shift_field(B)
    uconvert(u"µs^-1", (intensity / B) * sum(w .* κ))
end

"""
Evaluates a stored quadrupole-shift entry: the fine-structure form scales the
coefficients by intensity over `B`; the field-resolved hyperfine form scales by
the intensity alone, after checking `B` against the field the coefficients were
computed for.
"""
entry_shift(entry::QuadrupoleShiftEntry, intensity, w, B) =
    quadrupole_shift_at(entry.κ, intensity, w, B)

function entry_shift(entry::HyperfineQuadrupoleShiftEntry, intensity, w, B)
    validate_shift_field(B)
    if !isapprox(B, entry.B; rtol=1e-6)
        throw(
            ArgumentError(
                "The near-resonant quadrupole shift coefficients were computed " *
                "for B = $(entry.B), not $B; rebuild the LightShiftCoefficients " *
                "for the new field",
            ),
        )
    end
    uconvert(u"µs^-1", intensity * sum(w .* entry.κ))
end

"""
Returns the electric-dipole shift of the transition between the basis states
with the given indices, for the [`polarisation_weights`](@ref) `w`.
"""
dipole_transition_shift(
    c::LightShiftCoefficients,
    lower::Integer,
    upper::Integer,
    intensity,
    w,
) = shift_at_intensity(
    state_polarisability(c, upper, w) - state_polarisability(c, lower, w),
    intensity,
)

"""
Returns the near-resonant electric-quadrupole shift of the transition between
the basis states with the given indices, for the [`quadrupole_weights`](@ref)
`w`.
"""
function quadrupole_transition_shift(
    c::LightShiftCoefficients,
    lower::Integer,
    upper::Integer,
    intensity,
    w,
    B,
)
    entry = get(c.quadrupole_shifts, (lower, upper), nothing)
    if isnothing(entry)
        levels(i, k) = (c.basis[i].level, c.basis[k].level)
        passed = levels(lower, upper)
        if any(levels(i, k) == passed for (i, k) in keys(c.quadrupole_shifts))
            # Other components of the same two levels do have coefficients, so
            # only the probed component itself can be at fault.
            Δm = c.basis[upper].m - c.basis[lower].m
            Δm = isinteger(Δm) ? Int(Δm) : Δm
            throw(
                ArgumentError(
                    "The Δm = $Δm component '$(c.basis[lower])' => " *
                    "'$(c.basis[upper])' of an electric-quadrupole transition " *
                    "cannot be driven (|Δm| ≤ 2), so there is no resonance whose " *
                    "shift could be observed",
                ),
            )
        elseif any(levels(k, i) == passed for (i, k) in keys(c.quadrupole_shifts))
            throw(
                ArgumentError(
                    "'$(c.basis[lower])' is of higher energy than " *
                    "'$(c.basis[upper])'; give the transition as lower => upper",
                ),
            )
        end
        throw(
            ArgumentError(
                "'$(c.basis[lower])' and '$(c.basis[upper])' are not connected by " *
                "an electric-quadrupole transition with a known Einstein A " *
                "coefficient",
            ),
        )
    end
    if !isapprox(c.photon_energy, entry.ħω; rtol=1e-3)
        wavelength(ħω) = round(u"nm", 2π * u"ħ" * u"c" / ħω; digits=3)
        throw(
            ArgumentError(
                "The near-resonant quadrupole shift presumes the laser to be tuned " *
                "to resonance with the probed transition, but the coefficients were " *
                "computed for $(wavelength(c.photon_energy)) rather than the " *
                "$(wavelength(entry.ħω)) of '$(c.basis[lower])' => '$(c.basis[upper])'",
            ),
        )
    end
    entry_shift(entry, intensity, w, B)
end

"""
Raises an error unless the beam direction `n` and the magnetic flux density `B`
keyword arguments are given either both or not at all.
"""
function check_geometry_keywords(n, B)
    if isnothing(n) != isnothing(B)
        throw(
            ArgumentError(
                "The beam direction n and the magnetic flux density B must be " *
                "given together",
            ),
        )
    end
end

"""
    light_shift(coefficients::LightShiftCoefficients, state::StateSpec, intensity, ε)
    light_shift(coefficients::LightShiftCoefficients, transition::Pair, intensity, ε[; n, B])

Returns the ac Stark (light) shift of a state, or of a `lower => upper`
transition, from the precomputed `coefficients` (in angular units).

For a transition this is the shift of the upper state minus that of the lower
one, i.e. the amount by which the frequency at which resonant Rabi flopping is
observed exceeds the unperturbed transition frequency.

Without the keyword arguments, only the electric-dipole contributions are
included. Those are what shifts an isolated state, and they do not depend on
the beam direction — an E1 light shift is set by the polarisation `ε` alone,
whose overall scale is irrelevant as the intensity is given separately.

Passing the beam direction `n` and the magnetic flux density `B` (always
together) additionally accounts for the near-resonant coupling to the other
Zeeman components of the probed quadrupole transition itself, which is only
defined for a transition and does depend on the beam direction; see
[`quadrupole_light_shift`](@ref) for the model and its range of validity.
"""
light_shift(c::LightShiftCoefficients, state::StateSpec, intensity, ε) =
    shift_at_intensity(
        state_polarisability(c, stateindex(c.basis, state), polarisation_weights(ε)),
        intensity,
    )

function light_shift(
    c::LightShiftCoefficients,
    transition::Pair,
    intensity,
    ε;
    n=nothing,
    B=nothing,
)
    check_geometry_keywords(n, B)
    lower = stateindex(c.basis, transition.first)
    upper = stateindex(c.basis, transition.second)
    shift = dipole_transition_shift(c, lower, upper, intensity, polarisation_weights(ε))
    isnothing(n) && return shift
    shift +
    quadrupole_transition_shift(c, lower, upper, intensity, quadrupole_weights(ε, n), B)
end

"""
    light_shift(species, state::StateSpec, laser, intensity, ε)
    light_shift(species, lower::StateSpec, upper::StateSpec, laser, intensity, ε[; n, B])

Returns the ac Stark (light) shift of a single state or transition for a laser
of the given wavelength or angular frequency, intensity and polarisation.

This is a convenience wrapper that does the full sum over intermediate levels on
every call; use [`LightShiftCoefficients`](@ref) to hoist that work out of a
loop over intensities or polarisations.

As for the precomputed form, giving the beam direction `n` and the magnetic flux
density `B` as keywords adds the near-resonant
[`quadrupole_light_shift`](@ref).
"""
light_shift(species, state::StateSpec, laser, intensity, ε) = shift_at_intensity(
    sum(
        polarisation_weights(ε) .*
        state_polarisabilities(species, state, photon_energy(laser)),
    ),
    intensity,
)

function light_shift(
    species,
    lower::StateSpec,
    upper::StateSpec,
    laser,
    intensity,
    ε;
    n=nothing,
    B=nothing,
)
    check_geometry_keywords(n, B)
    ħω = photon_energy(laser)
    w = polarisation_weights(ε)
    α_lower = sum(w .* state_polarisabilities(species, lower, ħω))
    α_upper = sum(w .* state_polarisabilities(species, upper, ħω))
    shift = shift_at_intensity(α_upper - α_lower, intensity)
    isnothing(n) && return shift
    shift + quadrupole_light_shift(species, lower, upper, intensity, ε; n, B)
end

"""
    quadrupole_light_shift(coefficients::LightShiftCoefficients, transition::Pair, intensity, ε; n, B)
    quadrupole_light_shift(species, lower::StateSpec, upper::StateSpec, intensity, ε; n, B)

Returns the near-resonant electric-quadrupole contribution to the light shift of
a `lower => upper` transition (in angular units), assuming the laser is tuned to
resonance with that transition. When evaluating from precomputed
[`LightShiftCoefficients`](@ref), an error is raised if the laser frequency they
were computed for is inconsistent with that premise (i.e. belongs to a different
transition).

Driving one Zeeman component of a quadrupole transition also couples the two
states involved, off resonantly, to every other component sharing one of them;
the resulting ac Stark shifts move the observed resonance. This is the shift for
which ⁸⁸Sr⁺ clock evaluations quote an "E2 ac Stark shift" (`[Lindvall2025]`,
Sec. III F 2): with a perfectly linear polarisation the shifts of the two
components of a ``±m`` Zeeman pair are equal and opposite and hence average out,
whereas an elliptical polarisation leaves a net shift.

The detunings involved are the Zeeman splittings of the two manifolds, so the
result scales as the intensity over the magnetic flux density `B` — the signed
component along the quantisation axis. Both the polarisation `ε` and the beam
direction `n` enter, through the [`quadrupole_weights`](@ref); the laser
frequency does not, as it is pinned to the probed transition.

This is second-order perturbation theory in the coupling, so it holds as long as
every Rabi frequency involved is small compared with the Zeeman splittings, and
breaks down for near-degenerate components at very low fields. Counter-rotating
terms and the far-detuned quadrupole channels to other levels are neglected;
both are smaller by many orders of magnitude, and the latter are in any case
dwarfed by the dipole contributions to [`light_shift`](@ref).

# References

- `[Lindvall2025]`: T. Lindvall, A. E. Wallin, K. J. Hanhijärvi, and T. Fordell,
  "⁸⁸Sr⁺ Optical Clock with 7.9 × 10⁻¹⁹ Systematic Uncertainty and Measurement of
  Its Absolute Frequency", Phys. Rev. Applied **24**, 044082 (2025),
  [doi:10.1103/cztf-bfvp](https://doi.org/10.1103/cztf-bfvp).
"""
function quadrupole_light_shift(
    c::LightShiftCoefficients,
    transition::Pair,
    intensity,
    ε;
    n,
    B,
)
    quadrupole_transition_shift(
        c,
        stateindex(c.basis, transition.first),
        stateindex(c.basis, transition.second),
        intensity,
        quadrupole_weights(ε, n),
        B,
    )
end

# For a hyperfine species, the coefficients are only defined at a given static
# field; direct the caller to the four-argument form.
function quadrupole_shift_coefficients(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec,
    upper::StateSpec,
)
    throw(
        ArgumentError(
            "The near-resonant quadrupole shift of a hyperfine species requires " *
            "the static field: pass B (cf. quadrupole_shift_coefficients(species, " *
            "lower, upper, B))",
        ),
    )
end

function quadrupole_light_shift(
    species::HyperfineOneElectronSpecies,
    lower::StateSpec,
    upper::StateSpec,
    intensity,
    ε;
    n,
    B,
)
    κ = quadrupole_shift_coefficients(species, lower, upper, B)
    if isnothing(κ)
        lo = fine_structure(parse_level(lower.level))
        hi = fine_structure(parse_level(upper.level))
        e2 = multipole_rank(lo, hi) == 2 && abs(lo.j - hi.j) <= 2 <= lo.j + hi.j
        if e2 && !isnothing(einstein_a(species, lo, hi))
            Δm = upper.m - lower.m
            Δm = isinteger(Δm) ? Int(Δm) : Δm
            throw(
                ArgumentError(
                    "The Δm = $Δm component '$lower' => '$upper' of an " *
                    "electric-quadrupole transition cannot be driven (|Δm| ≤ 2), " *
                    "so there is no resonance whose shift could be observed",
                ),
            )
        elseif e2 && !isnothing(einstein_a(species, hi, lo))
            throw(
                ArgumentError(
                    "'$lower' is of higher energy than '$upper'; give the " *
                    "transition as lower => upper",
                ),
            )
        end
        throw(
            ArgumentError(
                "'$lower' and '$upper' are not connected by an electric-quadrupole " *
                "transition with a known Einstein A coefficient",
            ),
        )
    end
    uconvert(u"µs^-1", intensity * sum(quadrupole_weights(ε, n) .* κ))
end

function quadrupole_light_shift(
    species,
    lower::StateSpec,
    upper::StateSpec,
    intensity,
    ε;
    n,
    B,
)
    κ = quadrupole_shift_coefficients(species, lower, upper)
    if isnothing(κ)
        lo = convert(NoHyperfineNumberSpec, lower.level)
        hi = convert(NoHyperfineNumberSpec, upper.level)
        e2 = multipole_rank(lo, hi) == 2 && abs(lo.j - hi.j) <= 2 <= lo.j + hi.j
        if e2 && !isnothing(einstein_a(species, lo, hi))
            # The levels are connected in the given order, so only the probed
            # component itself can be at fault.
            Δm = upper.m - lower.m
            Δm = isinteger(Δm) ? Int(Δm) : Δm
            throw(
                ArgumentError(
                    "The Δm = $Δm component '$lower' => '$upper' of an " *
                    "electric-quadrupole transition cannot be driven (|Δm| ≤ 2), " *
                    "so there is no resonance whose shift could be observed",
                ),
            )
        elseif e2 && !isnothing(einstein_a(species, hi, lo))
            throw(
                ArgumentError(
                    "'$lower' is of higher energy than '$upper'; give the " *
                    "transition as lower => upper",
                ),
            )
        end
        throw(
            ArgumentError(
                "'$lower' and '$upper' are not connected by an electric-quadrupole " *
                "transition with a known Einstein A coefficient",
            ),
        )
    end
    quadrupole_shift_at(κ, intensity, quadrupole_weights(ε, n), B)
end

"""
Returns the level and its light-shift data, or raises an error if there is none.
"""
function polarisability_data(species, level)
    spec = convert(NoHyperfineNumberSpec, level)
    data = level_polarisability(species, spec)
    if isnothing(data)
        throw(ArgumentError("No light-shift data known for level '$level'"))
    end
    spec, data
end

"""
Returns the per-channel polarisabilities of the stretched state ``m = F`` of a
hyperfine level, from which the ``α_0``/``α_1``/``α_2`` decomposition is
extracted (``α_{σ^-} = α_0 - α_1/2 - α_2/2``, ``α_π = α_0 + α_2``,
``α_{σ^+} = α_0 + α_1/2 - α_2/2`` there, as the tensor angular factor is unity
at ``m = F``).
"""
stretched_channels(
    species::HyperfineOneElectronSpecies,
    spec::HyperfineNumberSpec,
    ħω,
) = state_polarisabilities(species, StateSpec(spec, spec.f), ħω)

function stretched_channels(species, spec, ħω)
    throw(
        ArgumentError(
            "Level '$spec' carries hyperfine structure, but the species has none",
        ),
    )
end

"""
    scalar_polarisability(species, level, laser)

Returns the dynamic scalar polarisability ``α_0(ω)`` of the given level at the
given laser wavelength or angular frequency.

Together with [`vector_polarisability`](@ref) and
[`tensor_polarisability`](@ref), this decomposes the polarisability of the
sublevel ``m`` for polarisation `ε` as

``α = α_0 + 𝒜 \\frac{m}{2j} α_1 +
\\frac{3 |ε_0|^2 - 1}{2} \\frac{3 m^2 - j(j+1)}{j(2j-1)} α_2``,

where ``𝒜 = |ε_{-1}|^2 - |ε_{+1}|^2`` is the degree of circular polarisation
about the quantisation axis and ``ε_0 = ε · ẑ``. [`light_shift`](@ref) evaluates
the equivalent sum directly, without going through this decomposition.

For a hyperfine ``F`` level (with ``j`` replaced by ``F`` in the decomposition),
the components are extracted from the per-channel polarisabilities of the
stretched state; cf. [`Levels.state_polarisabilities`](@ref) for the
approximations involved.
"""
function scalar_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        α = stretched_channels(species, parsed, ħω)
        return (α[1] + α[2] + α[3]) / 3
    end
    spec, data = polarisability_data(species, level)
    e_level = species.energies[spec]
    α = data.static_scalar
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        α += uconvert(u"C*m^2/V", d^2 / (3 * (2 * spec.j + 1)) * 2Δ / (Δ^2 - ħω^2))
    end
    α
end

"""
    vector_polarisability(species, level, laser)

Returns the dynamic vector polarisability ``α_1(ω)`` of the given level at the
given laser wavelength or angular frequency (see
[`scalar_polarisability`](@ref) for the convention).

The vector polarisability vanishes in the static limit, so the far-detuned
remainder lumped into [`LevelPolarisability`](@ref)`.static_scalar` contributes
nothing to it.
"""
function vector_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f > 0 || return 0.0u"C*m^2/V"
        α = stretched_channels(species, parsed, ħω)
        return α[3] - α[1]
    end
    spec, data = polarisability_data(species, level)
    e_level = species.energies[spec]
    j = spec.j
    α = 0.0u"C*m^2/V"
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        weight = dipole_cg(j, j, 1, upper.j)^2 - dipole_cg(j, j, -1, upper.j)^2
        α += uconvert(u"C*m^2/V", d^2 / (2 * upper.j + 1) * weight * 2ħω / (Δ^2 - ħω^2))
    end
    α
end

"""
    tensor_polarisability(species, level, laser)

Returns the dynamic tensor polarisability ``α_2(ω)`` of the given level at the
given laser wavelength or angular frequency (see
[`scalar_polarisability`](@ref) for the convention).

Zero for ``j ≤ 1/2``, which has no oriented sublevels to distinguish.
"""
function tensor_polarisability(species, level, laser)
    ħω = photon_energy(laser)
    parsed = parse_level(level)
    if parsed isa HyperfineNumberSpec
        parsed.f >= 1 || return 0.0u"C*m^2/V"
        α = stretched_channels(species, parsed, ħω)
        return (2α[2] - α[1] - α[3]) / 3
    end
    spec, data = polarisability_data(species, level)
    j = spec.j
    j > 1//2 || return 0.0u"C*m^2/V"
    e_level = species.energies[spec]
    α = data.static_tensor
    for (upper, d) in data.reduced_dipoles
        Δ = species.energies[upper] - e_level
        # The π-polarised shift of the stretched state m = j is exactly
        # α_0 + α_2, as the tensor angular factor is unity there.
        weight = dipole_cg(j, j, 0, upper.j)^2 / (2 * upper.j + 1) - 1 / (3 * (2j + 1))
        α += uconvert(u"C*m^2/V", d^2 * weight * 2Δ / (Δ^2 - ħω^2))
    end
    α
end

export LightShiftCoefficients,
    light_shift,
    quadrupole_light_shift,
    scalar_polarisability,
    tensor_polarisability,
    vector_polarisability
public photon_energy, polarisation_weights, quadrupole_weights
