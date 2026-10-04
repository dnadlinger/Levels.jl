# The Lindblad model of the internal states: rotating-frame Hamiltonian (with
# any harmonic beat terms), spontaneous-emission and laser-dephasing jump
# operators, all assembled from the atomic data.

"""
Label of a spontaneous-emission jump operator: the fine-structure level pair,
the spherical component ``q = m_\\mathrm{upper} - m_\\mathrm{lower}`` of the emitted
photon and, for component-resolved operators, the single `lower => upper` state
pair it connects (`nothing` for the polarisation-grouped form).
"""
struct DecayLabel
    lower::NoHyperfineNumberSpec
    upper::NoHyperfineNumberSpec
    q::Int
    component::Union{Nothing,Pair{StateSpec,StateSpec}}
end

"""
Label of a laser-phase-diffusion jump operator: the index of the beam.
"""
struct DephasingLabel
    beam::Int
end

"""
The Lindblad master equation of the internal states of a [`LaserScheme`](@ref)
in the rotating frame of its [`RotatingFrame`](@ref),

```math
\\dot ρ = -i [H(t), ρ] + \\sum_k \\left( L_k ρ L_k^† - \\tfrac{1}{2} \\{ L_k^† L_k, ρ \\} \\right),
\\qquad
H(t) = H_0 + \\sum_j \\left( M_j e^{-i w_j t} + M_j^† e^{+i w_j t} \\right),
```

over the states of `basis`: `hamiltonian` is the static ``H_0`` (Hermitian,
µs⁻¹) — the frame diagonal plus the couplings ``(C + C^†)/2`` of all
time-independent beam components, with `C` from [`coupling_matrix`](@ref) — and
`harmonics` the ``(w_j, M_j)`` pairs of the components left with a beat note
(empty for a static model). `jump_operators` (µs⁻¹ᐟ²) are labelled by
`jump_labels` ([`DecayLabel`](@ref), [`DephasingLabel`](@ref)).

For a hyperfine species the basis states denote the adiabatically-labelled
eigenstates at the static field, with exact at-field amplitudes throughout.
Constructed via [`lindblad_model`](@ref).
"""
struct LindbladModel{B<:StateBasis,F<:RotatingFrame,H<:Quantity,W<:Quantity,J<:Quantity}
    basis::B
    frame::F
    hamiltonian::Matrix{H}
    harmonics::Vector{Tuple{W,Matrix{H}}}
    jump_operators::Vector{Matrix{J}}
    jump_labels::Vector{Union{DecayLabel,DephasingLabel}}
end

# Energies of the basis states at the static field relative to the zero-field
# centroids of their fine-structure levels — the one place the species kind is
# dispatched on: first-order Zeeman shifts for a fine-structure species, exact
# hyperfine + Zeeman eigen-energies (adiabatic labels) for a hyperfine one.
function centroid_energies(species, basis::StateBasis{NoHyperfineNumberSpec}, B)
    [uconvert(ANGULAR_UNIT, zeeman_shift(species, state, B)) for state in basis]
end

function centroid_energies(species, basis::StateBasis{HyperfineNumberSpec}, B)
    manifolds = Dict(
        fs => hyperfine_manifold(species, fs, B) for fs in fine_structure_levels(basis)
    )
    [
        uconvert(
            ANGULAR_UNIT,
            state_energy(manifolds[fine_structure(state.level)], state),
        ) for state in basis
    ]
end

# All (lower, upper) fine-structure pairs of the basis levels connected by a
# known decay, plus — per upper level — the decay channels leaving the basis.
function decay_pairs(species, levels)
    pairs = Tuple{NoHyperfineNumberSpec,NoHyperfineNumberSpec}[]
    open = Dict{NoHyperfineNumberSpec,Vector{NoHyperfineNumberSpec}}()
    for hi in levels
        for ((lo, up), _) in species.einstein_as
            up == hi || continue
            if lo in levels
                push!(pairs, (lo, hi))
            else
                push!(get!(open, hi, NoHyperfineNumberSpec[]), lo)
            end
        end
    end
    sort!(pairs; by=p -> (string(p[2]), string(p[1])))
    pairs, open
end

"""
    lindblad_model(species, scheme::LaserScheme; decay = :grouped, open_channels = :error)

Assembles the [`LindbladModel`](@ref) of the scheme from the atomic data.

**Hamiltonian.** The frame diagonal holds each state's energy at the static
field relative to its level's zero-field centroid (Zeeman shift, or the exact
hyperfine + Zeeman eigen-energy) minus its frame frequency of
[`rotating_frame`](@ref); every beam adds ``(C + C^†)/2`` with the complex Rabi
frequencies of [`coupling_matrix`](@ref) over the two fine-structure manifolds it
connects (all their ``F`` levels for a hyperfine species) — into the static part
for its time-independent components, as harmonic terms for those with a beat
note.

**Spontaneous emission.** For every pair of basis levels with an Einstein A
coefficient (rank ``k`` from [`multipole_rank`](@ref)), one jump operator per
emitted spherical component ``q = -k…k``,
``L_q = \\sqrt{A} \\sum_{m' - m = q} ⟨\\mathrm{up}|T^k_q|\\mathrm{lo}⟩_\\mathrm{rel} \\,
|\\mathrm{lo}⟩⟨\\mathrm{up}|`` with the relative amplitudes of
[`Levels.transition_amplitude`](@ref) (exact at the static field) — the
grouping appropriate when the Zeeman (hyperfine) splittings are unresolved
against the photon bandwidth, ``\\sum_q L_q^† L_q = A \\, P_\\mathrm{upper}``.
`decay = :resolved` instead gives one operator per state pair (appropriate for
well-resolved hyperfine structure; in a per-state frame — `frame_kind =
:states` — it is also what makes an operator whose components see different
frame differences time-independent, at the price of the coherence-transfer
terms between them, and is then required). An upper level with decay channels
to levels missing from the basis raises an error naming them, unless
`open_channels = :drop` discards those channels (the level then lives longer
than its physical lifetime); a basis holding only part of a level's states
likewise silently drops the channels into the missing states.

**Laser linewidth.** A beam with non-zero `linewidth` ``γ`` (FWHM) adds the
phase-diffusion operator ``\\sqrt{γ} \\sum_i n_{b,i} P_i`` with ``P_i`` the state
projectors and ``n_{b,i}`` the signed traversal counts of the beam on the frame
tree, so that a coherence between two states decays at half the summed
linewidths of the beams whose frequencies define their relative frame — a
result independent of the frame reference. Linewidths of beams outside the
frame tree are not representable this way and raise an error.
"""
function lindblad_model(
    species,
    scheme::LaserScheme;
    decay::Symbol=:grouped,
    open_channels::Symbol=:error,
)
    decay in (:grouped, :resolved) ||
        throw(ArgumentError("decay must be :grouped or :resolved, got $decay"))
    open_channels in (:error, :drop) || throw(
        ArgumentError("open_channels must be :error or :drop, got $open_channels"),
    )

    basis = scheme.basis
    B = scheme.static_field
    num_states = length(basis)
    couplings = beam_couplings(species, scheme)
    frame = rotating_frame(species, scheme, couplings)
    levels = frame.levels

    # --- Hamiltonian ---------------------------------------------------------
    energies = centroid_energies(species, basis, B)
    H = zeros(typeof(complex(1.0) * ANGULAR_UNIT), num_states, num_states)
    for i in 1:num_states
        H[i, i] = complex(1.0) * (energies[i] - frame.offsets[i])
    end
    harmonics = Tuple{typeof(1.0 * ANGULAR_UNIT),typeof(H)}[]
    for (b, C) in enumerate(couplings)
        static, parts = split_beam_coupling(frame, b, C)
        H .+= (static .+ static') ./ 2
        for (w, part) in parts
            w > zero(w) ? push!(harmonics, (w, part ./ 2)) :
            push!(harmonics, (-w, Matrix(part') ./ 2))
        end
    end
    harmonics = merge_harmonics(harmonics)

    # --- Spontaneous emission ------------------------------------------------
    J = typeof(complex(1.0) * JUMP_UNIT)
    jump_operators = Matrix{J}[]
    jump_labels = Union{DecayLabel,DephasingLabel}[]
    pairs, open = decay_pairs(species, levels)
    if !isempty(open) && open_channels == :error
        desc = join(
            ("'$hi' → " * join(("'$lo'" for lo in los), ", ") for (hi, los) in open),
            "; ",
        )
        throw(
            ArgumentError(
                "Decay channels leave the basis: $desc. Add the levels, or pass " *
                "open_channels = :drop to discard them.",
            ),
        )
    end
    for (lo, hi) in pairs
        rank = multipole_rank(lo, hi)
        a = einstein_a(species, lo, hi)
        lower_levels = unique!([l for l in basis.levels if fine_structure(l) == lo])
        upper_levels = unique!([l for l in basis.levels if fine_structure(l) == hi])
        amplitude = amplitude_evaluator(species, lower_levels, upper_levels, B)
        for q in (-rank):rank
            L = zeros(J, num_states, num_states)
            components = Tuple{Int,Int}[]
            for (i, s_lo) in enumerate(basis)
                fine_structure(s_lo.level) == lo || continue
                for (k, s_hi) in enumerate(basis)
                    fine_structure(s_hi.level) == hi || continue
                    s_hi.m - s_lo.m == q || continue
                    amp = amplitude(s_lo, s_hi)
                    iszero(amp) && continue
                    value = uconvert(JUMP_UNIT, sqrt(a) * complex(amp))
                    if decay == :resolved
                        Lc = zeros(J, num_states, num_states)
                        Lc[i, k] = value
                        push!(jump_operators, Lc)
                        push!(jump_labels, DecayLabel(lo, hi, q, s_lo => s_hi))
                    else
                        L[i, k] = value
                        push!(components, (i, k))
                    end
                end
            end
            if decay == :grouped && !isempty(components)
                # All components of one operator must see the same frame
                # difference, or the operator would be time-dependent.
                diffs = [frame.offsets[i] - frame.offsets[k] for (i, k) in components]
                if !all(d -> isapprox(d, diffs[1]; atol=1e-9 * oneunit(d)), diffs)
                    throw(
                        ArgumentError(
                            "The per-state frame makes the grouped decay operator " *
                            "'$lo' ← '$hi' (q = $q) time-dependent; use " *
                            "decay = :resolved (dropping the coherence transfer " *
                            "between its components) or frame_kind = :levels",
                        ),
                    )
                end
                push!(jump_operators, L)
                push!(jump_labels, DecayLabel(lo, hi, q, nothing))
            end
        end
    end

    # --- Laser phase diffusion -----------------------------------------------
    for (b, beam) in enumerate(scheme.beams)
        iszero(beam.linewidth) && continue
        if !(b in frame.tree_beams)
            throw(
                ArgumentError(
                    "Beam $b has a non-zero linewidth but does not define the " *
                    "rotating frame (it beats against it); its phase noise is not " *
                    "representable as level dephasing",
                ),
            )
        end
        L = zeros(J, num_states, num_states)
        γ = uconvert(JUMP_UNIT, sqrt(beam.linewidth))
        for i in 1:num_states
            L[i, i] = complex(frame.traversals[i, b]) * γ
        end
        push!(jump_operators, L)
        push!(jump_labels, DephasingLabel(b))
    end

    LindbladModel(basis, frame, H, harmonics, jump_operators, jump_labels)
end

"""
    decay_operators(model::LindbladModel) -> (operators, labels)
    dephasing_operators(model::LindbladModel) -> (operators, labels)

Return the subsets of the model's jump operators (with their labels) that
describe spontaneous emission ([`DecayLabel`](@ref)) or laser phase diffusion
([`DephasingLabel`](@ref)).
"""
function decay_operators(model::LindbladModel)
    idx = findall(l -> l isa DecayLabel, model.jump_labels)
    model.jump_operators[idx], model.jump_labels[idx]
end

@doc (@doc decay_operators) function dephasing_operators(model::LindbladModel)
    idx = findall(l -> l isa DephasingLabel, model.jump_labels)
    model.jump_operators[idx], model.jump_labels[idx]
end

"""
    level_projector(model::LindbladModel, level)

Returns the projector (dimensionless matrix over the model basis) onto all
states of the given level — a fine-structure level covering all its hyperfine
levels for a hyperfine species.
"""
function level_projector(model::LindbladModel, level)
    spec = parse_level(level)
    num_states = length(model.basis)
    P = zeros(num_states, num_states)
    for (i, state) in enumerate(model.basis)
        matches =
            spec isa HyperfineNumberSpec ? state.level == spec :
            fine_structure(state.level) == fine_structure(spec)
        matches && (P[i, i] = 1.0)
    end
    P
end

export DecayLabel, DephasingLabel, LindbladModel, lindblad_model
export decay_operators, dephasing_operators, level_projector
