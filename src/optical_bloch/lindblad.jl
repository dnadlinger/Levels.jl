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
The Lindblad master equation of the internal states of an
[`Eigenstructure`](@ref Levels.Eigenstructure) under a [`LaserScheme`](@ref),
in the rotating frame of its [`RotatingFrame`](@ref),

```math
\\dot ρ = -i [H(t), ρ] + \\sum_k \\left( L_k ρ L_k^† - \\tfrac{1}{2} \\{ L_k^† L_k, ρ \\} \\right),
\\qquad
H(t) = H_0 + \\sum_j \\left( M_j e^{-i w_j t} + M_j^† e^{+i w_j t} \\right),
```

over the states of `basis` (the structure's basis, in its field eigenbasis):
`hamiltonian` is the static ``H_0`` (Hermitian, rad µs⁻¹) — the frame diagonal
plus the couplings ``(C + C^†)/2`` of all time-independent beam components,
with `C` from [`coupling_matrix`](@ref) — and `harmonics` the ``(w_j, M_j)``
pairs of the components left with a beat note (empty for a static model).
`couplings` keeps the complex coupling matrix ``C_b`` of each beam (rad µs⁻¹,
upper⟨row|lower⟩⟨col|, exact at the static field) as it entered the
Hamiltonian, beat-note components included — the per-beam split that the
motional layer needs, e.g. to attach each beam's displacement operator in
[`motional_model`](@ref). `jump_operators` (µs⁻¹ᐟ²) are labelled by
`jump_labels` ([`DecayLabel`](@ref), [`DephasingLabel`](@ref)).

For a hyperfine species the basis states denote the adiabatically-labelled
eigenstates at the static field, with exact at-field amplitudes throughout.
Constructed via [`lindblad_model`](@ref); `structure` and `scheme` are kept
for the motional layer.
"""
struct LindbladModel{
    S<:Eigenstructure,
    B<:StateBasis,
    F<:RotatingFrame,
    H<:Quantity,
    W<:Quantity,
    J<:Quantity,
}
    structure::S
    scheme::LaserScheme
    basis::B
    frame::F
    hamiltonian::Matrix{H}
    couplings::Vector{Matrix{H}}
    harmonics::Vector{Tuple{W,Matrix{H}}}
    jump_operators::Vector{Matrix{J}}
    jump_labels::Vector{Union{DecayLabel,DephasingLabel}}
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
    lindblad_model(s::Eigenstructure, scheme::LaserScheme; decay = :grouped, open_channels = :error)

Assembles the [`LindbladModel`](@ref) of the scheme on the structure `s` —
which fixes the species, the static field and the states kept — from the
atomic data. All matrices use the structure's field eigenbasis and ordering.

**Hamiltonian.** The frame diagonal holds each state's energy at the static
field relative to its level's zero-field centroid ([`Levels.offsets`](@ref))
minus the corresponding frame offset of [`rotating_frame`](@ref); every
beam adds ``(C + C^†)/2`` with the complex Rabi frequencies of
[`coupling_matrix`](@ref) over the two fine-structure manifolds it connects
(all their ``F`` levels for a hyperfine species) — into the static part for
its time-independent components, as harmonic terms for those with a beat
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
terms between them, and is then required). Omitting states truncates the
*model* to the retained eigenstates: an upper level with decay channels to
levels missing from the structure raises an error naming them, unless
`open_channels = :drop` discards those channels (the level then lives longer
than its physical lifetime); missing F levels or Zeeman states are treated by
the same policy.

**Laser linewidth.** A beam with non-zero `linewidth` ``γ`` (FWHM) adds the
phase-diffusion operator ``\\sqrt{γ} \\sum_i n_{b,i} P_i`` with ``P_i`` the state
projectors and ``n_{b,i}`` the signed traversal counts of the beam on the frame
tree. The coherence between states ``i`` and ``j`` decays at
``(1/2) ∑_b γ_b (n_{b,i} - n_{b,j})^2``, independently of the frame reference.
Linewidths of beams outside the
frame tree are not representable this way and raise an error.
"""
function lindblad_model(
    s::Eigenstructure,
    scheme::LaserScheme;
    decay::Symbol=:grouped,
    open_channels::Symbol=:error,
)
    decay in (:grouped, :resolved) ||
        throw(ArgumentError("decay must be :grouped or :resolved, got $decay"))
    open_channels in (:error, :drop) || throw(
        ArgumentError("open_channels must be :error or :drop, got $open_channels"),
    )

    species = s.species
    basis = s.basis
    num_states = length(basis)
    couplings = beam_couplings(s, scheme)
    frame = rotating_frame(s, scheme, couplings)
    levels = frame.levels

    # --- Hamiltonian ---------------------------------------------------------
    energies = offsets(s)
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
        rank = multipole_rank(lo => hi)
        a = einstein_a(species, lo => hi)
        amplitude = amplitude_lookup(s, lo, hi, rank)
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

    if open_channels == :error
        for (i, state) in enumerate(basis)
            lifetime_known = lifetime(species, state.level)
            isnothing(lifetime_known) && continue
            expected = 1 / lifetime_known
            actual = sum(
                sum(abs2, view(L, :, i)) for L in jump_operators;
                init=zero(expected),
            )
            if actual < (1 - 1e-10) * expected
                throw(
                    ArgumentError(
                        "Decay channels from $state leave the retained states; add the missing states or pass open_channels = :drop",
                    ),
                )
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
        # The angular FWHM doubles as the phase-diffusion rate (its rad marker is
        # dropped by the conversion, like an Einstein A's angular linewidth).
        γ = uconvert(JUMP_UNIT, sqrt(beam.linewidth))
        for i in 1:num_states
            L[i, i] = complex(frame.traversals[i, b]) * γ
        end
        push!(jump_operators, L)
        push!(jump_labels, DephasingLabel(b))
    end

    LindbladModel(
        s,
        scheme,
        basis,
        frame,
        H,
        couplings,
        harmonics,
        jump_operators,
        jump_labels,
    )
end

"""
    decay_operators(model::LindbladModel) -> (operators, labels)
    dephasing_operators(model::LindbladModel) -> (operators, labels)

Returns the subsets of the model's jump operators (with their labels) that
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

Base.show(io::IO, model::LindbladModel) = print(
    io,
    "LindbladModel(",
    length(model.basis),
    " states, ",
    length(model.couplings),
    " beams, ",
    length(model.jump_operators),
    " jumps, ",
    length(model.harmonics),
    " beat harmonics)",
)

export DecayLabel, DephasingLabel, LindbladModel, lindblad_model
export decay_operators, dephasing_operators, level_projector
