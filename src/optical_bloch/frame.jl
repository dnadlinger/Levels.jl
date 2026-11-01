# The rotating frame of a laser scheme: one frame frequency per fine-structure
# level (or, optionally, per state), fixed along a spanning tree of the beam
# graph such that the tree couplings are time-independent.

"""
The rotating frame of a [`LaserScheme`](@ref).

Every basis state ``i`` rotates at a frame frequency ``ω_i``; the rotating-frame
Hamiltonian then carries ``E_i/ħ - ω_i`` on its diagonal, and a beam of frequency
``ω`` coupling states ``i`` (lower) and ``k`` (upper) does so with the time
dependence ``e^{-i (ω - ω_k + ω_i) t}``. The frame frequencies are assigned
along a breadth-first spanning tree rooted at the scheme's `frame_reference`
so that this exponent vanishes on the tree, in one of two ways
(`LaserScheme(…; frame_kind)`):

- `:levels` (default): all Zeeman/hyperfine states of a fine-structure level
  share one frame frequency; the graph has the levels as nodes and the beams as
  edges (a beam with several polarisation components is one edge). This keeps
  every spontaneous-emission operator time-independent, whatever its grouping.
  Beams outside the tree — a second beam on an already connected level pair,
  or one closing a loop of beams — beat at the frequency mismatch.
- `:states`: the graph has the basis states as nodes and the individual
  non-vanishing coupling components of the beams as edges, so beams of
  different frequency on the same level pair can still be absorbed into the
  frame as long as they couple disjoint state pairs (the σ⁺-pump/π-probe EIT
  configuration on ``J = 1/2 → 1/2``). The Hamiltonian is then static
  whenever the component graph has no frequency-inconsistent loop. A
  spontaneous-emission operator grouping several components stays
  time-independent only if all its components see the same frame difference;
  [`lindblad_model`](@ref) checks this and asks for `decay = :resolved`
  otherwise.

`offsets` store ``ω_i`` relative to the zero-field centroid of the state's
level, so that the tabulated line centres cancel: across a tree edge the offset
grows by the beam's detuning plus the zero-field hyperfine shifts of its
reference levels. `traversals[i, b]` counts, signed (upwards positive), how
often the tree path from the reference to state ``i`` runs along beam ``b`` —
the coefficients with which each beam's frequency (and hence its phase noise)
enters the state's frame. `beats` lists the coupling components
`(beam, lower, upper, w)` left with a beat note ``w = ω - (ω_k - ω_i)``
(`w ≠ 0`), which the [`LindbladModel`](@ref) keeps as harmonic terms; states no
beam reaches keep their centroid as frame frequency.
"""
struct RotatingFrame{E<:Quantity}
    "Frame construction: `:levels` or `:states`."
    kind::Symbol

    "Fine-structure levels of the scheme, in basis order."
    levels::Vector{NoHyperfineNumberSpec}

    "The reference level (whose frame frequency is its zero-field centroid)."
    reference::NoHyperfineNumberSpec

    "Frame frequency of each basis state relative to its level's zero-field centroid."
    offsets::Vector{E}

    "Signed traversal count of each beam on the tree path to each state (state × beam)."
    traversals::Matrix{Int}

    "Indices of the beams with at least one tree edge."
    tree_beams::Vector{Int}

    "Coupling components `(beam, lower, upper, w)` with a non-zero beat frequency."
    beats::Vector{Tuple{Int,Int,Int,E}}
end

"""
    frame_offset(frame::RotatingFrame, level)
    frame_offset(frame::RotatingFrame, state::StateSpec, basis)

Returns the frame offset of a level (all its states share it in a `:levels`
frame; an error is raised if they differ) or of one state.
"""
function frame_offset(frame::RotatingFrame, basis::StateBasis, level)
    spec = fine_structure(parse_level(level))
    idx = [i for (i, s) in enumerate(basis) if fine_structure(s.level) == spec]
    isempty(idx) && throw(ArgumentError("Level '$level' is not part of the basis"))
    values = frame.offsets[idx]
    all(v -> isapprox(v, values[1]; atol=1e-9 * oneunit(v)), values) || throw(
        ArgumentError(
            "The states of level '$level' do not share one frame frequency " *
            "(per-state frame); query individual states instead",
        ),
    )
    values[1]
end

frame_offset(frame::RotatingFrame, basis::StateBasis, state::StateSpec) =
    frame.offsets[stateindex(basis, state)]

"""
Returns the angular frequency of a beam relative to the zero-field centroid
interval of its fine-structure level pair: its detuning plus the zero-field
hyperfine shifts of its reference levels (zero for a fine-structure species).
"""
function beam_centroid_offset(species, beam::LaserBeam)
    uconvert(
        ANGULAR_UNIT,
        beam.frequency.offset + hyperfine_shift(species, beam.frequency.upper) -
        hyperfine_shift(species, beam.frequency.lower),
    )
end

# The complex coupling matrix of each beam in the structure's field
# eigenbasis (rad µs⁻¹), exact at the static field.
function beam_couplings(s::Eigenstructure, scheme::LaserScheme)
    Matrix{typeof(complex(1.0) * ANGULAR_UNIT)}[
        coupling_matrix(
            s,
            beam_levels(beam)[1] => beam_levels(beam)[2],
            beam.intensity,
            beam.ε,
            beam.direction,
        ) for beam in scheme.beams
    ]
end

"""
Validates the beams and frame reference of a scheme against the structure —
every beam connects two manifolds of the structure, its reference levels are
of the structure's level kind, and the frame reference (defaulting to the
upper level of the first beam, or the first manifold) is a manifold of the
structure — returning the resolved reference level.
"""
function validate_scheme(s::Eigenstructure, scheme::LaserScheme)
    basis = s.basis
    levels = fine_structure_levels(basis)
    level_kind = eltype(basis.levels)
    for (i, beam) in enumerate(scheme.beams)
        lo, hi = beam_levels(beam)
        for (which, level) in (("lower", lo), ("upper", hi))
            if !(level in levels)
                throw(
                    ArgumentError(
                        "Beam $i: $which level '$level' is not part of the structure",
                    ),
                )
            end
        end
        if !(beam.frequency.lower isa level_kind)
            throw(
                ArgumentError(
                    "Beam $i: the reference levels must be of the structure's kind " *
                    "('$(beam.frequency.lower)' for a basis of $level_kind levels)",
                ),
            )
        end
    end
    reference = something(scheme.frame_reference, first(levels))
    if !(reference in levels)
        throw(
            ArgumentError(
                "Frame reference level '$reference' is not part of the structure",
            ),
        )
    end
    reference
end

"""
    rotating_frame(s::Eigenstructure, scheme::LaserScheme) -> RotatingFrame

Assigns the frame frequencies of the scheme on the structure `s` (cf.
[`RotatingFrame`](@ref)), validating that every beam connects two of its
manifolds.
"""
rotating_frame(s::Eigenstructure, scheme::LaserScheme) =
    rotating_frame(s, scheme, beam_couplings(s, scheme))

function rotating_frame(s::Eigenstructure, scheme::LaserScheme, couplings)
    species = s.species
    basis = s.basis
    reference = validate_scheme(s, scheme)
    num_states = length(basis)
    num_beams = length(scheme.beams)
    levels = fine_structure_levels(basis)
    offsets_beam = [beam_centroid_offset(species, beam) for beam in scheme.beams]
    E = typeof(1.0 * ANGULAR_UNIT)
    frame_offsets = zeros(E, num_states)
    traversals = zeros(Int, num_states, num_beams)
    tree_beams = Int[]
    beats = Tuple{Int,Int,Int,E}[]
    level_of = [fine_structure(s.level) for s in basis]
    # Non-vanishing coupling components of each beam, as (lower, upper) index
    # pairs (numerical residue of the geometry, e.g. 1e-17 π components of a
    # σ⁺ beam, does not count).
    components = map(couplings) do C
        scale = maximum(abs, C; init=zero(real(eltype(C))))
        [
            (i, k) for i in 1:num_states for
            k in 1:num_states if abs(C[k, i]) > 1e-12 * scale
        ]
    end

    if scheme.frame_kind == :levels
        # Spanning tree over the levels.
        level_offsets = Dict{NoHyperfineNumberSpec,E}(l => zero(E) for l in levels)
        level_traversals = Dict{NoHyperfineNumberSpec,Vector{Int}}(
            l => zeros(Int, num_beams) for l in levels
        )
        visited = Set([reference])
        queue = [reference]
        while !isempty(queue)
            level = popfirst!(queue)
            for (b, beam) in enumerate(scheme.beams)
                lo, hi = beam_levels(beam)
                if level == lo && !(hi in visited)
                    next, sign = hi, +1
                elseif level == hi && !(lo in visited)
                    next, sign = lo, -1
                else
                    continue
                end
                level_offsets[next] = level_offsets[level] + sign * offsets_beam[b]
                level_traversals[next] = copy(level_traversals[level])
                level_traversals[next][b] += sign
                push!(visited, next)
                push!(tree_beams, b)
                push!(queue, next)
            end
        end
        for i in 1:num_states
            frame_offsets[i] = level_offsets[level_of[i]]
            traversals[i, :] .= level_traversals[level_of[i]]
        end
        for (b, beam) in enumerate(scheme.beams)
            b in tree_beams && continue
            lo, hi = beam_levels(beam)
            w = offsets_beam[b] - (level_offsets[hi] - level_offsets[lo])
            iszero(w) && continue
            for (i, k) in components[b]
                push!(beats, (b, i, k, w))
            end
        end
    else
        # Spanning tree over the states, rooted at the first state of the
        # reference level.
        root = findfirst(==(reference), level_of)
        visited = falses(num_states)
        visited[root] = true
        queue = [root]
        in_tree = Set{Tuple{Int,Int,Int}}()
        while !isempty(queue)
            s = popfirst!(queue)
            for b in 1:num_beams, (i, k) in components[b]
                if s == i && !visited[k]
                    next, from, sign = k, i, +1
                elseif s == k && !visited[i]
                    next, from, sign = i, k, -1
                else
                    continue
                end
                frame_offsets[next] = frame_offsets[from] + sign * offsets_beam[b]
                traversals[next, :] .= traversals[from, :]
                traversals[next, b] += sign
                visited[next] = true
                push!(in_tree, (b, i, k))
                b in tree_beams || push!(tree_beams, b)
                push!(queue, next)
            end
        end
        for b in 1:num_beams, (i, k) in components[b]
            (b, i, k) in in_tree && continue
            w = offsets_beam[b] - (frame_offsets[k] - frame_offsets[i])
            abs(w) <= 1e-9 * abs(offsets_beam[b]) && continue
            push!(beats, (b, i, k, w))
        end
    end
    RotatingFrame(
        scheme.frame_kind,
        levels,
        reference,
        frame_offsets,
        traversals,
        tree_beams,
        beats,
    )
end

"""
Splits a beam's coupling matrix into its static part and its beat-note
components `(w, C_component)` (signed `w`), per the frame's `beats`.
"""
function split_beam_coupling(frame::RotatingFrame, b::Int, C::AbstractMatrix)
    static = copy(C)
    parts = Tuple{eltype(frame.offsets),typeof(static)}[]
    for (bb, i, k, w) in frame.beats
        bb == b || continue
        part = zero(C)
        part[k, i] = C[k, i]
        static[k, i] = zero(eltype(C))
        push!(parts, (w, part))
    end
    static, parts
end

"""
Merges harmonic terms `(w > 0, M)` at (numerically) equal frequencies.
"""
function merge_harmonics(harmonics)
    merged = similar(harmonics, 0)
    for (w, M) in harmonics
        k = findfirst(t -> isapprox(t[1], w; rtol=1e-9), merged)
        if isnothing(k)
            push!(merged, (w, M))
        else
            merged[k] = (merged[k][1], merged[k][2] .+ M)
        end
    end
    merged
end

export RotatingFrame, rotating_frame, frame_offset
