# The rotating frame of a laser scheme: one frame frequency per fine-structure
# level, fixed along a spanning tree of the beam graph such that all tree-beam
# couplings are time-independent.

"""
The rotating frame of a [`LaserScheme`](@ref).

Every fine-structure level ``L`` (all its Zeeman/hyperfine states alike) rotates
at one frame frequency ``ω_L``; the rotating-frame Hamiltonian then carries
``E_i/ħ - ω_{L(i)}`` on its diagonal, and a beam of frequency ``ω`` connecting
levels ``L`` and ``U`` couples them with the time dependence
``e^{-i (ω - ω_U + ω_L) t}``. The frame frequencies are chosen level by level
along a breadth-first spanning tree of the beam graph (levels as nodes, beams as
edges) rooted at the scheme's `frame_reference`, so that this exponent vanishes
for every tree beam; a beam with several polarisation components still has one
frequency, so it is one edge. Beams outside the tree — a second beam on an
already connected level pair, or one closing a loop of beams — retain a beat
note at the frequency mismatch (`beats`), which the [`LindbladModel`](@ref)
keeps as harmonic terms. Levels not reached by any beam keep their centroid as
frame frequency.

`frame_offsets` store ``ω_L`` relative to the zero-field centroid of each level,
so that the tabulated line centres cancel: for a tree beam, the offset of the
upper level exceeds that of the lower by the beam's detuning plus the zero-field
hyperfine shifts of its reference levels. `traversals[L]` counts, per beam, the
signed number of times the tree path from the reference to ``L`` runs along that
beam (upwards positive) — the coefficients with which each beam's frequency
(and hence its phase noise) enters the level's frame.
"""
struct RotatingFrame{E<:Quantity}
    "Fine-structure levels of the scheme, in basis order."
    levels::Vector{NoHyperfineNumberSpec}

    "The level whose frame frequency is its zero-field centroid."
    reference::NoHyperfineNumberSpec

    "Frame frequency of each level relative to its zero-field centroid."
    frame_offsets::Dict{NoHyperfineNumberSpec,E}

    "Indices of the beams forming the spanning tree."
    tree_beams::Vector{Int}

    "Signed traversal count per beam on the tree path to each level."
    traversals::Dict{NoHyperfineNumberSpec,Vector{Int}}

    "Beams outside the tree, with the beat frequency of their coupling, ``ω - (ω_U - ω_L)``."
    beats::Vector{Tuple{Int,E}}
end

"""
Returns the angular frequency of a beam relative to the zero-field centroid
interval of its fine-structure level pair: its detuning plus the zero-field
hyperfine shifts of its reference levels (zero for a fine-structure species).
"""
function beam_centroid_offset(species, beam::LaserBeam)
    uconvert(
        u"µs^-1",
        beam.frequency.offset + hyperfine_shift(species, beam.frequency.upper) -
        hyperfine_shift(species, beam.frequency.lower),
    )
end

"""
    rotating_frame(species, scheme::LaserScheme) -> RotatingFrame

Assigns the level frame frequencies of the scheme (cf.
[`RotatingFrame`](@ref)).
"""
function rotating_frame(species, scheme::LaserScheme)
    levels = fine_structure_levels(scheme.basis)
    nbeams = length(scheme.beams)
    offsets_beam = [beam_centroid_offset(species, beam) for beam in scheme.beams]
    E = nbeams == 0 ? typeof(1.0u"µs^-1") : eltype(offsets_beam)

    frame_offsets = Dict{NoHyperfineNumberSpec,E}(l => zero(E) for l in levels)
    traversals =
        Dict{NoHyperfineNumberSpec,Vector{Int}}(l => zeros(Int, nbeams) for l in levels)
    tree_beams = Int[]
    visited = Set([scheme.frame_reference])
    queue = [scheme.frame_reference]
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
            frame_offsets[next] = frame_offsets[level] + sign * offsets_beam[b]
            traversals[next] = copy(traversals[level])
            traversals[next][b] += sign
            push!(visited, next)
            push!(tree_beams, b)
            push!(queue, next)
        end
    end

    beats = Tuple{Int,E}[]
    for (b, beam) in enumerate(scheme.beams)
        b in tree_beams && continue
        lo, hi = beam_levels(beam)
        w = offsets_beam[b] - (frame_offsets[hi] - frame_offsets[lo])
        push!(beats, (b, w))
    end
    RotatingFrame(
        levels,
        scheme.frame_reference,
        frame_offsets,
        tree_beams,
        traversals,
        beats,
    )
end

export RotatingFrame, rotating_frame
