"""
An ordered set of Zeeman states drawn from the levels of a species, fixing the
index order for matrix representations of operators.

The state for a given index is available by direct indexing, and iteration yields
the states in index order. The reverse lookups are [`stateindex`](@ref) for a
single state and [`staterange`](@ref) for the index block of a whole level.
"""
Base.@kwdef struct StateBasis{L<:LevelSpec}
    "The states in index order."
    states::Vector{StateSpec{L}}

    "Reverse lookup from state to index."
    indices::Dict{StateSpec{L},Int}

    "The levels with states in the basis, in order of first appearance."
    levels::Vector{L}

    "Index ranges for the levels whose states are contiguous in the basis."
    level_ranges::Dict{L,UnitRange{Int}}
end

"""
    StateBasis(states::AbstractVector{<:StateSpec})

Creates a basis from an explicit list of states, which are indexed in the given
order.

Level specifications are converted to their canonical
[`NoHyperfineNumberSpec`](@ref) form.
"""
function StateBasis(states::AbstractVector{<:StateSpec})
    if isempty(states)
        throw(ArgumentError("At least one state is required to make up a basis"))
    end
    states = [convert(StateSpec{NoHyperfineNumberSpec}, s) for s in states]

    indices = Dict{StateSpec{NoHyperfineNumberSpec},Int}()
    for (i, state) in enumerate(states)
        j = state.level.j
        if abs(state.m) > j || !isinteger(j - state.m)
            throw(
                ArgumentError(
                    "Invalid projection m = $(state.m) for level with j = $j: $state",
                ),
            )
        end
        if haskey(indices, state)
            throw(ArgumentError("Duplicate state in basis: $state"))
        end
        indices[state] = i
    end

    levels = unique!([s.level for s in states])
    level_ranges = Dict{NoHyperfineNumberSpec,UnitRange{Int}}()
    for level in levels
        idxs = [i for (i, s) in enumerate(states) if s.level == level]
        if maximum(idxs) - minimum(idxs) + 1 == length(idxs)
            level_ranges[level] = minimum(idxs):maximum(idxs)
        end
    end

    StateBasis(; states, indices, levels, level_ranges)
end

"""
    StateBasis(levels::AbstractVector)
    StateBasis(levels...)

Creates a basis from a list of levels, where each level contributes all its
``2j + 1`` Zeeman sublevels in order of increasing `m`, with the level blocks
arranged in the given order.
"""
function StateBasis(levels::AbstractVector)
    specs = [convert(NoHyperfineNumberSpec, level) for level in levels]
    if !allunique(specs)
        throw(ArgumentError("Duplicate level in basis specification"))
    end
    StateBasis([StateSpec(level, m) for level in specs for m in (-level.j):(level.j)])
end

StateBasis(levels...) = StateBasis(collect(levels))

"""
    stateindex(basis::StateBasis, state::StateSpec)
    stateindex(basis::StateBasis, level, m)

Returns the index of the given state in the basis.

An error is raised if the state is not part of the basis.
"""
function stateindex(basis::StateBasis{L}, state::StateSpec) where {L}
    index = get(basis.indices, convert(StateSpec{L}, state), nothing)
    if isnothing(index)
        throw(ArgumentError("State $state is not part of the basis"))
    end
    index
end

function stateindex(basis::StateBasis{L}, level, m) where {L}
    stateindex(basis, StateSpec(convert(L, level), m))
end

"""
    staterange(basis::StateBasis, level)

Returns the range of basis indices of the states belonging to the given level.

An error is raised if the level has no states in the basis, or if its states are
not arranged contiguously.
"""
function staterange(basis::StateBasis{L}, level) where {L}
    spec = convert(L, level)
    range = get(basis.level_ranges, spec, nothing)
    if isnothing(range)
        if spec in basis.levels
            throw(
                ArgumentError(
                    "States of level '$level' are not contiguous in the basis",
                ),
            )
        else
            throw(ArgumentError("Level '$level' has no states in the basis"))
        end
    end
    range
end

Base.length(basis::StateBasis) = length(basis.states)
Base.iterate(basis::StateBasis) = iterate(basis.states)
Base.iterate(basis::StateBasis, state) = iterate(basis.states, state)
Base.eltype(::Type{StateBasis{L}}) where {L} = StateSpec{L}
Base.getindex(basis::StateBasis, i::Integer) = basis.states[i]
Base.firstindex(::StateBasis) = 1
Base.lastindex(basis::StateBasis) = length(basis.states)
Base.in(state::StateSpec, basis::StateBasis{L}) where {L} =
    haskey(basis.indices, convert(StateSpec{L}, state))

export StateBasis, stateindex, staterange
