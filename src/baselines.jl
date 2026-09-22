"""
    BaselineLookup

The baseline coordinates of a Measurement Set, arranged for lookup rather than
for storage.

MSv4 states the baselines as two antenna-name coordinates over `baseline_id`,
which answers "which baseline is this pair" only by walking them. A solver asks
once per visibility, so the accessors answer from this instead.

`ids`, `antenna1` and `antenna2` run in step: position `i` of each describes one
baseline. `bypair` maps an antenna pair, in either order, to that position, and
`byantenna` one antenna name to every position naming it, both in `baseline_id`
order. A pair the store repeats resolves to its first position, which is the
baseline the coordinates give it.

A `baseline_id` needs none of this — the `baseline_id` dimension is itself the
lookup for that, so `ids` is here to answer with, not to search.
"""
struct BaselineLookup{I, S}
    ids::Vector{I}
    antenna1::Vector{S}
    antenna2::Vector{S}
    bypair::Dict{Tuple{S, S}, Int}
    byantenna::Dict{S, Vector{Int}}
end

function BaselineLookup(
        ids::AbstractVector, antenna1::AbstractVector, antenna2::AbstractVector
    )
    length(ids) == length(antenna1) == length(antenna2) || throw(
        DimensionMismatch(
            "a Measurement Set names one antenna pair per baseline; this one has " *
                "$(length(ids)) baseline ids but $(length(antenna1)) and " *
                "$(length(antenna2)) antenna names"
        )
    )
    S = promote_type(eltype(antenna1), eltype(antenna2))
    a1 = convert(Vector{S}, antenna1)
    a2 = convert(Vector{S}, antenna2)
    id = collect(ids)

    bypair = Dict{Tuple{S, S}, Int}()
    byantenna = Dict{S, Vector{Int}}()
    for i in eachindex(id, a1, a2)
        get!(bypair, (a1[i], a2[i]), i)
        get!(bypair, (a2[i], a1[i]), i)
        push!(get!(() -> Int[], byantenna, a1[i]), i)
        a1[i] == a2[i] || push!(get!(() -> Int[], byantenna, a2[i]), i)
    end
    return BaselineLookup(id, a1, a2, bypair, byantenna)
end

Base.length(bl::BaselineLookup) = length(bl.ids)

"The antennas any baseline of `bl` names, sorted."
antenna_names(bl::BaselineLookup) = sort(unique(vcat(bl.antenna1, bl.antenna2)))
