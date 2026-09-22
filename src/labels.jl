"""
    LabelLookup

A coordinate that names one label per integration, arranged for lookup rather
than for storage.

MSv4 states `scan_name` and `field_name` over `time`, so both "which scans are
these" and "which integrations is this scan" are a pass over the time axis.
This holds the answers instead: `labels` are the distinct values in
first-occurrence order, and `positions` maps each to the positions along `time`
that it covers, ascending. `kind` is the coordinate the labels came from, and
is what the errors name.

Positions rather than times, because a caller wanting times indexes the
dimension with them and one wanting neither pays for neither.
"""
struct LabelLookup{S}
    kind::Symbol
    labels::Vector{S}
    positions::Dict{S, Vector{Int}}
end

function LabelLookup(kind::Symbol, values::AbstractVector)
    S = eltype(values)
    labels = S[]
    positions = Dict{S, Vector{Int}}()
    for i in eachindex(values)
        v = values[i]
        found = get!(positions, v) do
            push!(labels, v)
            return Int[]
        end
        push!(found, i)
    end
    return LabelLookup(kind, labels, positions)
end

Base.length(ll::LabelLookup) = length(ll.labels)

"The distinct labels of `ll`, in first-occurrence order."
labels(ll::LabelLookup) = ll.labels

"""
    positions(ll::LabelLookup, label) -> Vector{Int}

Where along `time` the integrations `label` covers sit, ascending. A label `ll`
does not carry is an error naming the ones it does.
"""
function positions(ll::LabelLookup, label)
    found = get(ll.positions, label, nothing)
    found === nothing && throw(
        ArgumentError(
            "no $(ll.kind) `$label` in this Measurement Set; it names " *
                join(ll.labels, ", ")
        )
    )
    return found
end
