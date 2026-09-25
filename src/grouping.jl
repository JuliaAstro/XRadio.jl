"""
    groupby(ps::ProcessingSet, f) -> OrderedDict{key, ProcessingSet}

The time samples of `ps` gathered by label. `f(ms)` labels each time sample of
a Measurement Set: a vector as long as `times(ms)`, or an array over `time`
alone whose times are those of `ms`, such as a coordinate the Measurement Set
states per integration. A value that is not an array labels the whole
Measurement Set.

Each label maps to a processing set holding, for every Measurement Set with
samples under that label, those samples under the Measurement Set's own name,
with `ps`'s attributes. Labels come in the order they first occur, across the
Measurement Sets in order and along `time` within each, and a label need not
name contiguous samples, so a field observed twice is one group.

```julia
groupby(ps, ms -> ms[:scan_name])                    # per scan name
groupby(ps, ByScan())                                # per (source, scan)
groupby(ps, ms -> floor.(Int, times(ms) ./ 600))     # per ten-minute bin
groupby(ps, spectralwindow)                          # per spectral window
```

A member is a `view` of its Measurement Set, or the Measurement Set itself
where one label covers all of it, so grouping a lazily opened set reads no
data; [`read`](@ref Base.read(::Union{MeasurementSet, Dataset})) brings a
member into memory.
"""
function DimensionalData.groupby(ps::ProcessingSet, f)
    groups = OrderedDict{Any, OrderedDict{Symbol, MeasurementSet}}()
    add!(label, name, ms) =
        get!(() -> OrderedDict{Symbol, MeasurementSet}(), groups, label)[name] = ms
    for (name, ms) in pairs(ps)
        given = f(ms)
        if given isa AbstractArray
            t = _time_values(ms)
            split = LabelLookup(:label, _time_labels(given, ms, t))
            for label in labels(split)
                pos = positions(split, label)
                piece = length(pos) == length(t) ? ms : view(ms, dimtype(:time)(At(t[pos])))
                add!(label, name, piece)
            end
        else
            add!(given, name, ms)
        end
    end
    K = isempty(groups) ? Any : mapreduce(typeof, typejoin, keys(groups))
    return OrderedDict{K, ProcessingSet}(
        label => ProcessingSet(sets, copy(metadata(ps))) for (label, sets) in groups
    )
end

function _time_labels(given::DimensionalData.AbstractDimArray, ms, t)
    d = dims(given)
    (length(d) == 1 && only(d) isa dimtype(:time) && parent(val(only(d))) == t) || throw(
        ArgumentError(
            "the labels are stated over " * join(map(_dimlabel, d), ", ") *
                ", not over the $(length(t)) times of this Measurement Set"
        )
    )
    return collect(parent(given))
end

function _time_labels(given::AbstractVector, ms, t)
    length(given) == length(t) || throw(
        ArgumentError(
            "$(length(given)) labels for the $(length(t)) times of this Measurement Set; " *
                "`f` must give one label per integration"
        )
    )
    return collect(given)
end

_time_labels(given::AbstractArray, ms, t) = throw(
    ArgumentError(
        "`f` gave a $(join(size(given), '×')) array; it must give a vector of one " *
            "label per integration, $(length(t)) for this Measurement Set, or a " *
            "single label for the whole Measurement Set"
    )
)

"""
    ByScan()

The scan criterion for [`groupby`](@ref DimensionalData.groupby(::ProcessingSet, ::Any)):
`ByScan()(ms)` labels each time sample of `ms` with `(source = …, scan = …)`,
the source observed and the scan name, so `groupby(ps, ByScan())` gathers every
spectral window of a scan into one group.

The scan is the `scan_name` coordinate. The source comes from the base data
group's [`field_and_source`](@ref) dataset: its `source_name` stated per time,
or per field and matched through the `field_name` coordinate. A Measurement
Set that cannot name one source and one scan per time sample is an error.
"""
struct ByScan end

function (::ByScan)(ms::MeasurementSet)
    t = _time_values(ms)
    haskey(ms, :scan_name) || throw(
        ArgumentError("this Measurement Set states no `scan_name`, so its samples have no scan")
    )
    scan = _time_labels(ms[:scan_name], ms, t)
    source = _sources_per_time(ms, t)
    return [(source = String(s), scan = String(c)) for (s, c) in zip(source, scan)]
end

function _sources_per_time(ms::MeasurementSet, t)
    fs = field_and_source(ms)
    haskey(fs, :source_name) || throw(
        ArgumentError("the `field_and_source` dataset states no `source_name`")
    )
    stated = fs[:source_name]
    d = dims(stated)
    length(d) == 1 && only(d) isa dimtype(:field_name) || return _time_labels(stated, ms, t)
    haskey(ms, :field_name) || throw(
        ArgumentError(
            "sources are stated per field, but this Measurement Set states no " *
                "`field_name` per time"
        )
    )
    by_field = Dict(zip(String.(lookup(only(d))), String.(parent(stated))))
    return map(_time_labels(ms[:field_name], ms, t)) do field
        get(by_field, String(field)) do
            throw(
                ArgumentError(
                    "field `$field` observed at some time is not in the " *
                        "`field_and_source` dataset, which lists " *
                        join(sort!(collect(keys(by_field))), ", ")
                )
            )
        end
    end
end
