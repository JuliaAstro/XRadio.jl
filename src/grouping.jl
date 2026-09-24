"""
    splitby(ps::ProcessingSet, f) -> ProcessingSet

`ps` with each Measurement Set cut along `time` by the labels `f(ms)` gives it,
one per integration: the integrations sharing a label become one Measurement
Set. Labels come in the order they first occur, and each piece keeps every
integration its label names, contiguous or not, so a field observed twice is
one piece.

```julia
splitby(ps, ms -> ms[:scan_name])                    # one Measurement Set per scan
splitby(ps, ms -> ms[:field_name])                   # per field
splitby(ps, ms -> floor.(Int, times(ms) ./ 600))     # per ten-minute bin
```

`f(ms)` is a vector as long as `times(ms)`, or an array over `time` alone whose
times are those of `ms`, such as a coordinate the Measurement Set states per
integration. Anything else is an error.

Each piece is a `view` of its Measurement Set, so splitting a lazily opened set
reads no data; [`read`](@ref Base.read(::Union{MeasurementSet, Dataset}))
brings one piece into memory. A piece is named after its Measurement Set and
its label, `name_label`, with the parts of a tuple label joined by `_`. Two
pieces whose names coincide are an error.

See also [`groupby`](@ref DimensionalData.groupby(::ProcessingSet, ::Any)),
which gathers Measurement Sets without cutting them.
"""
function splitby(ps::ProcessingSet, f)
    pieces = OrderedDict{Symbol, MeasurementSet}()
    for (name, ms) in pairs(ps)
        t = _time_values(ms)
        split = LabelLookup(:label, _time_labels(f(ms), ms, t))
        for label in labels(split)
            key = Symbol(name, '_', _label_name(label))
            haskey(pieces, key) && throw(
                ArgumentError(
                    "splitting gives two Measurement Sets named `$key`; " *
                        "labels whose names coincide cannot both be kept"
                )
            )
            pieces[key] = view(ms, dimtype(:time)(At(t[positions(split, label)])))
        end
    end
    return ProcessingSet(pieces, copy(metadata(ps)))
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

_time_labels(given, ms, t) = throw(
    ArgumentError(
        "`f` gave a $(nameof(typeof(given))); it must give a vector of one label " *
            "per integration, $(length(t)) for this Measurement Set"
    )
)

_label_name(label) = string(label)
_label_name(label::Tuple) = join(map(_label_name, label), '_')

"""
    groupby(ps::ProcessingSet, f) -> OrderedDict{key, ProcessingSet}

The Measurement Sets of `ps` gathered by the key `f(ms)`: each key maps to a
processing set of the Measurement Sets that give it, under their own names and
`ps`'s attributes. Keys come in the order they first occur, and so do the
Measurement Sets within a group.

```julia
groupby(ps, spectralwindow)
groupby(splitby(ps, ms -> ms[:scan_name]), ms -> (only(sources(ms)), only(scans(ms))))
```

`groupby` neither cuts a Measurement Set nor reads one; a Measurement Set
holding several scans stays whole under whatever key `f` gives it.
[`splitby`](@ref) cuts first where a key must describe each piece.
"""
function DimensionalData.groupby(ps::ProcessingSet, f)
    keys_ = [f(ms) for ms in ps]
    groups = OrderedDict{eltype(keys_), ProcessingSet}()
    for (key, (name, ms)) in zip(keys_, pairs(ps))
        group = get!(() -> ProcessingSet(OrderedDict{Symbol, MeasurementSet}(), copy(metadata(ps))), groups, key)
        getfield(group, :sets)[name] = ms
    end
    return groups
end
