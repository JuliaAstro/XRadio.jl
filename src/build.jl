"""
Building MSv4 datasets from arrays.

The schema already says what each dataset holds: which coordinates and
variables, which dimensions each spans, and which attributes the standard
attaches to them. [`XRadio.subdataset`](@ref) builds a dataset from that
description, so a reader supplies values and units and nothing else — no
converter needs to know that MSv4 spells weather's time axis `time_weather`,
that a station position is written against `cartesian_pos_label`, or that a
temperature carries `units = "K"`.

Two rules of the standard's own structure make this possible. An index
coordinate is named after the one dimension it spans, so it *is* that
dimension; and a variable that admits several shapes takes the one whose
dimensions the dataset was actually given.
"""

# The attributes the standard fixes for an array. Optional ones are left out —
# `long_name` is decoration, and `ellipsoid` only applies to a geodetic
# position — and a required attribute the standard does not fix, such as a
# direction's frame, comes from the caller.
function _schema_attrs(a::ArraySpec, overrides)
    meta = Dict{Symbol, Any}()
    for attr in a.attrs
        (attr.optional || attr.default === nothing) && continue
        meta[attr.name] = attr.default
    end
    haskey(overrides, a.name) || return meta
    for (key, value) in pairs(overrides[a.name])
        meta[key] = value
    end
    return meta
end

# An index coordinate spans exactly the dimension it is named after.
_is_index_coord(a::ArraySpec) = length(a.dims) == 1 && only(a.dims) == [a.name]

# Labels are categorical and everything else is sampled, which is the
# difference between a name and a measurement.
_lookup(values::AbstractVector{<:AbstractString}, metadata) =
    Categorical(values; metadata)
_lookup(values, metadata) = Sampled(values; metadata)

# The dimensions a variable spans, in the order a Julia array holds them, which
# is the reverse of the order the standard states them in. Several shapes are
# admitted where the standard admits several — a temperature on the correlated
# data's time axis or on the weather station's own — and the one the dataset
# was given the dimensions for is the one meant.
function _layer_dims(a::ArraySpec, dimensions, name::Symbol)
    for d in a.dims
        all(n -> haskey(dimensions, n), d) &&
            return Tuple(dimensions[n] for n in Iterators.reverse(d))
    end
    shapes = join(("(" * join(d, ", ") * ")" for d in a.dims), " or ")
    throw(ArgumentError(
        "`$(a.diskname)` of the `$name` dataset spans $shapes, and the " *
            "dimensions given are $(join(sort(collect(string.(keys(dimensions)))), ", "))"
    ))
end

"""
    XRadio.subdataset(name; metadata, array_metadata, arrays...) -> Dataset

Build the MSv4 dataset `name` — `:antenna`, `:weather`, `:gain_curve` and the
rest of [`XRadio.SCHEMAS`](@ref) — from the arrays it holds.

Each keyword names a coordinate or variable of that dataset, in XRadio's
lower-case spelling, and one given as `nothing` is left out: that is how an
optional variable a file does not carry stays absent rather than being invented
as zeros. Arrays are in **Julia order**, fastest-varying dimension first, which
is the reverse of the order the standard states them in.

Dimensions, axis order and the attributes the standard fixes all come from the
schema. `metadata` adds attributes to the dataset itself and `array_metadata`
to individual arrays, for what the standard requires but does not fix — a
direction's frame, a spectral window's name.

```julia
subdataset(
    :weather;
    station_name = ["BR", "FD"],
    cartesian_pos_label = ["x", "y", "z"],
    time_weather = [1.77e9, 1.77e9 + 300],
    station_position = positions,          # (3, station)
    temperature = temperatures,            # (time, station)
    pressure = nothing,                    # not measured, so not written
)
```

The result is a plain [`Dataset`](@ref); nothing here checks it against the
standard, which is [`check`](@ref)'s work.
"""
function subdataset(name::Symbol; metadata = (;), array_metadata = (;), arrays...)
    spec = schema(name)
    return _from_schema(
        Dataset, spec, name, Dict{Symbol, Any}(:type => first(spec.type), pairs(metadata)...),
        array_metadata, arrays
    )
end

# A tree of type `T` holding `arrays`, laid out and annotated as `spec` says.
function _from_schema(::Type{T}, spec, name, metadata, array_metadata, arrays) where {T}
    specs = Dict(a.name => a for a in XRadio.arrays(spec))
    for key in keys(arrays)
        haskey(specs, key) || throw(ArgumentError(
            "`$key` is not part of the `$name` dataset, which holds " *
                join(sort(collect(string.(keys(specs)))), ", ")
        ))
    end
    given = [(specs[k], v) for (k, v) in pairs(arrays) if v !== nothing]

    dimensions = Dict{Symbol, Any}()
    for (a, values) in given
        _is_index_coord(a) || continue
        dimensions[a.name] = dimtype(a.name)(
            _lookup(values, _schema_attrs(a, array_metadata))
        )
    end

    tree = T(; metadata)
    for (a, values) in given
        _is_index_coord(a) && continue
        tree[a.name] = DimArray(
            values, _layer_dims(a, dimensions, name);
            metadata = _schema_attrs(a, array_metadata)
        )
    end
    return tree
end

"""
    XRadio.measurement_set(type = :visibility; metadata, array_metadata,
                           subdatasets, arrays...) -> MeasurementSet

Build a Measurement Set from the arrays of its correlated dataset, `type`
being `:visibility` or `:spectrum`, and the sub-datasets it holds.

Arrays are given as for [`XRadio.subdataset`](@ref): by XRadio's lower-case
names, in Julia order, `nothing` for one left out, with the dimensions and
fixed attributes taken from the schema. `subdatasets` is a `NamedTuple` of
datasets, such as `(antenna = ..., field_and_source_base = ...)`, each linked
to the Measurement Set as assigning it would.

The Measurement Set's own attributes default to what needs no knowledge of the
data: its `type`, `schema_version`, `creation_date` (now), `creator`
(XRadio.jl), and a `base` data group naming the correlated data, flag, weight
and UVW given and the `field_and_source_base` sub-dataset if there is one.
`metadata` adds to or replaces these, and must supply what the standard
requires and the data alone does not say: `observation_info` and
`processor_info`. Coordinates' required attributes the standard does not fix,
such as `time`'s `integration_time` and the spectral window's name, come
through `array_metadata`.

```julia
ms = XRadio.measurement_set(;
    time = times, frequency = freqs, polarization = ["RR", "LL"],
    baseline_id = 0:2, baseline_antenna1_name = a1, baseline_antenna2_name = a2,
    uvw_label = ["u", "v", "w"],
    visibility = vis, flag = flags, weight = weights, uvw = uvws,   # Julia order
    field_name = fill("f0", length(times)), scan_name = fill("1", length(times)),
    metadata = (; observation_info, processor_info),
    array_metadata = (; time = (; integration_time), frequency = spw_attrs,
                        uvw = (; frame = "icrs")),
    subdatasets = (; antenna, field_and_source_base),
)
```

Nothing here checks the result against the standard; that is
[`check`](@ref)'s work.
"""
function measurement_set(
        type::Symbol = :visibility;
        metadata = (;), array_metadata = (;), subdatasets = (;), arrays...
    )
    type in (:visibility, :spectrum) || throw(ArgumentError(
        "a Measurement Set's correlated dataset is `:visibility` or `:spectrum`; got `$type`"
    ))
    spec = schema(type)
    given(key) = get(arrays, key, nothing) !== nothing
    now_ = _iso8601(now(UTC))
    base = Dict{Symbol, Any}(:description => "", :date => now_)
    for (role, key) in ((:correlated_data, type), (:flag, :flag), (:weight, :weight), (:uvw, :uvw))
        given(key) && (base[role] = spec[key].diskname)
    end
    haskey(subdatasets, :field_and_source_base) &&
        (base[:field_and_source] = nodename(:field_and_source_base))
    defaults = Dict{Symbol, Any}(
        :type => first(spec.type),
        :schema_version => string(MSV4_SCHEMA_VERSION),
        :creation_date => now_,
        :creator => Dict{Symbol, Any}(
            :software_name => "XRadio.jl", :version => string(pkgversion(XRadio)),
        ),
        :data_groups => Dict{Symbol, Any}(:base => base),
    )
    ms = _from_schema(
        MeasurementSet, spec, type, merge(defaults, Dict{Symbol, Any}(pairs(metadata))),
        array_metadata, arrays
    )
    for (key, xds) in pairs(subdatasets)
        setproperty!(ms, key, xds)
    end
    return ms
end
