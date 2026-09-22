"""
The containers an MSv4 store is read into.

MSv4 is an xarray data model: labelled arrays gathered into datasets, gathered
into a tree. `DimensionalData` provides the same three levels, so a
[`MeasurementSet`](@ref) is a `DimTree` whose layers are the correlated data
and whose branches are the sub-datasets. A [`ProcessingSet`](@ref) is the
collection of Measurement Sets a store holds; it carries no arrays of its own,
so it is an ordered collection rather than a tree.
"""

"""
    Measure(value, metadata)

A value together with the attributes needed to interpret it.

MSv4 attaches these to arrays — a frequency axis carries its channel width, a
time axis its integration time — and stores them as attributes rather than as
variables. `value` is the number or array itself; `metadata` holds its units,
reference frame, and whatever else the standard attaches.

```jldoctest
julia> m = XRadio.Measure(6.048, Dict(:units => "s", :type => "quantity"));

julia> XRadio.value(m)
6.048
```
"""
struct Measure{T, M}
    value::T
    metadata::M
end

"""
    value(m::Measure)

The number or array a [`Measure`](@ref) wraps, discarding its attributes.
"""
value(m::Measure) = m.value
value(x) = x

DimensionalData.metadata(m::Measure) = m.metadata

Base.show(io::IO, m::Measure) = print(io, "Measure(", m.value, ", ", _unit_summary(m), ")")
function _unit_summary(m::Measure)
    u = get(metadata(m), :units, nothing)
    return u === nothing ? "no units" : string(u)
end

"""
    MeasurementSet <: AbstractDimTree

One Measurement Set v4: the correlated data of a single observation together
with everything needed to interpret it.

Its layers are the correlated arrays — `:visibility` or `:spectrum`, `:flag`,
`:weight`, `:uvw` — and the auxiliary coordinates that label them but are not
dimensions, such as `:field_name` and `:scan_name`. Its branches are the
sub-datasets: `:antenna`, `:pointing`, `:weather`, and the rest. Its metadata
is the dataset's attributes, including the `data_groups` that say which flags
and weights belong to which correlated array.

Being an `AbstractDimTree`, it supports the whole of `DimensionalData`'s
selection machinery:

```julia
ms[:visibility]                       # one layer
ms[Ti(At(t))]                         # select a time
ms[Frequency(1:64), Polarization(1)]  # slice
ms.antenna                            # a sub-dataset
```

A Measurement Set covers one observation, spectral window, polarization setup,
observation mode, processor and beam per antenna. Anything coarser is several
Measurement Sets in a [`ProcessingSet`](@ref).
"""
@kwdef mutable struct MeasurementSet <: DimensionalData.AbstractDimTree
    data::DimensionalData.DataDict = DimensionalData.DataDict()
    dims::Tuple = ()
    refdims::Tuple = ()
    layerdims::DimensionalData.TupleDict = DimensionalData.TupleDict()
    layermetadata::DimensionalData.DataDict = DimensionalData.DataDict()
    metadata::Any = DimensionalData.NoMetadata()
    branches::DimensionalData.TreeDict = DimensionalData.TreeDict()
    tree::Union{Nothing, DimensionalData.AbstractDimTree} = nothing
    baseline_lookup::Union{Nothing, BaselineLookup} = nothing
    scan_lookup::Union{Nothing, LabelLookup} = nothing
    field_lookup::Union{Nothing, LabelLookup} = nothing
end

# Every lookup is built from coordinates that are themselves layers of the
# tree, so any layer that comes or goes can change them; the lookups are
# dropped and rebuilt on the next query.
function _drop_lookups!(ms::MeasurementSet)
    setfield!(ms, :baseline_lookup, nothing)
    setfield!(ms, :scan_lookup, nothing)
    setfield!(ms, :field_lookup, nothing)
    return ms
end

function Base.setindex!(ms::MeasurementSet, A::DimensionalData.AbstractDimArray, key::Symbol)
    _drop_lookups!(ms)
    return invoke(
        setindex!,
        Tuple{DimensionalData.AbstractDimTree, DimensionalData.AbstractDimArray, Symbol},
        ms, A, key
    )
end

function Base.delete!(ms::MeasurementSet, key::Symbol)
    _drop_lookups!(ms)
    return invoke(delete!, Tuple{DimensionalData.AbstractDimTree, Symbol}, ms, key)
end

function Base.pop!(ms::MeasurementSet, key::Symbol)
    _drop_lookups!(ms)
    return invoke(pop!, Tuple{DimensionalData.AbstractDimTree, Symbol}, ms, key)
end

"""
    Dataset <: AbstractDimTree

A sub-dataset of a [`MeasurementSet`](@ref) — the antennas, the pointing, the
weather, and so on. It holds labelled arrays and attributes exactly as the
Measurement Set does, but is not itself correlated data.
"""
@kwdef mutable struct Dataset <: DimensionalData.AbstractDimTree
    data::DimensionalData.DataDict = DimensionalData.DataDict()
    dims::Tuple = ()
    refdims::Tuple = ()
    layerdims::DimensionalData.TupleDict = DimensionalData.TupleDict()
    layermetadata::DimensionalData.DataDict = DimensionalData.DataDict()
    metadata::Any = DimensionalData.NoMetadata()
    branches::DimensionalData.TreeDict = DimensionalData.TreeDict()
    tree::Union{Nothing, DimensionalData.AbstractDimTree} = nothing
end

# `DimensionalData` reads a tree's layer metadata but has no `metadata` method
# for the tree itself, so its generic fallback would report `NoMetadata` and the
# dataset attributes — the whole of MSv4's metadata — would be unreachable.
DimensionalData.metadata(t::MeasurementSet) = getfield(t, :metadata)
DimensionalData.metadata(t::Dataset) = getfield(t, :metadata)

"""
    ProcessingSet

The Measurement Sets held in one store, in the order the store lists them.

A processing set carries no data of its own beyond its own attributes; it
exists so that related observations can be processed together. It behaves as an
ordered collection, indexable by position or by name:

```julia
ps = open(ProcessingSet, "track.ps.zarr")
length(ps)
ps[1]
ps[:track_00]
filter(ms -> "CALIBRATE_BANDPASS#ON_SOURCE" in scan_intents(ms), ps)
```
"""
struct ProcessingSet
    sets::OrderedDict{Symbol, MeasurementSet}
    metadata::Dict{Symbol, Any}
end

ProcessingSet(sets::OrderedDict{Symbol, MeasurementSet}) =
    ProcessingSet(sets, Dict{Symbol, Any}(:type => "processing_set"))

DimensionalData.metadata(ps::ProcessingSet) = getfield(ps, :metadata)

Base.length(ps::ProcessingSet) = length(getfield(ps, :sets))
Base.isempty(ps::ProcessingSet) = isempty(getfield(ps, :sets))
Base.keys(ps::ProcessingSet) = keys(getfield(ps, :sets))
Base.values(ps::ProcessingSet) = values(getfield(ps, :sets))
Base.pairs(ps::ProcessingSet) = pairs(getfield(ps, :sets))
Base.haskey(ps::ProcessingSet, name::Symbol) = haskey(getfield(ps, :sets), name)

Base.getindex(ps::ProcessingSet, name::Symbol) = getfield(ps, :sets)[name]
Base.getindex(ps::ProcessingSet, name::AbstractString) = ps[Symbol(name)]
Base.getindex(ps::ProcessingSet, i::Integer) = getfield(ps, :sets).vals[i]

Base.iterate(ps::ProcessingSet) = iterate(values(ps))
Base.iterate(ps::ProcessingSet, state) = iterate(values(ps), state)
Base.eltype(::Type{ProcessingSet}) = MeasurementSet

Base.IteratorSize(::Type{ProcessingSet}) = Base.HasLength()

# Tab completion over the Measurement Set names.
Base.propertynames(ps::ProcessingSet) = Tuple(keys(ps))
function Base.getproperty(ps::ProcessingSet, name::Symbol)
    haskey(getfield(ps, :sets), name) && return getfield(ps, :sets)[name]
    return getfield(ps, name)
end

"""
    filter(f, ps::ProcessingSet) -> ProcessingSet

The Measurement Sets of `ps` for which `f` returns true, keeping their names
and the set's own attributes.
"""
function Base.filter(f, ps::ProcessingSet)
    kept = OrderedDict{Symbol, MeasurementSet}(
        k => v for (k, v) in pairs(ps) if f(v)
    )
    return ProcessingSet(kept, copy(metadata(ps)))
end

"""
    merge(ps::ProcessingSet, others::ProcessingSet...) -> ProcessingSet

The Measurement Sets of every argument in one processing set, in the order the
arguments list them, under the union of their attributes.

Two kinds of disagreement are errors rather than a value quietly winning:

  - A **name** two of them share. Unlike `merge` on a dictionary, the later
    value does not win: two Measurement Sets under one name are two different
    observations, and keeping one would drop the other.
  - An **attribute** two of them state differently. A processing set's
    attributes are its store's root attributes, which `write` puts back at the
    root verbatim, so a merged set carrying one argument's would claim that
    store's attributes for every Measurement Set in the result. An attribute
    only one of them states is carried through.

```julia
merge(select(ps; field = "M87"), select(ps; field = "3C273"))
```

Sets taken from one store agree by construction, which is the case this serves.
To combine sets that genuinely disagree, state the attributes the result should
have: `ProcessingSet(sets, attrs)`.
"""
function Base.merge(ps::ProcessingSet, others::ProcessingSet...)
    sets = OrderedDict{Symbol, MeasurementSet}(pairs(ps))
    attrs = copy(metadata(ps))
    for other in others
        for (name, ms) in pairs(other)
            haskey(sets, name) && throw(
                ArgumentError(
                    "more than one of these processing sets names a Measurement " *
                        "Set `$name`; rename one before merging, since the two " *
                        "are different observations"
                )
            )
            sets[name] = ms
        end
        for (key, value) in metadata(other)
            if !haskey(attrs, key)
                attrs[key] = value
            elseif !isequal(attrs[key], value)
                throw(
                    ArgumentError(
                        "these processing sets disagree about `$key`: " *
                            "$(repr(attrs[key])) against $(repr(value)). A merged " *
                            "set carries one value per attribute, so pass the " *
                            "attributes the result should have to `ProcessingSet`."
                    )
                )
            end
        end
    end
    return ProcessingSet(sets, attrs)
end

"""
    transform(f, ps::ProcessingSet) -> ProcessingSet

`f` applied to each Measurement Set of `ps`, keeping their names and the set's
own attributes, so the result is a processing set again and composes with
[`filter`](@ref) and anything else taking one.

```julia
transform(ms -> delete!(deepcopy(ms), :weight), ps)
```

`f` must return a `MeasurementSet`, since that is all a processing set holds.
To read a value out of each Measurement Set rather than build a new one, map
over the Measurement Sets themselves: `map(gaincurves, values(ps))`.
"""
function transform(f, ps::ProcessingSet)
    made = OrderedDict{Symbol, MeasurementSet}()
    for (name, ms) in pairs(ps)
        got = f(ms)
        got isa MeasurementSet || throw(
            ArgumentError(
                "`transform` builds a ProcessingSet, so `f` must return a " *
                    "MeasurementSet; for `$name` it returned a $(nameof(typeof(got))). " *
                    "Use `map(f, values(ps))` to read a value out of each " *
                    "Measurement Set instead."
            )
        )
        made[name] = got
    end
    return ProcessingSet(made, copy(metadata(ps)))
end
