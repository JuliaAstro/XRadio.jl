"""
Reading MSv4 from a Zarr store.

A processing set is a Zarr group whose `type` attribute is `processing_set`.
Its child groups are the Measurement Sets, and their child groups are the
sub-datasets. Each array in a group is one labelled variable.

Both Zarr versions occur in the wild — MSv4 stores written before 2026 are
Zarr v2 — and the two record an array's dimension names differently, so
[`dimension_names`](@ref) handles both.
"""

"""
    open(ProcessingSet, store; lazy = true) -> ProcessingSet

Open the MSv4 processing set at `store`, which may be a path or an
already-opened Zarr group.

The arrays are left on disk and read on demand, so opening a large store costs
only its metadata. Pass `lazy = false`, or use `read`, to bring the data into
memory instead.

A store carrying consolidated metadata is opened from it, which is one read
rather than one per node. [`write`](@ref) produces one; a store without it opens
the slow way.

```julia
ps = open(ProcessingSet, "track.ps.zarr")
```
"""
function Base.open(::Type{ProcessingSet}, store; lazy::Bool = true)
    group = _zarr_group(store)
    type = get(group.attrs, "type", nothing)
    if type != "processing_set"
        throw(
            ArgumentError(
                "`$store` is not an MSv4 processing set: its `type` attribute is \
                 $(repr(type)), expected \"processing_set\"."
            )
        )
    end
    sets = OrderedDict{Symbol, MeasurementSet}()
    for name in sort!(collect(keys(group.groups)))
        sets[Symbol(name)] = _read_measurement_set(group.groups[name]; lazy)
    end
    return ProcessingSet(sets, _attributes(group.attrs))
end

"""
    read(ProcessingSet, store) -> ProcessingSet

Open the processing set at `store` and bring every array into memory.

Equivalent to `open(ProcessingSet, store; lazy = false)`.
"""
Base.read(::Type{ProcessingSet}, store) = open(ProcessingSet, store; lazy = false)

_zarr_group(g::Zarr.ZGroup) = g

function _zarr_group(path::AbstractString)
    store, subpath = Zarr.storefromstring(path, false)
    consolidated = _use_consolidated(store, subpath)
    return Zarr.zopen(store, "r"; path = subpath, consolidated)
end

# Whether this store's metadata should be read from its consolidated copy.
#
# Consolidated metadata is optional: a store written by another tool, or by an
# older version of this package, has none, and its absence is not an error.
#
# Only Zarr v2's `.zmetadata` is believed. Zarr v3 keeps the same information in
# the root `zarr.json`, and `Zarr.jl` reads that back as a group with no
# attributes and no children, so a v3 store must be opened the slow way until
# that is fixed.
function _use_consolidated(store::Zarr.AbstractStore, path::AbstractString)
    Zarr.ZarrFormat(store, path) isa Zarr.ZarrFormat{2} || return false
    return store[path, ".zmetadata"] !== nothing
end

"""
    dimension_names(array::Zarr.ZArray) -> Vector{Symbol}

The MSv4 dimension names of a stored array.

Zarr v2 has no place for them in its own metadata, so xarray writes them into
the `_ARRAY_DIMENSIONS` attribute. Zarr v3 has a `dimension_names` field in the
array metadata, which `Zarr.jl` accepts but does not retain, so it is read back
from the store.
"""
function dimension_names(array::Zarr.ZArray)
    names = get(array.attrs, "_ARRAY_DIMENSIONS", nothing)
    names === nothing || return Symbol.(names)
    names = _v3_dimension_names(array)
    names === nothing || return names
    throw(
        ArgumentError(
            "array `$(array.path)` declares no dimension names, so it cannot be \
             labelled. A conforming MSv4 store names every dimension."
        )
    )
end

function _v3_dimension_names(array::Zarr.ZArray)
    raw = try
        array.storage[array.path, "zarr.json"]
    catch
        nothing
    end
    raw === nothing && return nothing
    json = JSON.parse(String(copy(raw)))
    names = get(json, "dimension_names", nothing)
    return names === nothing ? nothing : Symbol.(names)
end

# Attributes that describe how an array is stored rather than what it means.
# `_ARRAY_DIMENSIONS` is where Zarr v2 keeps dimension names, `_FillValue` is
# xarray's encoding of the fill, and `coordinates` lists the non-index
# coordinates that label a variable. All three are re-derived on write, so
# keeping them in the metadata would only let a stale copy travel with the data.
const STORAGE_ATTRS = (:_ARRAY_DIMENSIONS, :_FillValue, :coordinates)

# Attribute dictionaries come back from Zarr keyed by string. MSv4 attribute
# names are identifiers, so they are more usefully symbols; nested dictionaries
# and the xarray encoding of a labelled array are decoded on the way through.
function _attributes(attrs::AbstractDict)
    out = Dict{Symbol, Any}()
    for (k, v) in attrs
        key = Symbol(k)
        key in STORAGE_ATTRS && continue
        out[key] = _attribute(v)
    end
    return out
end

function _attribute(v::AbstractDict)
    # xarray writes an array-valued attribute as a dict of data, dims and
    # attrs. Anything else is a plain nested dictionary.
    if haskey(v, "data") && haskey(v, "attrs")
        return Measure(_attribute(v["data"]), _attributes(v["attrs"]))
    end
    return _attributes(v)
end
_attribute(v::AbstractVector) = map(_attribute, v)
_attribute(v::AbstractString) = String(v)
_attribute(v) = v

# Zarr stores fixed-length strings as a padded numeric type. Nothing downstream
# wants that, so string arrays are converted on read.
_unwrap_strings(a::AbstractArray{<:Zarr.MaxLengthStrings.MaxLengthString}) = String.(a)
_unwrap_strings(a::AbstractArray) = a

function _read_measurement_set(group::Zarr.ZGroup; lazy::Bool)
    return _read_branches!(_read_into(MeasurementSet, group; lazy), group; lazy)
end

# Sub-datasets nest: a field and source dataset carries its ephemeris as a node
# of its own, and a dataset type of a caller's own may carry anything. Reading
# only the first level would drop the rest on the floor.
#
# Each is linked to its parent, as `xr.open_datatree` links it, and so has to
# agree with the parent on every dimension both state.
function _read_branches!(tree, group::Zarr.ZGroup; lazy::Bool)
    for name in sort!(collect(keys(group.groups)))
        sub = group.groups[name]
        node = _read_into(Dataset, sub; lazy, inherited = dims(tree))
        _link!(tree, nodekey(name), node, sub.path)
        _read_branches!(node, sub; lazy)
    end
    return tree
end

# Build a tree from one Zarr group: its arrays become layers, and its index
# coordinates become the dimension lookups those layers are labelled by. A
# dimension the group gives no coordinate but one of its parents does is
# labelled by the parent's, which is where xarray reads it from too.
function _read_into(::Type{T}, group::Zarr.ZGroup; lazy::Bool, inherited::Tuple = ()) where {T}
    arrays = group.arrays
    names = sort!(collect(keys(arrays)))
    dnames = Dict(name => axisnames(arrays[name]) for name in names)

    # An index coordinate is a one-dimensional array named after the dimension
    # it spans; it labels that dimension rather than being a layer of its own.
    # Index lookups are always read eagerly: they are small, and a lazy lookup
    # would fetch from disk on every selection.
    isindex(name) = dnames[name] == [Symbol(name)]
    lookups = Dict{Symbol, Any}(
        Symbol(name) => _lookup(arrays[name]) for name in names if isindex(name)
    )

    tree = T(; metadata = _attributes(group.attrs))
    for name in names
        isindex(name) && continue
        tree[julianame(name)] = _dimarray(arrays[name], dnames[name], lookups, inherited; lazy)
    end
    return tree
end

"""
    axisnames(array::Zarr.ZArray) -> Vector{Symbol}

The dimension names of a stored array, in the order its Julia axes run.

Zarr stores arrays in C order, so the first dimension named in the metadata is
the *slowest* varying one. `Zarr.jl` maps a store onto Julia's column-major
arrays without copying, which puts the fastest-varying dimension first — the
reverse of [`dimension_names`](@ref).

`XRadio` keeps that order rather than permuting back. A permutation would
either copy the whole array or destroy its contiguity, and neither is worth it
when every dimension is labelled: `vis[Frequency(At(f))]` selects the same data
whichever end of the tuple frequency sits at. It does mean `size` differs from
the shape the reference implementation reports, in exactly this reversal.
"""
axisnames(array::Zarr.ZArray) = reverse(dimension_names(array))

"""
    _lookup(array::Zarr.ZArray)

An index coordinate as a `DimensionalData` lookup, carrying its attributes.

An MSv4 coordinate's attributes are not decoration: the `time` axis holds the
integration time, and the `frequency` axis holds the spectral window's name and
channel width. Reading the values alone would silently discard them, so they
travel on the lookup, where `metadata(lookup(array, Frequency))` reaches them.
"""
function _lookup(array::Zarr.ZArray)
    values = _unwrap_strings(Array(array))
    meta = _attributes(array.attrs)
    # A categorical lookup for labels, a sampled one for anything ordered.
    return eltype(values) <: AbstractString ?
        Lookups.Categorical(values; metadata = meta) :
        Lookups.Sampled(values; metadata = meta)
end

function _dimarray(array::Zarr.ZArray, dnames, lookups, inherited = (); lazy::Bool)
    data = _unwrap_strings(lazy ? array : Array(array))
    dims = map(enumerate(dnames)) do (i, d)
        D = dimtype(d)
        haskey(lookups, d) && return D(lookups[d])
        above = DimensionalData.dims(inherited, D)
        above !== nothing && length(above) == size(data, i) && return above
        # A dimension with no index coordinate of its own — MSv4 leaves
        # `uvw_label` and the position labels implicit in some stores — is
        # numbered instead.
        return D(Base.OneTo(size(data, i)))
    end
    return DimArray(data, Tuple(dims); metadata = _attributes(array.attrs))
end
