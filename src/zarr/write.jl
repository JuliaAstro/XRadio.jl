"""
Writing MSv4 to a Zarr store.

The store mirrors what [`open`](@ref) reads: a group whose `type` is
`processing_set`, one child group per Measurement Set, one child group per
sub-dataset, and one array per variable. Dimension names are written the way
each Zarr version expects them, and every array is checked against the standard
before anything is written.
"""

"""
    write(path, ps::ProcessingSet; kw...) -> path

Write `ps` to `path` as an MSv4 processing set, and return the path.

The data is validated against the standard first, so a store this produces
conforms or the call fails. Pass `validate = false` to skip that — useful only
when deliberately writing a store that does not conform.

# Keywords
- `zarr_format = 2`: the Zarr version to write. Version 3 is what current
  reference tooling writes and is read back correctly, but `Zarr.jl` cannot yet
  write the fixed-width strings the standard specifies for v3, so those arrays
  come out variable-length and a strict check calls their dtype wrong. Version
  2 has no such gap and is what this package defaults to until it closes.
- `compressor`: the Zarr compressor for the data arrays. Defaults to
  `Zarr.BloscCompressor(cname = "lz4", clevel = 5)`, which decompresses
  markedly faster than zstd on visibility data for a small size cost.
- `chunks`: a dimension-name-to-size mapping fixing the chunk shape. The
  default is one chunk per array, which is right for small datasets and
  wrong for large ones — a store meant for partial reads needs chunking.
- `overwrite = false`: whether to replace an existing store at `path`.
- `validate = true`: whether to check conformance before writing.
- `schemas = ()`: dataset specifications of the caller's own, for sub-dataset
  types the standard does not define. They decide both what conforms and what
  each layer is named on disk; see [`schema_types`](@ref).

```julia
write("track.ps.zarr", ps; chunks = Dict(:time => 64, :frequency => 128))
```
"""
function Base.write(
        path::AbstractString, ps::ProcessingSet;
        zarr_format::Integer = 2,
        compressor = Zarr.BloscCompressor(cname = "lz4", clevel = 5),
        chunks = Dict{Symbol, Int}(),
        overwrite::Bool = false,
        validate::Bool = true,
        schemas = (),
    )
    validate && XRadio.validate(ps; schemas)

    table = schema_types(schemas)
    root = _create_store(path, metadata(ps); zarr_format, overwrite)
    for (name, ms) in pairs(ps)
        _write_tree(root, string(name), ms; zarr_format, compressor, chunks, table)
    end
    _consolidate(root)
    return path
end

# Gather every node's metadata into one blob at the root, so that opening the
# store costs one read rather than one per node. A processing set partitioned
# per scan runs to hundreds of nodes and thousands of metadata files, which is
# what makes this worth a finalize step over the whole store rather than
# something done per Measurement Set.
#
# Consolidated metadata is not part of the Zarr v3 specification. It is written
# there anyway, because the reference implementation reads it; `open` does not
# read it back on v3, for the reason given at `_use_consolidated`.
#
# `Zarr.consolidate_metadata(::ZGroup)` always writes v2's `.zmetadata`, whatever
# the group's own format, so the format is passed explicitly.
function _consolidate(root::Zarr.ZGroup)
    if root.zarr_format isa Zarr.ZarrFormat{3}
        Zarr.consolidate_metadata(root.storage, root.path, Zarr.ZarrFormat(3))
    else
        Zarr.consolidate_metadata(root.storage, root.path)
    end
    return root
end

# The root group of a new store, holding only the processing set's own
# attributes. Measurement Sets are written into it one at a time, so a
# converter reading a file larger than memory never has to hold more than the
# partition it is building.
function _create_store(path::AbstractString, meta; zarr_format::Integer, overwrite::Bool)
    if ispath(path)
        overwrite || throw(
            ArgumentError("`$path` exists; pass `overwrite = true` to replace it.")
        )
        rm(path; recursive = true)
    end
    store = Zarr.DirectoryStore(path)
    return Zarr.zgroup(store, "", zarr_format; attrs = _to_attrs(meta))
end

function _write_tree(
        parent::Zarr.ZGroup, name::AbstractString, tree;
        zarr_format, compressor, chunks, table = MSV4_SCHEMA_TYPES
    )
    group = Zarr.zgroup(parent, name; attrs = _to_attrs(metadata(tree)))

    # Index coordinates are dimensions of the tree rather than layers, so they
    # have to be written back out as arrays of their own, carrying the
    # attributes that travel on their lookups. A sub-dataset's inherited
    # dimensions are its parent's, and are written there.
    for d in getfield(tree, :dims)
        coord = DimArray(collect(val(d)), (d,); metadata = metadata(lookup(d)))
        _write_array(group, string(msv4name(d)), coord; zarr_format, compressor, chunks)
    end

    # Each data variable names the non-index coordinates of its dataset, which
    # is how a reader learns they are coordinates and not data of their own.
    # Promotion is dataset-wide, so every variable names all of them: listing
    # only those sharing a variable's dimensions would leave a coordinate such
    # as `polarization_type`, whose `receptor_label` no data variable spans,
    # mentioned nowhere and read back as a variable.
    auxiliary = _auxiliary_coords(tree, table)
    allcoords = sort!(collect(values(auxiliary)))
    for key in keys(tree)
        coords = haskey(auxiliary, key) ? String[] : allcoords
        _write_array(
            group, _diskname(tree, key, table), tree[key];
            zarr_format, compressor, chunks, coords
        )
    end
    for (subkey, sub) in branches(tree)
        _write_tree(group, nodename(subkey), sub; zarr_format, compressor, chunks, table)
    end
    return group
end

# The layers that are coordinates rather than data variables: the schema calls
# them coordinates, but they span a dimension they are not named after, so
# `DimensionalData` holds them as layers. `field_name` and `scan_name` label
# the time axis this way.
function _auxiliary_coords(tree, table)
    spec = _spec(tree, table)
    spec === nothing && return Dict{Symbol, String}()
    return Dict{Symbol, String}(
        c.name => c.diskname for c in spec.coords if haskey(tree, c.name)
    )
end

# The specification a tree's `type` attribute names, or `nothing` when it names
# none. A dataset whose type nothing describes is still written; it just gets
# no help with its names.
function _spec(tree, table)
    type = get(metadata(tree), :type, nothing)
    return type isa AbstractString ? get(table, type, nothing) : nothing
end

# The on-disk name of a layer. The schema knows it; a variable the schema does
# not define — an extra version of the visibilities, say — is written upper
# case, which is how MSv4 spells its data variables.
function _diskname(tree, key::Symbol, table)
    spec = _spec(tree, table)
    spec !== nothing && haskey(spec, key) && return spec[key].diskname
    return uppercase(string(key))
end

function _write_array(
        group::Zarr.ZGroup, name::AbstractString, array;
        zarr_format, compressor, chunks, coords = String[]
    )
    data = _to_storable(parent(array), zarr_format)
    # MSv4 states dimensions slowest-varying first, the reverse of the order a
    # Julia array's axes run in.
    dnames = reverse(collect(map(msv4name, dims(array))))

    attrs = _to_attrs(metadata(array))
    # Zarr v2 has nowhere in its own metadata for dimension names, so xarray
    # puts them in an attribute; v3 has a field for them, patched in below.
    zarr_format == 2 && (attrs["_ARRAY_DIMENSIONS"] = string.(dnames))
    isempty(coords) || (attrs["coordinates"] = join(coords, " "))

    z = Zarr.zcreate(
        eltype(data), group, name, size(data)...;
        chunks = _chunk_shape(array, chunks),
        compressor = _isstring(eltype(data)) ? Zarr.NoCompressor() : compressor,
        attrs, zarr_format,
    )
    z[axes(data)...] = data
    zarr_format == 3 && _finalize_v3_metadata(z, dnames)
    return z
end

# Chunk sizes are given per dimension name, in MSv4 terms; the array's own axis
# order decides where each lands.
function _chunk_shape(array, chunks)
    isempty(chunks) && return size(array)
    return ntuple(ndims(array)) do i
        d = msv4name(dims(array)[i])
        min(get(chunks, d, size(array, i)), size(array, i))
    end
end

# Corrections to a freshly created Zarr v3 array's metadata.
#
# `Zarr.jl` does not expose `dimension_names` when creating an array, so it is
# added afterwards. Without it the store is not self-describing and nothing —
# this package included — can label the array on read.
#
# It also serializes complex arrays in a form nothing can read back, including
# itself: the `data_type` is the lower-cased Julia type name (`complexf32`)
# where the Zarr v3 specification and `Zarr.jl`'s own reader both say
# `complex64`, and the `fill_value` is an object where the specification says a
# two-element array. MSv4 visibilities are complex, so both are corrected here.
const _V3_COMPLEX_NAMES = Dict("complexf32" => "complex64", "complexf64" => "complex128")

function _finalize_v3_metadata(z::Zarr.ZArray, dnames)
    key = Zarr._concatpath(z.path, "zarr.json")
    json = JSON.parse(String(copy(z.storage[key])), Dict{String, Any})
    json["dimension_names"] = string.(dnames)

    dtype = get(json, "data_type", nothing)
    if dtype isa AbstractString && haskey(_V3_COMPLEX_NAMES, dtype)
        json["data_type"] = _V3_COMPLEX_NAMES[dtype]
        fill = get(json, "fill_value", nothing)
        fill isa AbstractDict && haskey(fill, "re") && haskey(fill, "im") &&
            (json["fill_value"] = [fill["re"], fill["im"]])
    end

    z.storage[key] = Vector{UInt8}(JSON.json(json))
    return z
end

"""
    _to_storable(array, zarr_format)

An array in the form its Zarr version can store.

MSv4 specifies fixed-width unicode for its string arrays, and that is what Zarr
v2 gets: the strings are padded to the longest present. `Zarr.jl` can read
fixed-width strings from a v3 store but not write them — its codec picker
selects the variable-length UTF-8 codec for a padded string type and then fails
to encode it — so v3 stores get variable-length strings instead. They read back
correctly everywhere, but a strict check against the standard will call their
dtype wrong.
"""
function _to_storable(a::AbstractArray{<:AbstractString}, zarr_format)
    zarr_format == 2 || return String.(a)
    width = isempty(a) ? 1 : maximum(length, a)
    return Zarr.MaxLengthStrings.MaxLengthString{width, UInt32}.(a)
end
_to_storable(a::AbstractArray, zarr_format) = a
_to_storable(a, zarr_format) = collect(a)

# Attributes go back out as JSON, so symbol keys become strings again and a
# `Measure` is written in the dict form xarray reads it from.
_to_attrs(meta::AbstractDict) = Dict{String, Any}(string(k) => _to_attr(v) for (k, v) in meta)
_to_attrs(::DimensionalData.NoMetadata) = Dict{String, Any}()
_to_attrs(meta) = Dict{String, Any}()

_to_attr(m::Measure) = Dict{String, Any}(
    "data" => _to_attr(value(m)),
    "attrs" => _to_attrs(metadata(m)),
    "dims" => value(m) isa AbstractArray ? String[] : String[],
)
_to_attr(v::AbstractDict) = _to_attrs(v)
_to_attr(v::AbstractVector) = map(_to_attr, v)
_to_attr(v) = v

# Blosc shuffles by element size, which is meaningless for a string, and the
# string codecs do their own encoding.
_isstring(::Type{T}) where {T} = T <: AbstractString || T <: Zarr.MaxLengthStrings.MaxLengthString
