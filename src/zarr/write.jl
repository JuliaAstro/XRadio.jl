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
  Writing into an existing store keeps the store's own version.
- `compressor`: the Zarr compressor for the data arrays. Defaults to
  `Zarr.BloscCompressor(cname = "lz4", clevel = 5)`, which decompresses
  markedly faster than zstd on visibility data for a small size cost.
- `chunks = (;)`: chunk sizes for the dimensions named, by their MSv4 names,
  such as `(; baseline_id = 1)`. A name that is not a dimension of the data is
  an error.
- `chunk_bytes = 16 MiB`: the uncompressed size the dimensions `chunks` leaves
  open are fitted to; see [Chunk shapes](@ref).
- `mode = "w-"`: `"w-"` creates a store and refuses a `path` that exists;
  `"w"` replaces whatever is at `path`; `"a"` writes into the store at `path`,
  creating it if there is none. See [Writing into an existing store](@ref).
- `replace = false`: in mode `"a"`, whether a layer the store holds may be
  rewritten when `ps` holds something else there.
- `data = true`: `false` writes a template: coordinates, attributes and
  sub-datasets as usual, and each Measurement Set's data variables created
  empty, to be filled later through `open(ProcessingSet, path; mode = "r+")`.
  See [Writing a template, then filling it](@ref).
- `validate = true`: whether to check conformance before writing.
- `schemas = ()`: dataset specifications of the caller's own, for sub-dataset
  types the standard does not define. They decide both what conforms and what
  each layer is named on disk; see [`schema_types`](@ref).

```julia
write("track.ps.zarr", ps)                              # scan-sized chunks
write("track.ps.zarr", ps; chunks = (; baseline_id = 1))  # one baseline per chunk
```

To write a single Measurement Set, see
[`write(path, name => ms)`](@ref Base.write(::AbstractString, ::Pair{Symbol, MeasurementSet})).
"""
function Base.write(
        path::AbstractString, ps::ProcessingSet;
        mode::AbstractString = "w-", replace::Bool = false, kw...
    )
    meta = delete!(copy(metadata(ps)), UNREAD_ROOT_ARRAYS)
    return _write_members(path, pairs(ps), meta; mode, replace, kw...)
end

"""
    write(path, name => ms; mode = "w-", replace = false, kw...) -> path

Write the Measurement Set `ms` as the Measurement Set `name` of the store at
`path`. In `mode = "a"` this adds one Measurement Set to a store, or adds to
one it holds, leaving the store's other Measurement Sets and its own attributes
as they are; see [Writing into an existing store](@ref). In `"w-"` and `"w"` it
creates a store holding `ms` alone.

This is how calibrated visibilities join the store they were computed from,
without rewriting the visibilities they came from:

```julia
ps = open(ProcessingSet, path)
ms = copy(ps[:track_0])
ms[:visibility_corrected] = rebuild(ms[:visibility], corrected)
add_datagroup!(ms, :corrected; correlated_data = "VISIBILITY_CORRECTED")
write(path, :track_0 => ms; mode = "a")
```

The keywords are as for writing a processing set.
"""
function Base.write(
        path::AbstractString, (name, ms)::Pair{Symbol, MeasurementSet};
        mode::AbstractString = "w-", replace::Bool = false, kw...
    )
    return _write_members(path, (name => ms,), nothing; mode, replace, kw...)
end

const _WRITE_MODES = ("w-", "w", "a")

# Write the `name => ms` pairs of `members` to the store at `path`. The store's
# own attributes become `meta`; `nothing` keeps those of a store written into
# and gives a new store the bare `processing_set` type.
function _write_members(
        path::AbstractString, members, meta;
        mode::AbstractString, replace::Bool,
        zarr_format::Integer = 2,
        compressor = Zarr.BloscCompressor(cname = "lz4", clevel = 5),
        chunks::NamedTuple = (;),
        chunk_bytes::Integer = DEFAULT_CHUNK_BYTES,
        validate::Bool = true,
        schemas = (),
        data::Bool = true,
    )
    mode in _WRITE_MODES || throw(
        ArgumentError(
            "mode must be one of " * join(repr.(_WRITE_MODES), ", ") * "; got $(repr(mode))"
        )
    )
    replace && mode != "a" && throw(
        ArgumentError("`replace` applies to mode \"a\", which writes into a store; mode is $(repr(mode))")
    )
    replace && !data && throw(
        ArgumentError("`data = false` cannot replace stored layers; it only creates empty ones")
    )
    sets = OrderedDict{Symbol, MeasurementSet}(members)
    validate && XRadio.validate(
        meta === nothing ? ProcessingSet(sets) : ProcessingSet(sets, meta); schemas
    )
    chunks = ChunkRule(chunks, chunk_bytes, values(sets))
    table = schema_types(schemas)

    appending = mode == "a" && ispath(path)
    if appending
        root = _open_store(path)
        zarr_format = root.zarr_format isa Zarr.ZarrFormat{3} ? 3 : 2
        meta === nothing || Zarr.writeattrs(
            Zarr.ZarrFormat(zarr_format), root.storage, root.path, _to_attrs(meta)
        )
    else
        root = _create_store(
            path, something(meta, Dict{Symbol, Any}(:type => "processing_set"));
            zarr_format, overwrite = mode == "w"
        )
    end
    for (name, ms) in sets
        key = string(name)
        if appending && haskey(root.groups, key)
            group = root.groups[key]
            stored = _read_measurement_set(group; lazy = true)
            _append_tree(group, stored, ms; replace, zarr_format, compressor, chunks, table, data)
        else
            _write_tree(root, key, ms; zarr_format, compressor, chunks, table, data)
        end
    end
    _consolidate(root)
    return path
end

# An existing processing-set store, opened for writing.
function _open_store(path::AbstractString)
    isdir(path) || throw(ArgumentError("`$path` is not a directory store"))
    root = Zarr.zopen(path, "w")
    type = get(root.attrs, "type", nothing)
    type == "processing_set" || throw(
        ArgumentError(
            "`$path` is not an MSv4 processing set: its `type` attribute is $(repr(type))"
        )
    )
    return root
end

# Write into the stored node `group`, read back as `stored`, what `tree` holds
# and it does not.
function _append_tree(group, stored, tree; replace, zarr_format, compressor, chunks, table, data = true)
    where_ = strip(group.path, '/')
    for d in getfield(tree, :dims)
        dname = msv4name(d)
        have = dims(stored, DimensionalData.basetypeof(d))
        if have === nothing
            coord = DimArray(collect(val(d)), (d,); metadata = metadata(lookup(d)))
            _write_array(group, string(dname), coord; zarr_format, compressor, chunks = nothing)
        else
            isequal(collect(val(have)), collect(val(d))) || throw(
                ArgumentError(
                    "`$dname` of `$where_` differs from the stored one; a Measurement " *
                        "Set appended to a store must keep the store's coordinates"
                )
            )
        end
    end

    auxiliary = _auxiliary_coords(tree, table)
    allcoords = sort!(collect(values(auxiliary)))
    shape = chunk_sizes(chunks, tree)
    for key in keys(tree)
        diskname = _diskname(tree, key, table)
        write_as = n -> _write_array(
            group, n, tree[key]; zarr_format, compressor, chunks = shape,
            coords = haskey(auxiliary, key) ? String[] : allcoords,
            empty = !data && !haskey(auxiliary, key),
        )
        if !haskey(group.arrays, diskname)
            write_as(diskname)
        elseif !_is_stored_as(tree[key], group.arrays[diskname])
            replace || throw(
                ArgumentError(
                    "`$diskname` of `$where_` is not the stored array; pass " *
                        "`replace = true` to rewrite it"
                )
            )
            _replace_array(group, diskname, write_as)
        end
    end
    Zarr.writeattrs(
        Zarr.ZarrFormat(zarr_format), group.storage, group.path, _to_attrs(metadata(tree))
    )

    for (subkey, sub) in branches(tree)
        n = nodename(subkey)
        if haskey(group.groups, n)
            _append_tree(
                group.groups[n], branches(stored)[subkey], sub;
                replace, zarr_format, compressor, chunks, table
            )
        else
            _write_tree(group, n, sub; zarr_format, compressor, chunks, table)
        end
    end
    return group
end

# Whether `layer` holds the stored array `z` itself, as a lazy read of it does.
# Strings are read into memory even lazily, so a string layer is compared by
# value; they are labels, and small.
function _is_stored_as(layer, z::Zarr.ZArray)
    data = parent(layer)
    eltype(data) <: AbstractString &&
        return size(data) == size(z) && isequal(data, _unwrap_strings(Array(z)))
    data isa Zarr.ZArray || return false
    folder = _folder(data.storage)
    folder === nothing && return false
    return folder == _folder(z.storage) && strip(data.path, '/') == strip(z.path, '/')
end

_folder(s::Zarr.DirectoryStore) = realpath(s.folder)
_folder(s::Zarr.ConsolidatedStore) = _folder(s.parent)
_folder(s) = nothing

# Replace a stored array by writing the new one beside it and then moving it
# into place, so that data computed lazily from the old array is read before
# the old array is gone.
function _replace_array(group, name, write_as)
    folder = _folder(group.storage)
    folder === nothing && throw(
        ArgumentError("replacing an array is supported only in a directory store")
    )
    staged = name * ".replacing"
    ispath(joinpath(folder, group.path, staged)) &&
        rm(joinpath(folder, group.path, staged); recursive = true)
    write_as(staged)
    old = joinpath(folder, group.path, name)
    rm(old; recursive = true)
    mv(joinpath(folder, group.path, staged), old)
    return nothing
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
            ArgumentError("`$path` exists; pass `mode = \"w\"` to replace it, or `mode = \"a\"` to write into it")
        )
        rm(path; recursive = true)
    end
    store = Zarr.DirectoryStore(path)
    return Zarr.zgroup(store, "", zarr_format; attrs = _to_attrs(meta))
end

function _write_tree(
        parent::Zarr.ZGroup, name::AbstractString, tree;
        zarr_format, compressor, chunks, table = MSV4_SCHEMA_TYPES, data = true
    )
    group = Zarr.zgroup(parent, name; attrs = _to_attrs(metadata(tree)))

    # Index coordinates are dimensions of the tree rather than layers, so they
    # have to be written back out as arrays of their own, carrying the
    # attributes that travel on their lookups. A sub-dataset's inherited
    # dimensions are its parent's, and are written there.
    for d in getfield(tree, :dims)
        coord = DimArray(collect(val(d)), (d,); metadata = metadata(lookup(d)))
        _write_array(group, string(msv4name(d)), coord; zarr_format, compressor, chunks = nothing)
    end

    # Each data variable names the non-index coordinates of its dataset, which
    # is how a reader learns they are coordinates and not data of their own.
    # Promotion is dataset-wide, so every variable names all of them: listing
    # only those sharing a variable's dimensions would leave a coordinate such
    # as `polarization_type`, whose `receptor_label` no data variable spans,
    # mentioned nowhere and read back as a variable.
    auxiliary = _auxiliary_coords(tree, table)
    allcoords = sort!(collect(values(auxiliary)))
    shape = chunk_sizes(chunks, tree)
    for key in keys(tree)
        coords = haskey(auxiliary, key) ? String[] : allcoords
        _write_array(
            group, _diskname(tree, key, table), tree[key];
            zarr_format, compressor, chunks = shape, coords,
            empty = !data && !haskey(auxiliary, key),
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
        zarr_format, compressor, chunks, coords = String[], empty = false
    )
    data = empty ? parent(array) : _to_storable(parent(array), zarr_format)
    # MSv4 states dimensions slowest-varying first, the reverse of the order a
    # Julia array's axes run in.
    dnames = reverse(collect(map(msv4name, dims(array))))

    attrs = _to_attrs(metadata(array))
    # Zarr v2 has nowhere in its own metadata for dimension names, so xarray
    # puts them in an attribute; v3 has a field for them, patched in below.
    zarr_format == 2 && (attrs["_ARRAY_DIMENSIONS"] = String[string(n) for n in dnames])
    isempty(coords) || (attrs["coordinates"] = join(coords, " "))

    z = Zarr.zcreate(
        eltype(data), group, name, size(data)...;
        chunks = _chunk_shape(array, chunks),
        compressor = _isstring(eltype(data)) ? Zarr.NoCompressor() : compressor,
        attrs, zarr_format,
        (empty ? (; fill_value = _template_fill(eltype(data), name)) : (;))...,
    )
    empty || (z[axes(data)...] = data)
    zarr_format == 3 && _finalize_v3_metadata(z, dnames)
    zarr_format == 2 && empty && _finalize_v2_fill(z)
    return z
end

# What an unwritten sample of a template layer reads as: flagged, or not a
# number. Complex layers read zero: Zarr.jl cannot store a non-finite complex
# fill value.
_template_fill(::Type{Bool}, name) = true
_template_fill(::Type{T}, name) where {T <: AbstractFloat} = T(NaN)
_template_fill(::Type{T}, name) where {T <: Complex{<:AbstractFloat}} = zero(T)
_template_fill(::Type{T}, name) where {T} = throw(
    ArgumentError("a template layer must be Bool, floating-point or complex; `$name` is $T")
)

# Chunk sizes are given per dimension name, in MSv4 terms; the array's own axis
# order decides where each lands. `nothing` is one chunk.
_chunk_shape(array, ::Nothing) = size(array)
function _chunk_shape(array, sizes::AbstractDict)
    return ntuple(ndims(array)) do i
        min(get(sizes, msv4name(dims(array)[i]), size(array, i)), size(array, i))
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
    json["dimension_names"] = String[string(n) for n in dnames]

    dtype = get(json, "data_type", nothing)
    if dtype isa AbstractString && haskey(_V3_COMPLEX_NAMES, dtype)
        json["data_type"] = _V3_COMPLEX_NAMES[dtype]
        _complex_fill_as_array!(json)
    end

    z.storage[key] = Vector{UInt8}(JSON.json(json))
    return z
end

# Zarr v2 has the same complex `fill_value` defect, which only an explicit fill
# value (a template layer's) exposes.
function _finalize_v2_fill(z::Zarr.ZArray)
    key = Zarr._concatpath(z.path, ".zarray")
    json = JSON.parse(String(copy(z.storage[key])), Dict{String, Any})
    _complex_fill_as_array!(json)
    z.storage[key] = Vector{UInt8}(JSON.json(json))
    return z
end

function _complex_fill_as_array!(json)
    fill = get(json, "fill_value", nothing)
    fill isa AbstractDict && haskey(fill, "re") && haskey(fill, "im") &&
        (json["fill_value"] = [fill["re"], fill["im"]])
    return json
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
    zarr_format == 2 || return map(String, a)
    width = isempty(a) ? 1 : maximum(length, a)
    return map(Zarr.MaxLengthStrings.MaxLengthString{width, UInt32}, a)
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
