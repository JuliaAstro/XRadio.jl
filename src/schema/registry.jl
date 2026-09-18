"""
Lookup tables over the schemas.

A dataset in a store says which schema applies to it in its `type` attribute,
so that is the key [`schema`](@ref) dispatches on. Sub-datasets additionally
have an on-disk node name, which MSv4 suffixes with `_xds` after the xarray
type it is stored as; that suffix carries no information and is dropped in
memory.
"""

"""
    SCHEMAS

Every MSv4 dataset schema, keyed by its name.
"""
const SCHEMAS = let
    specs = (
        VISIBILITY_SCHEMA, SPECTRUM_SCHEMA, ANTENNA_SCHEMA, POINTING_SCHEMA,
        WEATHER_SCHEMA, SYSTEM_CALIBRATION_SCHEMA, GAIN_CURVE_SCHEMA,
        PHASE_CALIBRATION_SCHEMA, PHASED_ARRAY_SCHEMA, FIELD_SOURCE_SCHEMA,
        FIELD_SOURCE_EPHEMERIS_SCHEMA,
    )
    Dict{Symbol, DatasetSpec}(s.name => s for s in specs)
end

# Every `type` attribute value the standard defines, mapped to its schema. A
# schema may answer to more than one: `radiometer` data uses the visibility
# layout.
const MSV4_SCHEMA_TYPES = Dict{String, DatasetSpec}(
    t => s for s in values(SCHEMAS) for t in s.type
)

"""
    schema_types(extra) -> Dict{String, DatasetSpec}

The `type`-to-schema table [`check`](@ref) and [`write`](@ref) resolve a
dataset's `type` attribute against: the MSv4 schemas, together with the
`DatasetSpec`s in `extra`.

MSv4 is extensible at the level of whole sub-datasets. A store may carry a node
whose `type` the standard does not define, and it conforms as long as something
says what that type must look like; `extra` is where the caller says so, at the
point of use. Passing a specification whose `type` the standard already defines
replaces the standard one.

```julia
check(ps; schemas = [CALIBRATION_SCHEMA])
```
"""
function schema_types(extra)
    isempty(extra) && return MSV4_SCHEMA_TYPES
    return merge(
        MSV4_SCHEMA_TYPES,
        Dict{String, DatasetSpec}(t => s for s in extra for t in s.type)
    )
end

"""
    schema(type::AbstractString) -> DatasetSpec
    schema(name::Symbol) -> DatasetSpec

The schema for a dataset's `type` attribute, or for a schema name.

```jldoctest
julia> XRadio.schema("radiometer").name
:visibility
```
"""
schema(type::AbstractString) = get(MSV4_SCHEMA_TYPES, type) do
    throw(ArgumentError("`$type` is not an MSv4 dataset type. Known: $(join(sort(collect(keys(MSV4_SCHEMA_TYPES))), ", "))"))
end

schema(name::Symbol) = get(SCHEMAS, name) do
    throw(ArgumentError("`$name` is not an MSv4 schema. Known: $(join(sort(collect(string.(keys(SCHEMAS)))), ", "))"))
end

"""
    CORRELATED_TYPES

The `type` values marking a dataset as the correlated data of a Measurement
Set, as opposed to one of its sub-datasets.
"""
const CORRELATED_TYPES = Set(["visibility", "radiometer", "spectrum"])

"""
    nodekey(nodename) -> Symbol

The in-memory key for the sub-dataset stored under `nodename`.

MSv4 suffixes sub-dataset nodes with `_xds`, after the xarray type they are
stored as. The suffix says nothing about the data, so it is dropped in memory.

```jldoctest
julia> XRadio.nodekey("antenna_xds"), XRadio.nodekey("field_and_source_base_xds")
(:antenna, :field_and_source_base)
```
"""
function nodekey(nodename::AbstractString)
    return Symbol(endswith(nodename, "_xds") ? chop(nodename; tail = 4) : nodename)
end

"""
    nodename(key) -> String

The on-disk node name for the sub-dataset held in memory under `key` — the
inverse of [`nodekey`](@ref).

```jldoctest
julia> XRadio.nodename(:antenna)
"antenna_xds"
```
"""
nodename(key::Symbol) = string(key, "_xds")
