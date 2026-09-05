"""
Types describing the MSv4 schema.

A [`DatasetSpec`](@ref) records what one MSv4 dataset must look like: which
coordinates and data variables it carries, which of those are optional, what
attributes each must have, and which element types a conforming store may use.
The specs themselves live in `datasets.jl`; [`check`](@ref) compares data
against them.

MSv4 names its data variables in upper snake case (`VISIBILITY`, `UVW`) and its
coordinates in lower snake case (`time`, `baseline_id`). Those are the
identifiers written to disk. In memory `XRadio` uses the lower-case spelling
throughout, so a specification carries both: `diskname` is the standard's,
`name` is the Julia key derived from it.
"""

"""
    AttrSpec

One entry in a dataset's or array's attribute dictionary.

`values` lists the permitted values when the standard restricts them — it does
this for reference frames, time scales, and the like — and is `nothing` when
any value of `type` is allowed. `type` is the attribute's type as the standard
states it: `:str`, `:bool`, `:int`, `:float`, `:list`, `:dict`, or
`:dataarray` for an attribute that is itself a labelled array.
"""
struct AttrSpec
    name::Symbol
    type::Symbol
    values::Union{Nothing, Vector{String}}
    default::Any
    optional::Bool
    doc::String
end

function AttrSpec(
        name::Symbol, type::Symbol;
        values = nothing, default = nothing, optional::Bool = false, doc::AbstractString = ""
    )
    vals = values === nothing ? nothing : collect(String.(values))
    return AttrSpec(name, type, vals, default, optional, String(doc))
end

"""
    ArraySpec

One coordinate or data variable of a dataset.

`dims` holds the dimension tuples the array may span, in MSv4 naming. Several
are listed when the standard admits more than one shape — `SOURCE_LOCATION`,
for instance, gains a `field_name` dimension only for multi-field data.

`eltypes` is the set of element types a conforming *store* may use. It does not
constrain what an in-memory dataset may hold: see [`relaxed_eltype`](@ref) and
the `strict` keyword of [`check`](@ref).
"""
struct ArraySpec
    name::Symbol
    diskname::String
    dims::Vector{Vector{Symbol}}
    eltypes::Type
    attrs::Vector{AttrSpec}
    optional::Bool
    doc::String
end

function ArraySpec(
        diskname::AbstractString, dims, eltypes::Type;
        attrs::AbstractVector{AttrSpec} = AttrSpec[], optional::Bool = false,
        doc::AbstractString = ""
    )
    return ArraySpec(
        julianame(diskname), String(diskname), _dimsets(dims), eltypes,
        collect(attrs), optional, String(doc)
    )
end

"""
    julianame(diskname) -> Symbol

The in-memory key for the MSv4 identifier `diskname`.

MSv4's upper-snake-case data variables become lower-case Julia keys; its
coordinates are already lower case and pass through unchanged.

```jldoctest
julia> XRadio.julianame("VISIBILITY"), XRadio.julianame("baseline_id")
(:visibility, :baseline_id)
```
"""
julianame(diskname::AbstractString) = Symbol(lowercase(diskname))

# A single dimension tuple, or a collection of alternatives, normalized to a
# vector of vectors.
_dimsets(dims::Tuple{Vararg{Symbol}}) = [collect(dims)]
_dimsets(dims::AbstractVector{Symbol}) = [collect(dims)]
_dimsets(dims) = [collect(Symbol, d) for d in dims]

"""
    DatasetSpec

The schema of one MSv4 dataset.

`type` lists the values its `type` attribute may carry, which is how a dataset
found in a store identifies the specification that applies to it.
"""
struct DatasetSpec
    name::Symbol
    type::Vector{String}
    coords::Vector{ArraySpec}
    vars::Vector{ArraySpec}
    attrs::Vector{AttrSpec}
    doc::String
end

function DatasetSpec(
        name::Symbol, type;
        coords::AbstractVector{ArraySpec} = ArraySpec[],
        vars::AbstractVector{ArraySpec} = ArraySpec[],
        attrs::AbstractVector{AttrSpec} = AttrSpec[],
        doc::AbstractString = ""
    )
    types = type isa AbstractString ? [String(type)] : collect(String.(type))
    return DatasetSpec(name, types, collect(coords), collect(vars), collect(attrs), String(doc))
end

"""
    arrays(spec::DatasetSpec)

Iterate `spec`'s coordinates followed by its data variables.
"""
arrays(spec::DatasetSpec) = Iterators.flatten((spec.coords, spec.vars))

"""
    dimensions(spec::DatasetSpec) -> Vector{Symbol}

Every dimension `spec`'s coordinates and variables can span, in first-seen order.
"""
function dimensions(spec::DatasetSpec)
    out = Symbol[]
    for a in arrays(spec), ds in a.dims, d in ds
        d in out || push!(out, d)
    end
    return out
end

function Base.getindex(spec::DatasetSpec, name::Symbol)
    for a in arrays(spec)
        a.name === name && return a
    end
    throw(KeyError(name))
end

Base.haskey(spec::DatasetSpec, name::Symbol) = any(a -> a.name === name, arrays(spec))

"""
    relaxed_eltype(T::Type) -> Type

Widen a schema element-type set to every Julia type of the same kind.

The MSv4 dtypes are the concrete ones a conforming *store* may use. In memory
there is no reason to refuse a `BFloat16` visibility, a `Float32` antenna
position, or a dual number carried through a derivative, so a structural check
widens each set to its kind: reals to `AbstractFloat`, complex numbers to
`Complex{<:AbstractFloat}`, integers to `Integer`, strings to `AbstractString`.
Flag sets mixing `Bool` with unsigned bit fields keep both.

```jldoctest
julia> XRadio.relaxed_eltype(Union{Float32, Float64})
AbstractFloat

julia> XRadio.relaxed_eltype(Union{ComplexF32, ComplexF64})
Complex{<:AbstractFloat}
```
"""
function relaxed_eltype(T::Type)
    ts = Base.uniontypes(T)
    all(t -> t <: AbstractString, ts) && return AbstractString
    all(t -> t <: Complex, ts) && return Complex{<:AbstractFloat}
    all(t -> t <: AbstractFloat, ts) && return AbstractFloat
    # Bool is an Integer, so the flag case must be tested before the integer one.
    any(t -> t === Bool, ts) && return Union{Bool, Unsigned}
    all(t -> t <: Integer, ts) && return Integer
    return T
end

# MSv4 states dtypes in numpy's string spelling. Only the forms the standard
# uses appear here; anything else is a schema this version does not understand.
const NUMPY_DTYPES = Dict{String, Type}(
    "<f2" => Float16, "<f4" => Float32, "<f8" => Float64,
    "<c8" => ComplexF32, "<c16" => ComplexF64,
    "<i2" => Int16, "<i4" => Int32, "<i8" => Int64,
    "|u1" => UInt8, "<u2" => UInt16, "<u4" => UInt32, "<u8" => UInt64,
    "|b1" => Bool,
    "<U0" => String,
)

"""
    DictSpec

The schema of one of MSv4's metadata dictionaries — `observation_info`,
`processor_info`, `creator`, and the entries of `data_groups`. `fields`
reuses [`AttrSpec`](@ref), whose `optional` flag marks the keys the standard
allows to be absent.
"""
struct DictSpec
    name::Symbol
    fields::Vector{AttrSpec}
    doc::String
end

DictSpec(name::Symbol, fields::AbstractVector{AttrSpec}; doc::AbstractString = "") =
    DictSpec(name, collect(fields), String(doc))

function Base.getindex(spec::DictSpec, name::Symbol)
    for f in spec.fields
        f.name === name && return f
    end
    throw(KeyError(name))
end

Base.haskey(spec::DictSpec, name::Symbol) = any(f -> f.name === name, spec.fields)
