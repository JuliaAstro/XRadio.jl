"""
Checking data against the MSv4 standard.

[`check`](@ref) walks a processing set, a Measurement Set or a single dataset
and reports everything that does not conform, rather than stopping at the first
problem — a store written by an older tool usually has several small faults, and
seeing them together is more useful than finding them one run at a time.
[`validate`](@ref) is the same check with the failure raised as an error.
"""

"""
    SchemaIssue

One way in which data departs from the MSv4 standard. `path` locates it within
the store and `message` says what is wrong.
"""
struct SchemaIssue
    path::String
    message::String
end

Base.show(io::IO, issue::SchemaIssue) = print(io, issue.path, ": ", issue.message)

"""
    SchemaError(issues)

Raised by [`validate`](@ref) when data does not conform to the standard.
"""
struct SchemaError <: Exception
    issues::Vector{SchemaIssue}
end

function Base.showerror(io::IO, e::SchemaError)
    n = length(e.issues)
    println(io, "SchemaError: ", n, n == 1 ? " departure" : " departures", " from the MSv4 standard:")
    for issue in e.issues
        println(io, "  ", issue)
    end
    return
end

"""
    Checker(strict, schemas)

What one [`check`](@ref) call was asked to do: how strictly to judge element
types, and which schemas a dataset's `type` attribute may name.
"""
struct Checker
    strict::Bool
    schemas::Dict{String, DatasetSpec}
end

"""
    check(x; strict = true, schemas = ()) -> Vector{SchemaIssue}

Every way in which `x` departs from the MSv4 standard, where `x` is a
[`ProcessingSet`](@ref), a [`MeasurementSet`](@ref) or a [`Dataset`](@ref). An
empty result means it conforms.

The question this answers is whether the reference implementation would accept
the data, so it is no stricter than that implementation and no looser. MSv4 is
an extensible format: a dataset may carry attributes, coordinates and data
variables the standard does not define, and a store may carry whole
sub-datasets it does not define, provided their `type` names a schema passed
in `schemas`. The one bound on an extension is the
dataset's dimension set, which is checked.

`schemas` adds dataset specifications of the caller's own, for the sub-dataset
types a store carries that the standard does not define; see
[`schema_types`](@ref).

`strict` decides how element types are judged, and has no analogue in the
reference. Under `strict = true` an array must use one of the concrete types a
conforming *store* may hold — this is the check to run before writing. Under
`strict = false` any type of the same kind passes, so a `Float32` weight or a
`BFloat16` visibility is accepted; this is the check to run on data that has
been through a calculation.

```julia
issues = check(ps)
isempty(issues) || foreach(println, issues)
```
"""
function check(ps::ProcessingSet; strict::Bool = true, schemas = ())
    checker = Checker(strict, schema_types(schemas))
    issues = SchemaIssue[]
    if get(metadata(ps), :type, nothing) != "processing_set"
        push!(
            issues,
            SchemaIssue("/", "`type` attribute is $(repr(get(metadata(ps), :type, nothing))), expected \"processing_set\"")
        )
    end
    for (name, ms) in pairs(ps)
        append!(issues, _check_tree(ms, checker, string(name), Set{Symbol}()))
    end
    return issues
end

check(
    tree::Union{MeasurementSet, Dataset};
    strict::Bool = true, schemas = (), path::AbstractString = ""
) = _check_tree(tree, Checker(strict, schema_types(schemas)), path, Set{Symbol}())

_join(path, name) = isempty(path) ? string(name) : string(path, "/", name)

# Sub-datasets nest to any depth, so every node of the tree is checked against
# the schema its own `type` names.
function _check_tree(tree, c::Checker, path::AbstractString, parent::Set{Symbol})
    issues = SchemaIssue[]
    # A node holding neither variables nor dimensions is a branch of the tree
    # rather than a dataset, and has nothing to conform to.
    if !isempty(keys(tree)) || !isempty(dims(tree))
        append!(issues, _check_dataset(tree, c, path, parent))
    end
    own = Set(msv4name(d) for d in dims(tree))
    for (name, sub) in branches(tree)
        append!(issues, _check_tree(sub, c, _join(path, name), own))
    end
    return issues
end

function _check_dataset(tree, c::Checker, path::AbstractString, parent::Set{Symbol})
    issues = SchemaIssue[]
    meta = metadata(tree)
    type = get(meta, :type, nothing)

    spec = type isa AbstractString ? get(c.schemas, type, nothing) : nothing
    if spec === nothing
        push!(issues, SchemaIssue(path, "`type` attribute is $(repr(type)), which names no known dataset"))
        return issues
    end

    append!(issues, _check_dims(tree, spec, parent, path))
    append!(issues, _check_attrs(meta, spec.attrs, c, path))
    append!(issues, _check_arrays(tree, spec, c, path))
    return issues
end

# The dataset's dimension set is the one bound MSv4 places on an extension:
# extra arrays are welcome, but only over dimensions the schema already
# defines, or ones the enclosing dataset brought with it.
function _check_dims(tree, spec::DatasetSpec, parent::Set{Symbol}, path::AbstractString)
    issues = SchemaIssue[]
    found = Set(msv4name(d) for d in dims(tree))
    for d in sort!(collect(setdiff(required_dimensions(spec), found, parent)))
        push!(issues, SchemaIssue(path, "dimension `$d` is required, but missing"))
    end
    for d in sort!(collect(setdiff(found, dimensions(spec), parent)))
        push!(
            issues,
            SchemaIssue(path, "spans `$d`, which the `$(spec.name)` schema does not define")
        )
    end
    return issues
end

# Whether `name` labels a dimension of the tree with values of its own. A
# dimension a store gave no index coordinate is numbered instead, which MSv4
# counts as a coordinate that is not there.
function _hasindex(tree, name::Symbol)
    D = dimtype(name)
    hasdim(tree, D) || return false
    return !(lookup(dims(tree, D)) isa Lookups.NoLookup)
end

# The layers one specification applies to. An array a dataset may hold several
# versions of matches every layer whose name contains the specification's,
# which is how `VISIBILITY_CORRECTED` comes to be checked against `VISIBILITY`.
function _versions(tree, a::ArraySpec)
    any(s -> s.name === :allow_multiple_versions, a.attrs) ||
        return haskey(tree, a.name) ? [a.name] : Symbol[]
    return sort!([k for k in keys(tree) if occursin(string(a.name), string(k))])
end

function _check_arrays(tree, spec::DatasetSpec, c::Checker, path::AbstractString)
    issues = SchemaIssue[]
    for a in arrays(spec)
        names = _versions(tree, a)
        if !isempty(names)
            for name in names
                apath = _join(path, name === a.name ? a.diskname : uppercase(string(name)))
                append!(issues, _check_array(tree[name], a, c, apath))
            end
        elseif _hasindex(tree, a.name)
            # An index coordinate is a dimension of the tree rather than a
            # layer. Its attributes are no less required for that — the time
            # axis's scale and the frequency axis's spectral window live there.
            apath = _join(path, a.diskname)
            d = dims(tree, dimtype(a.name))
            append!(issues, _check_attrs(metadata(lookup(d)), a.attrs, c, apath))
            append!(issues, _check_eltype(eltype(val(d)), a, c, apath))
        elseif !a.optional
            push!(issues, SchemaIssue(_join(path, a.diskname), "required, but missing"))
        end
    end
    # Coordinates and data variables the schema does not define are always
    # allowed; `_check_dims` is what bounds them.
    return issues
end

function _check_eltype(T::Type, a::ArraySpec, c::Checker, path::AbstractString)
    permitted = c.strict ? a.eltypes : relaxed_eltype(a.eltypes)
    T <: permitted && return SchemaIssue[]
    return [
        SchemaIssue(
            path, "element type is $T, expected $(c.strict ? "one of " : "")$permitted"
        ),
    ]
end

function _check_array(array, a::ArraySpec, c::Checker, path::AbstractString)
    issues = _check_eltype(eltype(array), a, c, path)

    # Arrays are read with their axes in Julia order, which is the reverse of
    # the order the schema states them in.
    found = Tuple(reverse(map(msv4name, dims(array))))
    if !any(==(found), (Tuple(d) for d in a.dims))
        expected = join(("(" * join(d, ", ") * ")" for d in a.dims), " or ")
        push!(
            issues,
            SchemaIssue(path, "spans ($(join(found, ", "))), expected $expected")
        )
    end

    append!(issues, _check_attrs(metadata(array), a.attrs, c, path))
    return issues
end

function _check_attrs(meta, specs::AbstractVector{AttrSpec}, c::Checker, path::AbstractString)
    issues = SchemaIssue[]
    for s in specs
        # An attribute present with no value is missing, not wrong. Metadata
        # that is not a dictionary at all — a lookup carrying `NoMetadata` —
        # has every attribute missing.
        v = meta isa AbstractDict ? get(meta, s.name, nothing) : nothing
        if v === nothing
            s.optional || push!(
                issues, SchemaIssue(path, "attribute `$(s.name)` is required, but missing")
            )
            continue
        end
        append!(issues, _check_value(v, s, c, path))
    end
    # Extra attributes are always allowed.
    return issues
end

# One attribute value: that it is of the type the standard states, that a
# structured value's own contents conform, and that it is one of the values the
# standard permits.
function _check_value(v, s::AttrSpec, c::Checker, path::AbstractString)
    issues = SchemaIssue[]
    wrong(expected) = [SchemaIssue(path, "attribute `$(s.name)` is a $(typeof(v)), expected $expected")]

    if s.type === :str
        v isa AbstractString || return wrong("a string")
    elseif s.type === :bool
        v isa Bool || return wrong("a boolean")
    elseif s.type === :int
        v isa Integer || return wrong("an integer")
    elseif s.type === :float
        v isa AbstractFloat || return wrong("a float")
    elseif s.type === Symbol("list[str]")
        (v isa AbstractVector && all(x -> x isa AbstractString, v)) ||
            return wrong("a list of strings")
    elseif s.type === :dict
        v isa AbstractDict || return wrong("a dictionary")
        s.nested isa DictSpec &&
            append!(issues, _check_attrs(v, s.nested.fields, c, _join(path, s.name)))
    elseif s.type === :dataarray
        v isa Measure || return wrong("an array-valued attribute")
        s.nested isa ArraySpec &&
            append!(issues, _check_measure(v, s.nested, c, _join(path, s.name)))
    else
        throw(
            ArgumentError(
                "attribute `$(s.name)` is specified as `$(s.type)`, which is not an \
                 MSv4 attribute type."
            )
        )
    end

    if s.values !== nothing && !(v in s.values)
        push!(
            issues,
            SchemaIssue(path, "attribute `$(s.name)` is $(repr(v)), expected one of $(join(s.values, ", "))")
        )
    end
    return issues
end

# An array-valued attribute carries a measure: a value together with the
# attributes needed to read it. Every one the standard defines is a scalar.
function _check_measure(m::Measure, a::ArraySpec, c::Checker, path::AbstractString)
    v = value(m)
    n = v isa AbstractArray ? ndims(v) : 0
    issues = _check_eltype(v isa AbstractArray ? eltype(v) : typeof(v), a, c, path)
    if !any(d -> length(d) == n, a.dims)
        expected = join((string(length(d)) for d in a.dims), " or ")
        push!(issues, SchemaIssue(path, "has $n dimensions, expected $expected"))
    end
    append!(issues, _check_attrs(metadata(m), a.attrs, c, path))
    return issues
end

"""
    validate(x; strict = true, schemas = ()) -> x

Return `x` if it conforms to the MSv4 standard, and throw a [`SchemaError`](@ref)
listing every departure if it does not. See [`check`](@ref) for the meaning of
`strict` and `schemas`.
"""
function validate(x; strict::Bool = true, schemas = ())
    issues = check(x; strict, schemas)
    isempty(issues) || throw(SchemaError(issues))
    return x
end

"""
    isvalid(x; strict = true, schemas = ()) -> Bool

Whether `x` conforms to the MSv4 standard. Use [`check`](@ref) to find out how
it does not.
"""
Base.isvalid(
    x::Union{ProcessingSet, MeasurementSet, Dataset}; strict::Bool = true, schemas = ()
) = isempty(check(x; strict, schemas))
