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
    check(x; strict = true) -> Vector{SchemaIssue}

Every way in which `x` departs from the MSv4 standard, where `x` is a
[`ProcessingSet`](@ref), a [`MeasurementSet`](@ref) or a [`Dataset`](@ref). An
empty result means it conforms.

`strict` decides how element types are judged. Under `strict = true` an array
must use one of the concrete types a conforming *store* may hold — this is the
check to run before writing. Under `strict = false` any type of the same kind
passes, so a `Float32` weight or a `BFloat16` visibility is accepted; this is
the check to run on data that has been through a calculation.

```julia
issues = check(ps)
isempty(issues) || foreach(println, issues)
```
"""
function check(ps::ProcessingSet; strict::Bool = true)
    issues = SchemaIssue[]
    if get(metadata(ps), :type, nothing) != "processing_set"
        push!(
            issues,
            SchemaIssue("/", "`type` attribute is $(repr(get(metadata(ps), :type, nothing))), expected \"processing_set\"")
        )
    end
    isempty(ps) && push!(issues, SchemaIssue("/", "processing set holds no Measurement Sets"))
    for (name, ms) in pairs(ps)
        append!(issues, check(ms; strict, path = string(name)))
    end
    return issues
end

function check(ms::MeasurementSet; strict::Bool = true, path::AbstractString = "")
    issues = _check_dataset(ms, strict, path)

    type = get(metadata(ms), :type, nothing)
    if type isa AbstractString && type in CORRELATED_TYPES
        append!(issues, _check_data_groups(ms, path))
        # Only the correlated dataset is required to carry antennas.
        haskey(branches(ms), :antenna) || push!(
            issues, SchemaIssue(path, "no `antenna_xds` sub-dataset; MSv4 requires one")
        )
    else
        push!(
            issues,
            SchemaIssue(
                path,
                "`type` is $(repr(type)), expected one of $(join(sort(collect(CORRELATED_TYPES)), ", "))"
            )
        )
    end

    for (name, sub) in branches(ms)
        append!(issues, check(sub; strict, path = _join(path, name)))
    end
    return issues
end

check(ds::Dataset; strict::Bool = true, path::AbstractString = "") =
    _check_dataset(ds, strict, path)

_join(path, name) = isempty(path) ? string(name) : string(path, "/", name)

function _check_dataset(tree, strict::Bool, path::AbstractString)
    issues = SchemaIssue[]
    meta = metadata(tree)
    type = get(meta, :type, nothing)

    spec = try
        type isa AbstractString ? schema(type) : nothing
    catch
        nothing
    end
    if spec === nothing
        push!(issues, SchemaIssue(path, "`type` attribute is $(repr(type)), which names no MSv4 dataset"))
        return issues
    end

    append!(issues, _check_attrs(meta, spec.attrs, path))
    append!(issues, _check_arrays(tree, spec, strict, path))
    return issues
end

# The lookups of a tree's dimensions are index coordinates in MSv4's sense, so
# a coordinate is present if it is either a layer or a dimension of the tree.
_present(tree, name::Symbol) =
    haskey(tree, name) || any(d -> msv4name(d) === name, dims(tree))

function _check_arrays(tree, spec::DatasetSpec, strict::Bool, path::AbstractString)
    issues = SchemaIssue[]

    for a in arrays(spec)
        if !_present(tree, a.name)
            a.optional || push!(issues, SchemaIssue(_join(path, a.diskname), "required, but missing"))
            continue
        end
        apath = _join(path, a.diskname)
        if haskey(tree, a.name)
            append!(issues, _check_array(tree[a.name], a, strict, apath))
        else
            # An index coordinate is a dimension of the tree rather than a
            # layer. Its attributes are no less required for that — the time
            # axis's scale and the frequency axis's spectral window live there.
            d = dims(tree, dimtype(a.name))
            append!(issues, _check_attrs(metadata(lookup(d)), a.attrs, apath))
            if !(eltype(val(d)) <: (strict ? a.eltypes : relaxed_eltype(a.eltypes)))
                push!(
                    issues,
                    SchemaIssue(apath, "element type is $(eltype(val(d))), expected $(a.eltypes)")
                )
            end
        end
    end

    # Layers the standard does not define, except the extra versions of the
    # science arrays that data groups exist to carry.
    known = Set(a.name for a in arrays(spec))
    versioned = Set(a.name for a in arrays(spec) if any(x -> x.name === :allow_multiple_versions, a.attrs))
    for name in keys(tree)
        name in known && continue
        any(v -> startswith(string(name), string(v, "_")), versioned) && continue
        push!(
            issues,
            SchemaIssue(_join(path, string(name)), "not defined by the `$(spec.name)` schema")
        )
    end
    return issues
end

function _check_array(array, a::ArraySpec, strict::Bool, path::AbstractString)
    issues = SchemaIssue[]

    permitted = strict ? a.eltypes : relaxed_eltype(a.eltypes)
    if !(eltype(array) <: permitted)
        push!(
            issues,
            SchemaIssue(
                path,
                "element type is $(eltype(array)), expected $(strict ? "one of " : "")$permitted"
            )
        )
    end

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

    append!(issues, _check_attrs(metadata(array), a.attrs, path))
    return issues
end

function _check_attrs(meta, specs::AbstractVector{AttrSpec}, path::AbstractString)
    issues = SchemaIssue[]
    meta isa AbstractDict || return issues
    for s in specs
        if !haskey(meta, s.name)
            s.optional || push!(issues, SchemaIssue(path, "attribute `$(s.name)` is required, but missing"))
            continue
        end
        s.values === nothing && continue
        v = meta[s.name]
        v isa AbstractString || continue
        if !(v in s.values)
            push!(
                issues,
                SchemaIssue(path, "attribute `$(s.name)` is $(repr(v)), expected one of $(join(s.values, ", "))")
            )
        end
    end
    return issues
end

function _check_data_groups(ms::MeasurementSet, path::AbstractString)
    issues = SchemaIssue[]
    groups = get(metadata(ms), :data_groups, nothing)
    groups isa AbstractDict || return issues   # its absence is already reported as a missing attribute

    haskey(groups, :base) || push!(
        issues, SchemaIssue(path, "`data_groups` has no `base` group, which MSv4 requires")
    )

    for (name, group) in groups
        group isa AbstractDict || continue
        gpath = _join(path, "data_groups/$name")
        append!(issues, _check_attrs(group, DATA_GROUP.fields, gpath))
        # A group names the variables it ties together; they have to be there.
        for key in (:correlated_data, :flag, :weight, :uvw)
            haskey(group, key) || continue
            target = group[key]
            target isa AbstractString || continue
            haskey(ms, julianame(target)) || push!(
                issues, SchemaIssue(gpath, "`$key` names `$target`, which the dataset does not have")
            )
        end
        if haskey(group, :field_and_source)
            target = group[:field_and_source]
            target isa AbstractString && !haskey(branches(ms), nodekey(target)) && push!(
                issues,
                SchemaIssue(gpath, "`field_and_source` names `$target`, which is not a sub-dataset")
            )
        end
    end
    return issues
end

"""
    validate(x; strict = true) -> x

Return `x` if it conforms to the MSv4 standard, and throw a [`SchemaError`](@ref)
listing every departure if it does not. See [`check`](@ref) for the meaning of
`strict`.
"""
function validate(x; strict::Bool = true)
    issues = check(x; strict)
    isempty(issues) || throw(SchemaError(issues))
    return x
end

"""
    isvalid(x; strict = true) -> Bool

Whether `x` conforms to the MSv4 standard. Use [`check`](@ref) to find out how
it does not.
"""
Base.isvalid(x::Union{ProcessingSet, MeasurementSet, Dataset}; strict::Bool = true) =
    isempty(check(x; strict))
