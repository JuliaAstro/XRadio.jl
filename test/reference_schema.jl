# Compares the hand-written schema tables in `src/schema/` against the
# machine-readable schema exported from the reference implementation.
#
# The tables are written by hand so they read as Julia rather than as
# transcribed Python, but a hand-written table can drift from the standard
# silently. These helpers reduce both to the same normal form so that any
# difference — a missing variable, a wrong dtype, a permitted value that was
# not carried over — fails loudly.
#
# Regenerate `reference/*.json` against a new xradio release with:
#
#     python -c 'from xradio.schema.export import export_schema_json_file; \
#                import xradio.measurement_set.schema as s; \
#                export_schema_json_file(s.VisibilityXds, "VisibilityXds.json")'

using JSON

const REFERENCE_DIR = joinpath(@__DIR__, "reference")

# numpy dtype spellings, as they appear in the exported schema.
const REF_DTYPES = XRadio.NUMPY_DTYPES

"""
    reference_datasets() -> Dict{String, JSON.Object}

Every exported reference schema, keyed by the `type` attribute it declares.
"""
function reference_datasets()
    out = Dict{String, Any}()
    for file in readdir(REFERENCE_DIR; join = true)
        endswith(file, ".json") || continue
        json = JSON.parsefile(file)
        for t in reference_types(json)
            out[t] = json
        end
    end
    return out
end

"The values a reference dataset's `type` attribute may take."
function reference_types(json)
    for a in get(json, :attributes, ())
        a.name == "type" && return String.(get(a, :literal, String[]))
    end
    return String[]
end

"The element types a reference array permits, as a Julia `Union`."
function reference_eltype(arr)
    dtypes = get(arr, :dtypes, nothing)
    dtypes === nothing && return Any
    return Union{(REF_DTYPES[String(d)] for d in dtypes)...}
end

"The dimension tuples a reference array may span, as a comparable set."
reference_dims(arr) = Set(Tuple(Symbol.(g)) for g in get(arr, :dimensions, ()))

"The dimension tuples an `ArraySpec` may span, in the same normal form."
spec_dims(a::XRadio.ArraySpec) = Set(Tuple(d) for d in a.dims)

"""
    reference_attrs(node) -> Dict{Symbol, NamedTuple}

The attributes of a reference dataset or array, normalized to `(type, optional,
values)` so they can be compared with an `AttrSpec`. `values` is `nothing`
where the reference does not restrict them.
"""
function reference_attrs(node)
    out = Dict{Symbol, NamedTuple}()
    for a in get(node, :attributes, ())
        lits = get(a, :literal, nothing)
        out[Symbol(a.name)] = (
            type = Symbol(get(a, :type, "")),
            optional = Bool(get(a, :optional, false)),
            values = lits === nothing ? nothing : Set(String.(lits)),
        )
    end
    return out
end

"The same normal form, for a hand-written `AttrSpec`."
spec_attr(a::XRadio.AttrSpec) = (
    type = a.type,
    optional = a.optional,
    values = a.values === nothing ? nothing : Set(a.values),
)

"""
    reference_arrays(json) -> Dict{String, Any}

A reference dataset's coordinates and data variables, keyed by on-disk name.
"""
function reference_arrays(json)
    out = Dict{String, Any}()
    for kind in (:coordinates, :data_vars), a in get(json, kind, ())
        out[String(a.name)] = a
    end
    return out
end
