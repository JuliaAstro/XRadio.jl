"""
Selecting along a Measurement Set's own dimensions.

`DimensionalData` already indexes an `AbstractDimTree` by dimension, but three
of its defaults do not suit MSv4. Keywords resolve through `name2dim`, which
labels `baseline_id` as `Dim{:baseline_id}` and so matches nothing, where MSv4
names dimensions as the store does and that is the spelling to accept. Every
dimension is applied to every node, so a sub-dataset that legitimately carries
none of them warns. And a selector naming a single value reduces a coordinate
stated over that dimension alone to a bare element, which a tree cannot hold.
"""

const _SelectorDim = Dimension{<:Lookups.SelectorOrInterval}

# The keyword spelling is the store's own name, resolved through `dimtype` and
# so through `MSV4_DIMS`. A dimension the standard does not define keeps
# working: `dimtype` gives it `Dim{name}`, which is how a store spells it, and
# the keyword is its only spelling, since it has no type of its own.
for f in (:getindex, :view)
    for T in (:MeasurementSet, :Dataset)
        @eval Base.$f(tree::$T; kw...) = Base.$f(tree, _kw2dims(values(kw))...)

        @eval function Base.$f(tree::$T, D::_SelectorDim...)
            _require_dims(tree, D)
            return _select_tree(Base.$f, tree, D)
        end
    end
end

@inline _kw2dims(kw::NamedTuple{Keys}) where {Keys} = _kw2dims(Keys, values(kw))
@inline _kw2dims(names::Tuple, vals::Tuple) = (
    dimtype(first(names))(first(vals)),
    _kw2dims(Base.tail(names), Base.tail(vals))...,
)
_kw2dims(::Tuple{}, ::Tuple{}) = ()

# Each dimension cuts only the nodes that carry it. A sub-dataset legitimately
# has none of the correlated data's dimensions, so applying every dimension
# everywhere makes `DimensionalData` warn about nodes that are behaving
# correctly, and the warning that matters is lost among them.
function _select_tree(f, tree, D::Tuple)
    newlayers = [
        name => _select_layer(f, A, D)
            for (name, A) in pairs(DimensionalData.layers(tree))
    ]
    newbranches = DimensionalData.TreeDict(
        name => _select_tree(f, branch, D) for (name, branch) in pairs(branches(tree))
    )
    return DimensionalData.rebuild_from_arrays(tree, newlayers; branches = newbranches)
end

function _select_layer(f, A, D::Tuple)
    wanted = _carried(A, D)
    isempty(wanted) && return A
    out = f(A, wanted...)
    out isa DimensionalData.AbstractDimArray && return out
    # Reduced to a bare element. MSv4 keeps such a coordinate as a scalar, and
    # a zero-dimensional array is one; the tree holds arrays, not elements. The
    # value selected is not lost with the dimension — it is what the coordinate
    # now holds.
    return DimArray(fill(out), (); metadata = metadata(A))
end

_carried(A, D::Tuple) = filter(d -> hasdim(A, d), D)

function _require_dims(tree, D::Tuple)
    for d in D
        _anywhere(tree, d) || throw(
            ArgumentError(
                "this $(_treelabel(tree)) has no `$(_dimlabel(d))` dimension; it names " *
                    join(string.('`', sort(_dimnames(tree)), '`'), ", ")
            )
        )
    end
    return nothing
end

_anywhere(tree, d) = hasdim(tree, d) || any(b -> _anywhere(b, d), values(branches(tree)))

# Every dimension the tree names, at any depth, spelled as a store spells it
# and each reported once.
function _dimnames(tree)
    here = [_dimlabel(d) for d in dims(tree)]
    for branch in values(branches(tree))
        union!(here, _dimnames(branch))
    end
    return here
end

# The store's spelling where the standard defines one, the dimension's own name
# otherwise: a type of a caller's own has no MSv4 name, and `msv4name` throws
# rather than inventing one.
_dimlabel(d::Dimension) = _dimlabel(DimensionalData.basetypeof(d))
_dimlabel(::Type{<:Dim{S}}) where {S} = S
_dimlabel(D::Type{<:Dimension}) = get(_DIM_TO_NAME, Base.typename(D).wrapper, nameof(D))

_treelabel(::MeasurementSet) = "Measurement Set"
_treelabel(::Dataset) = "dataset"
