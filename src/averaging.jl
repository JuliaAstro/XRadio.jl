"""
    average(ms::MeasurementSet, by...) -> MeasurementSet
    average(ps::ProcessingSet, by...) -> ProcessingSet

`ms` with the samples sharing a label collapsed into one, along the axis each
criterion `by` labels ([`averaging_axis`](@ref)): time samples for a time
criterion, labeled as [`groupby`](@ref DimensionalData.groupby(::ProcessingSet, ::Any))
labels them, and channels for a channel criterion ([`ChannelBins`](@ref),
[`ByChannel`](@ref)). `by(ms)` gives one label per sample of its axis, or a
single label for the whole axis, and the averaged samples come in the order
their labels first occur. Several criteria apply in turn, the time criteria
first. A processing set is averaged Measurement Set by Measurement Set, keeping
their names.

```julia
average(ps, ByScan())                                # one sample per scan
average(ps, ms -> floor.(Int, times(ms) ./ 600))     # one per ten minutes
average(ms, _ -> :all)                               # one for the whole set
average(ps, ByScan(), ChannelBins(4))                # per scan, four channels to one
```

A datum contributes when it is unflagged, its weight is positive and finite and
its value is finite. Each averaged datum is the inverse-variance weighted mean
of those that contribute and its weight their sum; where none does it is `NaN`,
of weight zero and flagged. Averaging time then channels gives the same data as
averaging both at once.

Over time, `uvw` is averaged per baseline, and the time coordinate per label,
weighted by the contributing weight they carry (equally, where that is zero).
The `integration_time` of the result is the time each label spans, from its
first sample's start to its last sample's end. One Measurement Set states one
integration time, so labels spanning different lengths of time are an error:
average the groups of `groupby(ps, ByScan())` instead, or choose labels of equal
length.

Over channels, an averaged channel's frequency is the mean of its channels'
frequencies and its `channel_width` their summed widths. The channels under a
label need not be adjacent. One Measurement Set states one channel width, so
labels gathering different numbers of channels are an error.

Two labels averaging to the same time or frequency are an error. A coordinate
that is not numeric, such as `scan_name` or `field_name`, must hold one value
under each label. Any other variable stated along the averaged axis is an
error, as is more than one data group: remove them first. Sub-datasets sampled
on clocks of their own, such as the weather or the phase calibration, are kept
as they are.

The result is held in memory, whether or not `ms` is.
"""
average(ps::ProcessingSet, by, more...) = transform(ms -> average(ms, by, more...), ps)

function average(ms::MeasurementSet, by, more...)
    criteria = (by, more...)
    axes_ = map(averaging_axis, criteria)
    for (c, a) in zip(criteria, axes_)
        a in (:time, :frequency) || throw(
            ArgumentError("`averaging_axis` of $(repr(c)) is $(repr(a)); it must be :time or :frequency")
        )
    end
    order = [findall(==(:time), collect(axes_)); findall(==(:frequency), collect(axes_))]
    return foldl((m, i) -> _average(m, criteria[i], Val(axes_[i])), order; init = ms)
end

"""
    averaging_axis(by) -> Symbol

The axis the criterion `by` labels for [`average`](@ref): `:time`, the default,
so that any function of a Measurement Set is a time criterion, or `:frequency`.
A criterion type of your own that labels channels extends it.
"""
averaging_axis(by) = :time

"""
    ChannelBins(n)

The channel criterion for [`average`](@ref) that gathers each run of `n`
channels, in the order of the frequency coordinate, into one. The number of
channels must be a multiple of `n`.
"""
struct ChannelBins
    n::Int
    function ChannelBins(n::Integer)
        n > 0 || throw(ArgumentError("`ChannelBins` gathers a positive number of channels, got $n"))
        return new(n)
    end
end

averaging_axis(::ChannelBins) = :frequency
(b::ChannelBins)(ms::MeasurementSet) = (0:(length(frequencies(ms)) - 1)) .÷ b.n

"""
    ByChannel(f)

A channel criterion for [`average`](@ref) from the function `f`: `f(ms)` gives
one label per channel of `ms`, or a single label for all of them.

```julia
average(ms, ByChannel(ms -> round.(Int, frequencies(ms) ./ 8.0e6)))   # 8 MHz bins
```
"""
struct ByChannel{F}
    f::F
end

averaging_axis(::ByChannel) = :frequency
(b::ByChannel)(ms::MeasurementSet) = b.f(ms)

function _average(ms::MeasurementSet, by, ::Val{:time})
    ms = _read_one_group(ms)
    t = times(ms)
    given = by(ms)
    split = LabelLookup(
        :label, given isa AbstractArray ? _time_labels(given, ms, t) : fill(given, length(t))
    )
    members, label_of = _label_members(split, length(t))

    reduced, weights = _average_correlated(ms, Ti, label_of, length(members))
    if haskey(datagroup(ms, :base), :uvw)
        uvw_key = _group_layer(ms, :base, :uvw)
        U = ms[uvw_key]
        u = _label_mean(parent(U), _aligned_weights(weights, U), dimnum(U, Ti), members)
        reduced[uvw_key] = (U, u)
    end
    epochs, time_attrs = _averaged_epochs(ms, t, weights, members, labels(split))
    time = only(
        DimensionalData.Dimensions.format(
            (rebuild(dims(ms, Ti), Sampled(epochs; metadata = time_attrs)),), epochs
        )
    )
    ctx = (; dim = time, axis = "time", members, labels = labels(split))
    return _drop_lookups!(_average_tree(ms, ctx, reduced, ""))
end

function _average(ms::MeasurementSet, by, ::Val{:frequency})
    ms = _read_one_group(ms)
    D = dimtype(:frequency)
    f = frequencies(ms)
    split = LabelLookup(:label, _channel_labels(by(ms), f))
    members, label_of = _label_members(split, length(f))

    reduced, _ = _average_correlated(ms, D, label_of, length(members))
    centers, freq_attrs = _averaged_channels(ms, f, members, labels(split))
    freq = only(
        DimensionalData.Dimensions.format(
            (rebuild(dims(ms, D), Sampled(centers; metadata = freq_attrs)),), centers
        )
    )
    ctx = (; dim = freq, axis = "frequency", members, labels = labels(split))
    return _drop_lookups!(_average_tree(ms, ctx, reduced, ""))
end

function _read_one_group(ms::MeasurementSet)
    ms = read(ms)
    length(datagroups(ms)) == 1 || throw(
        ArgumentError(
            "`average` averages a Measurement Set with one data group; this one has " *
                join(keys(datagroups(ms)), ", ")
        )
    )
    return ms
end

# Each label's sample positions, and each sample's label number.
function _label_members(split, n::Int)
    members = [positions(split, label) for label in labels(split)]
    label_of = Vector{Int}(undef, n)
    for (l, p) in pairs(members)
        label_of[p] .= l
    end
    return members, label_of
end

function _channel_labels(given, f)
    given isa AbstractArray || return fill(given, length(f))
    if given isa DimensionalData.AbstractDimArray
        d = dims(given)
        (length(d) == 1 && only(d) isa dimtype(:frequency) && parent(val(only(d))) == f) || throw(
            ArgumentError(
                "the labels are stated over " * join(map(_dimlabel, d), ", ") *
                    ", not over the $(length(f)) channels of this Measurement Set"
            )
        )
        return collect(parent(given))
    end
    (given isa AbstractVector && length(given) == length(f)) || throw(
        ArgumentError(
            "the channel criterion gave a $(join(size(given), '×')) array; it must give a " *
                "vector of one label per channel, $(length(f)) for this Measurement Set, " *
                "or a single label for all of them"
        )
    )
    return collect(given)
end

# One node of the tree, and the sub-datasets below it. `reduced` holds layers
# already averaged by name.
function _average_tree(tree, ctx, reduced, path)
    newlayers = [
        name => _average_layer(tree, name, A, ctx, reduced, path)
            for (name, A) in pairs(DimensionalData.layers(tree))
    ]
    node = DimensionalData.rebuild_from_arrays(
        tree, newlayers; branches = DimensionalData.TreeDict()
    )
    for (name, branch) in pairs(branches(tree))
        subpath = _join(path, name)
        if _is_linked(branch, tree)
            _link!(node, name, _average_tree(branch, ctx, Dict{Symbol, Any}(), subpath), subpath)
        else
            hasdim(branch, DimensionalData.basetypeof(ctx.dim)) && throw(
                ArgumentError(
                    "sub-dataset `$subpath` states a `$(ctx.axis)` of its own rather than " *
                        "reading the Measurement Set's, so it cannot be averaged with it"
                )
            )
            branches(node)[name] = branch
        end
    end
    return node
end

function _average_layer(tree, name, A, ctx, reduced, path)
    D = DimensionalData.basetypeof(ctx.dim)
    hasdim(A, D) || return A
    haskey(reduced, name) && return _with_dim(reduced[name]..., ctx.dim)
    label = _join(path, name)
    eltype(A) <: Number && throw(
        ArgumentError(
            "`$label` is stated per $(ctx.axis), and `average` does not know how to " *
                "average it; remove it before averaging"
        )
    )
    return _with_dim(A, _constant_per_label(parent(A), dimnum(A, D), ctx, label), ctx.dim)
end

_with_dim(A, data, d) =
    rebuild(A; data, dims = DimensionalData.setdims(dims(A), d))

_resized(sz::Tuple, d::Int, n::Int) = Base.setindex(sz, n, d)

# The base data group's correlated data, flags and weights averaged along `D`,
# each as the layer it replaces and its averaged data, by layer name; and the
# contributing weight per datum.
function _average_correlated(ms, D, label_of, nlabel)
    data_key, flag_key, weight_key = (
        _group_layer(ms, :base, role) for role in (:correlated_data, :flag, :weight)
    )
    V = ms[data_key]
    F = _aligned_to(ms[flag_key], V)
    W = _aligned_to(ms[weight_key], V)
    eltype(F) == Bool || throw(
        ArgumentError(
            "`$flag_key` holds flags as bits, and `average` reads flags as true or false only"
        )
    )
    v, w, f, cw = _average_cells(parent(V), parent(F), parent(W), dimnum(V, D), label_of, nlabel)
    reduced = Dict{Symbol, Any}(data_key => (V, v), flag_key => (F, f), weight_key => (W, w))
    return reduced, DimArray(cw, dims(V))
end

# `A` with its dimensions in the order `like` holds them.
function _aligned_to(A, like)
    order = map(DimensionalData.basetypeof, dims(like))
    map(DimensionalData.basetypeof, dims(A)) == order && return A
    Set(map(DimensionalData.basetypeof, dims(A))) == Set(order) || throw(
        DimensionMismatch(
            "`$(DimensionalData.name(A))` spans " * join(map(_dimlabel, dims(A)), ", ") *
                ", not the dimensions of the data it belongs to, " *
                join(map(_dimlabel, dims(like)), ", ")
        )
    )
    return permutedims(A, dims(like))
end

# Each datum's weighted mean under its label, its summed weight, whether
# nothing contributed, and each input datum's contributing weight.
function _average_cells(v, f, w, adim::Int, label_of::Vector{Int}, nlabel::Int)
    Base.require_one_based_indexing(v, f, w)
    size(v) == size(f) == size(w) || throw(
        DimensionMismatch("data $(size(v)), flags $(size(f)) and weights $(size(w)) differ in size")
    )
    outsize = _resized(size(v), adim, nlabel)
    Tw = promote_type(Float64, eltype(w))
    num = zeros(promote_type(Float64, eltype(v)), outsize)
    den = zeros(Tw, outsize)
    cw = zeros(Tw, size(v))
    for I in CartesianIndices(v)
        vi, wi = v[I], w[I]
        (f[I] || !(wi > 0) || !isfinite(wi) || !isfinite(vi)) && continue
        J = CartesianIndex(Base.setindex(Tuple(I), label_of[I[adim]], adim))
        num[J] += wi * vi
        den[J] += wi
        cw[I] = wi
    end
    Tv = eltype(v)
    vout = map((n, d) -> iszero(d) ? _nan(Tv) : convert(Tv, n / d), num, den)
    return vout, convert(Array{eltype(w)}, den), iszero.(den), cw
end

_nan(::Type{T}) where {T <: Complex} = T(NaN, NaN)
_nan(::Type{T}) where {T} = T(NaN)

# The contributing weights summed over the dimensions `A` does not span, laid
# out to broadcast against `A`'s data.
function _aligned_weights(weights, A)
    shared = DimensionalData.commondims(A, weights)
    summed = otherdims(weights, shared)
    s = isempty(summed) ? weights : dropdims(sum(weights; dims = summed); dims = summed)
    s = permutedims(s, dims(s, shared))
    return reshape(parent(s), map(d -> hasdim(s, d) ? length(d) : 1, dims(A)))
end

# Σ wᵢaᵢ / Σ wᵢ along `tdim` over each label's samples; the equally weighted
# mean where every wᵢ is zero. A sample of zero weight adds nothing, even where
# its value is not finite.
function _label_mean(a, w, tdim::Int, members)
    out = similar(a, _resized(size(a), tdim, length(members)))
    for (l, p) in pairs(members)
        al = selectdim(a, tdim, p)
        wl = selectdim(w, tdim, p)
        num = sum(_weighted.(al, wl); dims = tdim)
        den = sum(wl; dims = tdim)
        empty = sum(al; dims = tdim) ./ length(p)
        selectdim(out, tdim, l:l) .= ifelse.(iszero.(den), empty, num ./ den)
    end
    return out
end

_weighted(a, w) = iszero(w) ? zero(a * w) : a * w

# A coordinate such as `scan_name` names something that does not change within
# an averaged sample, so each label must see one value of it.
function _constant_per_label(a, adim::Int, ctx, label)
    for (l, p) in pairs(ctx.members)
        first_ = selectdim(a, adim, first(p))
        for i in p
            isequal(selectdim(a, adim, i), first_) && continue
            seen = unique(vec(selectdim(a, adim, p)))
            throw(
                ArgumentError(
                    "`$label` takes more than one value under the label " *
                        "$(repr(ctx.labels[l])): " * join(repr.(seen), ", ") *
                        ". An averaged sample states one, so the labels must not " *
                        "mix samples that differ in `$label`"
                )
            )
        end
    end
    return collect(selectdim(a, adim, first.(ctx.members)))
end

# Each label's weighted mean epoch, and the time coordinate's attributes with
# `integration_time` the span each label covers.
function _averaged_epochs(ms, t, weights, members, names)
    attrs = metadata(lookup(ms, Ti))
    nominal = attrs[:integration_time]
    per_time = vec(_aligned_weights(weights, DimArray(t, (dims(ms, Ti),))))
    epochs = similar(t, float(eltype(t)), length(members))
    spans = similar(epochs)
    for (l, p) in pairs(members)
        t0 = t[first(p)]
        offset = t[p] .- t0
        w = per_time[p]
        total = sum(w)
        epochs[l] = t0 + (iszero(total) ? sum(offset) / length(p) : sum(w .* offset) / total)
        lo, hi = extrema(view(t, p))
        spans[l] = (hi - lo) + value(nominal)
    end
    for l in eachindex(spans)
        isapprox(spans[l], first(spans); rtol = 1.0e-9) && continue
        throw(
            ArgumentError(
                "the label $(repr(names[1])) spans $(first(spans)) s and the label " *
                    "$(repr(names[l])) spans $(spans[l]) s, but one Measurement Set " *
                    "states one `integration_time`. Split the data so each Measurement " *
                    "Set's labels span equal times first, for instance by scan with " *
                    "`groupby(ps, ByScan())`"
            )
        )
    end
    allunique(epochs) || throw(
        ArgumentError(
            "two labels average to the same epoch, and a time coordinate names each " *
                "sample once; label the samples so each label's weighted mean time differs"
        )
    )
    span = Measure(first(spans), metadata(nominal))
    return epochs, merge(attrs, Dict(:integration_time => span))
end

# Each label's mean channel frequency, and the frequency coordinate's attributes
# with `channel_width` the width each label gathers.
function _averaged_channels(ms, f, members, names)
    attrs = metadata(lookup(ms, dimtype(:frequency)))
    nominal = attrs[:channel_width]
    n = map(length, members)
    for l in eachindex(n)
        n[l] == first(n) && continue
        throw(
            ArgumentError(
                "the label $(repr(names[1])) gathers $(first(n)) channels and the label " *
                    "$(repr(names[l])) gathers $(n[l]), but one Measurement Set states one " *
                    "`channel_width`; label equal numbers of channels, for instance with " *
                    "`ChannelBins(n)` for an `n` dividing the $(length(f)) channels"
            )
        )
    end
    centers = [sum(view(f, p)) / length(p) for p in members]
    allunique(centers) || throw(
        ArgumentError(
            "two labels average to the same frequency, and a frequency coordinate names " *
                "each channel once; label the channels so each label's mean frequency differs"
        )
    )
    width = Measure(value(nominal) * first(n), metadata(nominal))
    return centers, merge(attrs, Dict(:channel_width => width))
end
