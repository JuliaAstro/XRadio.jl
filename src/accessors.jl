"""
Reaching the parts of a Measurement Set that MSv4 buries in attributes.

MSv4 keeps a good deal of what you actually want to know — which antennas, which
scans, which correlated array a data group points at — inside nested attribute
dictionaries. These accessors take it out again.
"""

"""
    datagroups(ms) -> Dict{Symbol, Any}

The data groups of `ms`: the named sets of correlated data, flags, weights and
baseline coordinates that belong together. A group named `base` is always
present in conforming data.

See also [`datagroup`](@ref).
"""
datagroups(ms::MeasurementSet) = get(metadata(ms), :data_groups, Dict{Symbol, Any}())

"""
    datagroup(ms, name = :base) -> Dict{Symbol, Any}

One data group of `ms` by name, defaulting to the mandatory `base` group.

```julia
datagroup(ms)[:correlated_data]   # "VISIBILITY"
```
"""
function datagroup(ms::MeasurementSet, name::Symbol = :base)
    groups = datagroups(ms)
    haskey(groups, name) && return groups[name]
    throw(
        KeyError(
            "no data group `$name`; this Measurement Set has $(join(sort(collect(string.(keys(groups)))), ", "))"
        )
    )
end

"""
    correlated(ms, group = :base) -> AbstractDimArray

The correlated data of one data group — the visibilities or the spectrum.

A Measurement Set may hold several versions of its science array, so which one
is meant depends on the data group. This resolves that indirection.

```julia
correlated(ms)              # the base group's VISIBILITY
correlated(ms, :imaging)    # whatever the imaging group points at
```
"""
correlated(ms::MeasurementSet, group::Symbol = :base) =
    ms[_group_layer(ms, group, :correlated_data)]

"""
    flags(ms, group = :base) -> AbstractDimArray

The flags belonging to one data group. True marks a datum as bad.
"""
flags(ms::MeasurementSet, group::Symbol = :base) = ms[_group_layer(ms, group, :flag)]

"""
    weights(ms, group = :base) -> AbstractDimArray

The weights belonging to one data group, equal to `1/σ²` per datum.
"""
weights(ms::MeasurementSet, group::Symbol = :base) = ms[_group_layer(ms, group, :weight)]

"""
    uvw(ms, group = :base) -> AbstractDimArray

The baseline coordinates belonging to one data group, in metres. Single-dish
data has none.
"""
uvw(ms::MeasurementSet, group::Symbol = :base) = ms[_group_layer(ms, group, :uvw)]

function _group_layer(ms::MeasurementSet, group::Symbol, role::Symbol)
    g = datagroup(ms, group)
    haskey(g, role) || throw(
        KeyError("data group `$group` names no `$role`")
    )
    return julianame(g[role])
end

"""
    antennas(ms) -> Vector{String}

The names of the antennas in `ms`, from its antenna sub-dataset.
"""
function antennas(ms::MeasurementSet)
    haskey(branches(ms), :antenna) || return String[]
    return _lookup_values(branches(ms)[:antenna], :antenna_name)
end

"""
    baselines(ms) -> Vector{Tuple{String, String}}
    baselines(lookup) -> Vector{Tuple{String, String}}

The antenna pair of every baseline, in `baseline_id` order.

MSv4 stores the two names as the `baseline_antenna1_name` and
`baseline_antenna2_name` coordinates over `baseline_id`; this pairs them up.

See also [`baseline`](@ref), [`XRadio.baseline_lookup`](@ref).
"""
baselines(ms::MeasurementSet) = baselines(baseline_lookup(ms))
baselines(bl::BaselineLookup) = collect(zip(bl.antenna1, bl.antenna2))

"""
    baseline(ms, key) -> BaselineID
    baseline(lookup, key) -> BaselineID

A selector for the one baseline that `key` names, for use anywhere a dimension
selector goes:

```julia
vis = correlated(ms)
vis[baseline(ms, ("BR", "LA"))]
uvw(ms)[baseline(ms, "BR-LA")]
```

The key is an antenna-name pair, the dash label `"BR-LA"` spelling that pair,
or an `Integer` `baseline_id` from the store. Pair matching ignores the order
of the two names, since MSv4 records each pair once in whichever order the
correlator gave.

The dash label is an input convenience only. MSv4 defines no such label, and a
separator that no antenna name may contain does not exist, so a name containing
a dash has to be passed as a pair.

Returning a selector rather than data leaves it to the caller which layers are
wanted, and composes with every other selector.

A name key is answered from a [`XRadio.baseline_lookup`](@ref), which a loop
over baselines can hold itself to make the query concrete. A `baseline_id` is
answered by the `baseline_id` dimension, which is already the lookup for it.

See also [`baselines`](@ref).
"""
baseline(ms::MeasurementSet, key) = baseline(baseline_lookup(ms), key)

function baseline(ms::MeasurementSet, id::Integer)
    d = _baseline_dim(ms)
    hasselection(d, At(id)) && return dimtype(:baseline_id)(At(id))
    throw(
        ArgumentError(
            "no baseline_id $id in this Measurement Set; it has $(length(d)) " *
                "baselines, with ids from $(minimum(val(d))) to $(maximum(val(d)))"
        )
    )
end

function baseline(bl::BaselineLookup, key::Tuple{AbstractString, AbstractString})
    a, b = key
    i = get(bl.bypair, (a, b), 0)
    i == 0 && throw(
        ArgumentError(
            "no baseline $a-$b in this Measurement Set; its baselines are over " *
                join(antenna_names(bl), ", ")
        )
    )
    return dimtype(:baseline_id)(At(bl.ids[i]))
end

baseline(bl::BaselineLookup, label::AbstractString) =
    baseline(bl, _split_baseline_label(label))

"""
    baselines(ms, antenna; auto = true) -> BaselineID
    baselines(lookup, antenna; auto = true) -> BaselineID

A selector for every baseline that `antenna` takes part in, in `baseline_id`
order. `auto = false` drops the autocorrelation.

```julia
correlated(ms)[baselines(ms, "LA"; auto = false)]
```

An antenna that no baseline names is an error; an antenna whose only baseline
is its autocorrelation selects nothing under `auto = false`.
"""
baselines(ms::MeasurementSet, antenna::AbstractString; auto::Bool = true) =
    baselines(baseline_lookup(ms), antenna; auto)

function baselines(bl::BaselineLookup, antenna::AbstractString; auto::Bool = true)
    found = get(bl.byantenna, antenna, nothing)
    found === nothing && throw(
        ArgumentError(
            "no baseline of this Measurement Set names antenna $antenna; its " *
                "baselines are over " * join(antenna_names(bl), ", ")
        )
    )
    wanted = auto ? bl.ids[found] :
        [bl.ids[i] for i in found if bl.antenna1[i] != bl.antenna2[i]]
    return dimtype(:baseline_id)(At(wanted))
end

"""
    baseline_antennas(ms, id) -> Tuple{String, String}

The antenna pair of the baseline `id` names, in the order the coordinates give
them — the inverse of [`baseline`](@ref), which goes from a pair to an id.

```julia
baseline_antennas(ms, 7)          # ("BR", "LA")
baseline(ms, ("BR", "LA"))        # BaselineID(At(7))
```

A `baseline_id` the Measurement Set does not have is an error, naming the range
it does have.

There is no method on a [`XRadio.baseline_lookup`](@ref), for the same reason
`baseline` has none taking an id: the `baseline_id` dimension is already the
lookup for an id, and searching the materialized ids would be slower. A loop
over every baseline wants [`baselines`](@ref), which pairs them all at once.

See also [`baselines`](@ref).
"""
function baseline_antennas(ms::MeasurementSet, id::Integer)
    sel = baseline(ms, id)
    antenna1, antenna2 = _baseline_antenna_coords(ms)
    return (antenna1[sel], antenna2[sel])
end

"""
    XRadio.baseline_lookup(ms) -> BaselineLookup

The lookup `ms` answers [`baseline`](@ref) and [`baselines`](@ref) from, built
on the first query and dropped whenever a layer of `ms` changes.

Reading it out of a Measurement Set costs a dynamic dispatch and an allocation
per call, because a Measurement Set holds its layers abstractly. Holding it
directly makes the query concrete, which is what a loop over baselines wants:

```julia
bl = XRadio.baseline_lookup(ms)
map(p -> uvw(ms)[baseline(bl, p)], pairs)
```

It describes `ms` as it stood when it was built. A coordinate mutated in place,
through the array `ms` hands back rather than by assigning a new layer, leaves
it describing baselines the Measurement Set no longer has.
"""
function baseline_lookup(ms::MeasurementSet)
    cached = getfield(ms, :baseline_lookup)
    cached === nothing || return cached
    d = _baseline_dim(ms)
    antenna1, antenna2 = _baseline_antenna_coords(ms)
    built = BaselineLookup(collect(val(d)), collect(antenna1), collect(antenna2))
    setfield!(ms, :baseline_lookup, built)
    return built
end

# The two coordinates MSv4 states the baseline antennas as. Asked for after the
# `baseline_id` dimension, whose absence is the more basic fact and has its own
# message.
function _baseline_antenna_coords(ms::MeasurementSet)
    (haskey(ms, :baseline_antenna1_name) && haskey(ms, :baseline_antenna2_name)) || throw(
        ArgumentError(
            "this Measurement Set names no baseline antennas; MSv4 states them as " *
                "`baseline_antenna1_name` and `baseline_antenna2_name`, and this one " *
                "has " * join(sort(string.(collect(keys(ms)))), ", ")
        )
    )
    return (ms[:baseline_antenna1_name], ms[:baseline_antenna2_name])
end

function _baseline_dim(ms::MeasurementSet)
    d = dims(ms, dimtype(:baseline_id))
    d === nothing && throw(
        ArgumentError(
            "this Measurement Set has no `baseline_id` dimension, so it has no " *
                "baselines; single-dish data is stated over antennas instead"
        )
    )
    return d
end

function _split_baseline_label(label::AbstractString)
    parts = split(label, '-')
    length(parts) == 2 || throw(
        ArgumentError(
            "`$label` is not a baseline label; the form is \"A1-A2\", two antenna " *
                "names separated by a dash. Pass the pair as a tuple if a name " *
                "contains a dash."
        )
    )
    return (String(parts[1]), String(parts[2]))
end

"""
    scans(ms) -> Vector{String}
    scans(ps) -> Vector{String}

The distinct scan names in `ms`, or across every Measurement Set of `ps`, in
the order they first occur.

A Measurement Set stating no `scan_name` gives an empty result, so that a
selection over a processing set whose members differ is not an error.

See also [`scan`](@ref), [`XRadio.scan_lookup`](@ref).
"""
scans(ms::MeasurementSet) = haskey(ms, :scan_name) ? copy(labels(scan_lookup(ms))) : String[]
scans(ps::ProcessingSet) = _distinct(scans, ps)

"""
    fields(ms) -> Vector{String}
    fields(ps) -> Vector{String}

The distinct field names in `ms`, or across every Measurement Set of `ps`, in
the order they first occur.

A Measurement Set stating no `field_name` gives an empty result, so that a
selection over a processing set whose members differ is not an error.

See also [`field`](@ref), [`XRadio.field_lookup`](@ref).
"""
fields(ms::MeasurementSet) = haskey(ms, :field_name) ? copy(labels(field_lookup(ms))) : String[]
fields(ps::ProcessingSet) = _distinct(fields, ps)

"""
    scan(ms, name) -> Ti
    scan(lookup, times, name) -> Ti

A selector for the integrations of the one scan `name` names, for use anywhere
a dimension selector goes:

```julia
correlated(ms)[scan(ms, "3")]
ms[scan(ms, "3")]
uvw(ms)[scan(ms, "3"), baseline(ms, ("BR", "LA"))]
```

Returning a selector rather than data leaves it to the caller which layers are
wanted, and composes with every other selector. The selector names times rather
than positions, so it keeps meaning against a Measurement Set that has since
been cut along another dimension.

A scan the Measurement Set does not name is an error listing the ones it does.

The second form answers from a [`XRadio.scan_lookup`](@ref) and the times a
loop over scans can hold itself, which is what makes the query concrete.

See also [`scans`](@ref), [`field`](@ref).
"""
scan(ms::MeasurementSet, name::AbstractString) =
    _label_selector(scan_lookup(ms), _time_values(ms), name)

"""
    field(ms, name) -> Ti
    field(lookup, times, name) -> Ti

A selector for the integrations `ms` spent on the one field `name` names, for
use anywhere a dimension selector goes:

```julia
correlated(ms)[field(ms, "M87")]
```

The field of a Measurement Set changes only where the observation mosaics or
tracks an ephemeris; a set holding one field gives every integration back.

A field the Measurement Set does not name is an error listing the ones it does.

See also [`fields`](@ref), [`scan`](@ref).
"""
field(ms::MeasurementSet, name::AbstractString) =
    _label_selector(field_lookup(ms), _time_values(ms), name)

scan(ll::LabelLookup, times::AbstractVector, name::AbstractString) =
    _label_selector(ll, times, name)
field(ll::LabelLookup, times::AbstractVector, name::AbstractString) =
    _label_selector(ll, times, name)

_label_selector(ll::LabelLookup, times::AbstractVector, name::AbstractString) =
    dimtype(:time)(At(times[positions(ll, name)]))

"""
    XRadio.scan_lookup(ms) -> LabelLookup
    XRadio.field_lookup(ms) -> LabelLookup

The lookup `ms` answers [`scans`](@ref) and [`scan`](@ref), or [`fields`](@ref)
and [`field`](@ref), from — built on the first query and dropped whenever a
layer of `ms` changes.

MSv4 names a scan and a field per integration, so answering either question
from the coordinate itself is a pass over the time axis. Holding the lookup and
the times keeps a loop over scans off that pass entirely:

```julia
sl = XRadio.scan_lookup(ms)
t = collect(DimensionalData.val(dims(ms, Ti)))
map(s -> correlated(ms)[scan(sl, t, s)], scans(ms))
```

It describes `ms` as it stood when it was built. A coordinate mutated in place,
through the array `ms` hands back rather than by assigning a new layer, leaves
it describing scans the Measurement Set no longer has.

A Measurement Set stating no such coordinate is an error; [`scans`](@ref) and
[`fields`](@ref) are the forms that give an empty result instead.
"""
function scan_lookup(ms::MeasurementSet)
    cached = getfield(ms, :scan_lookup)
    cached === nothing || return cached
    built = LabelLookup(:scan, parent(_time_coord(ms, :scan_name)))
    setfield!(ms, :scan_lookup, built)
    return built
end

"See [`XRadio.scan_lookup`](@ref), whose field twin this is."
function field_lookup(ms::MeasurementSet)
    cached = getfield(ms, :field_lookup)
    cached === nothing || return cached
    built = LabelLookup(:field, parent(_time_coord(ms, :field_name)))
    setfield!(ms, :field_lookup, built)
    return built
end

# A coordinate MSv4 states once per integration. Its positions index the time
# axis, so a coordinate stated over anything else would give a lookup whose
# answers silently name the wrong integrations.
function _time_coord(ms::MeasurementSet, key::Symbol)
    haskey(ms, key) || throw(
        ArgumentError(
            "this Measurement Set states no `$key`; MSv4 names one per integration, " *
                "and this one has " * join(sort(string.(collect(keys(ms)))), ", ")
        )
    )
    A = ms[key]
    d = dims(A)
    # A selection naming one time leaves the coordinate a scalar, which names
    # the one label the remaining integration carries.
    (isempty(d) || (length(d) == 1 && only(d) isa dimtype(:time))) || throw(
        ArgumentError(
            "`$key` is stated over $(join(string.(DimensionalData.name(d)), ", ")), " *
                "not over `time` alone, so it does not name one value per integration"
        )
    )
    return A
end

# The raw values, not the `Lookup` that indexing one gives back: a selector
# carries times, not a second lookup's worth of metadata.
_time_values(ms::MeasurementSet) = parent(val(_time_dim(ms)))

function _time_dim(ms::MeasurementSet)
    d = dims(ms, dimtype(:time))
    d === nothing && throw(
        ArgumentError(
            "this Measurement Set has no `time` dimension, so it has no integrations"
        )
    )
    return d
end

# A processing set partitions one observation, so the same scan or field is
# named by several of its Measurement Sets; `union` keeps each once, in the
# order the set lists them.
_distinct(accessor, ps::ProcessingSet) =
    isempty(ps) ? String[] : mapreduce(accessor, union, ps)

"""
    scan_intents(ms) -> Vector{Intent}

What the scans of `ms` were observed for, as [`Intent`](@ref)s — the store's
`"OBSERVE_TARGET#ON_SOURCE"` parsed into the ASDM enumeration members it names.

MSv4 states these on the `scan_name` coordinate, one list for the Measurement
Set rather than one per scan. The attribute is required, so a store stating
none is not conforming and [`XRadio.check`](@ref) is what says so; this gives
an empty result.

An intent outside the ASDM enumerations throws, naming the value. It cannot
have been written by this package, whose converter takes `Intent`s.

See also [`spectral_window_intents`](@ref), which stays strings because its
vocabulary is not enumerated.
"""
function scan_intents(ms::MeasurementSet)
    haskey(ms, :scan_name) || return Intent[]
    return Intent.(get(metadata(ms[:scan_name]), :scan_intents, String[]))
end

"""
    spectral_window_intents(ms) -> Vector{String}

What the spectral window of `ms` was set up to observe, such as
`"CONTINUUM"` or `"LINE"`.

MSv4 states these on the `frequency` coordinate, beside the spectral window's
name. Intents sit on two different axes and mean different things, which is why
this is separate from [`scan_intents`](@ref) rather than pooled with it.
"""
function spectral_window_intents(ms::MeasurementSet)
    d = dims(ms, dimtype(:frequency))
    d === nothing && return String[]
    return String.(get(metadata(lookup(d)), :spectral_window_intents, String[]))
end

"""
    frequencies(ms) -> Vector

The channel centre frequencies of `ms`, in hertz.
"""
frequencies(ms::MeasurementSet) = _lookup_values(ms, :frequency)

"""
    polarizations(ms) -> Vector{String}

The polarization products of `ms`, such as `["XX", "YY"]`.
"""
polarizations(ms::MeasurementSet) = _lookup_values(ms, :polarization)

"""
    times(ms) -> Vector

The mid-point of each integration in `ms`, in the scale and format its `time`
coordinate declares.
"""
times(ms::MeasurementSet) = _lookup_values(ms, :time)

"""
    gaincurves(ms) -> DimStack
    gaincurves(gain_curve_xds) -> DimStack

The gain curves a Measurement Set records, over the antennas and receptors it
states them for.

The result has two layers — `curve`, holding an
[`XRadio.AbstractGainCurve`](@ref), and `sensitivity`, in kelvin per jansky,
which turns a dimensionless curve into a flux scale. Both formats state the
sensitivity beside the curve rather than within it, so this does too.

```julia
curves = gaincurves(ms)
curves[:curve][AntennaName(At("BR")), ReceptorLabel(At("pol_0"))]
```

There is no frequency axis: MSv4's gain curve dataset has none, and one
Measurement Set covers one spectral window, which [`spectralwindow`](@ref)
names. This is the shape [`fitsidi_gain_curves`](@ref) returns for a whole
file, minus that axis.

Every curve is a [`GainPolynomial`](@ref), since a polynomial in elevation or
in zenith angle is all a conforming store can hold. An antenna and receptor the
store states no curve for has one whose coefficients are `NaN`, which is what
this package says of any quantity it does not have: the layer is one type
throughout rather than a curve beside a marker for its absence, and which
antennas the store does state a curve for is the `sensitivity` layer's own
`NaN`s. A Measurement Set with no gain curve dataset at all gives an empty
result.

Evaluating the whole array is a broadcast, and needs nothing of its own:

```julia
gain.(gaincurves(ms)[:curve], deg2rad(45))
```

keeps the antenna and receptor axes and is `NaN` where the store states no
curve. Across a processing set, `map(gaincurves, ps)` does the same.

The coefficients on disk are against the angle in degrees and come back against
radians, so `gain(curve, elevation_in_radians)` is the gain directly.
"""
function gaincurves(ms::MeasurementSet)
    haskey(branches(ms), :gain_curve) || return _no_gain_curves()
    return gaincurves(branches(ms)[:gain_curve])
end

function gaincurves(xds::Dataset)
    A, R = dimtype(:antenna_name), dimtype(:receptor_label)
    antenna, receptor = dims(xds, A), dims(xds, R)
    (antenna === nothing || receptor === nothing) && throw(ArgumentError(
        "a `gain_curve` dataset is stated over antenna_name and " *
            "receptor_label; this one has " *
            join(string.(DimensionalData.name(dims(xds))), ", ")
    ))
    # MSv4 makes `gain_curve_type` optional, but coefficients that name no
    # angle are not a curve: nothing says whether they are against elevation or
    # against zenith angle, and the two disagree everywhere but the horizon.
    haskey(xds, :gain_curve_type) || throw(ArgumentError(
        "this `gain_curve` dataset carries no `gain_curve_type`, so nothing " *
            "says which angle its coefficients are a polynomial in"
    ))

    coefficients, types = xds[:gain_curve], xds[:gain_curve_type]
    sensitivities = xds[:gain_curve_sensitivity]

    stations, receptors = val(antenna), val(receptor)
    # Every curve a conforming store can hold is a polynomial in one angle, and
    # a combination it states none for is one whose gain is unknown, which is
    # what `NaN` says of a measurement throughout this package. So the layer is
    # one concrete type rather than a curve beside a marker for its absence,
    # and evaluating it broadcasts without dispatching on every element.
    K = Vector{float(eltype(coefficients))}
    unknown = K([NaN])
    curve = Array{GainPolynomial{K}}(undef, length(stations), length(receptors))
    sensitivity = fill(
        oftype(zero(eltype(sensitivities)), NaN),
        length(stations), length(receptors),
    )
    for (i, a) in enumerate(stations), (j, r) in enumerate(receptors)
        s = sensitivities[A(At(a)), R(At(r))]
        k = collect(coefficients[A(At(a)), R(At(r))])
        # The dataset spans every antenna and receptor it names, so one it
        # states no curve for is the padding the writer left: no coefficients
        # and no sensitivity.
        stated = !(all(iszero, k) && isnan(s))
        curve[i, j] = gain_curve_from_msv4(types[A(At(a))], stated ? k : unknown)
        stated && (sensitivity[i, j] = s)
    end

    over = (antenna, receptor)
    return DimStack((;
        curve = DimArray(curve, over),
        sensitivity = DimArray(sensitivity, over; metadata = metadata(sensitivities)),
    ))
end

function _no_gain_curves()
    over = (dimtype(:antenna_name)(String[]), dimtype(:receptor_label)(String[]))
    empty = Array{GainPolynomial{Vector{Float64}}}(undef, 0, 0)
    return DimStack((;
        curve = DimArray(empty, over),
        sensitivity = DimArray(fill(NaN, 0, 0), over),
    ))
end

"""
    spectralwindow(ms) -> String

The name of the spectral window `ms` covers. MSv4 replaces MSv2's numeric data
description index with this name, and one Measurement Set covers exactly one.
"""
function spectralwindow(ms::MeasurementSet)
    d = dims(ms, dimtype(:frequency))
    d === nothing && return ""
    return String(get(metadata(lookup(d)), :spectral_window_name, ""))
end

# What each `select` keyword asks of a Measurement Set. The keywords are
# singular where the accessor is plural, since each names one value to match
# rather than the set being searched.
const _SELECT_BY = (
    antenna = antennas,
    field = fields,
    intent = scan_intents,
    scan = scans,
    spectralwindow = spectralwindow,
)

"""
    select(ps::ProcessingSet; kwargs...) -> ProcessingSet

The Measurement Sets of `ps` that state every attribute asked for, keeping
their names and the set's own attributes.

The keywords are `antenna`, `field`, `intent`, `scan` and `spectralwindow`,
spelled as the accessors are but singular, since each names one value to match
rather than the set being searched. A keyword takes one name or a collection of
them, and a collection matches a Measurement Set stating any of its values.
Several keywords together must all match.

```julia
select(ps; field = "M87")                       # one source
select(ps; field = "M87", spectralwindow = "spw_2")
select(ps; scan = ["3", "4", "5"])              # any of three scans
select(ps; intent = "CALIBRATE_BANDPASS")       # whatever the subintent
```

`intent` takes an [`Intent`](@ref), a `ScanIntents` member or the string form,
and matches against [`scan_intents`](@ref). An intent naming no subintent
matches a scan whatever subintent it states; `"CALIBRATE_BANDPASS#ON_SOURCE"`
matches only that one.

The result is a `ProcessingSet`, so selections chain. `filter(f, ps)` remains
for predicates these keywords do not reach.

An unrecognized keyword is an error rather than a selection that quietly keeps
everything.

This keeps or drops whole Measurement Sets and never cuts one. To cut inside
one, index it with its dimensions:

```julia
ms[baselines(ms, "LA")]                      # baselines touching an antenna
ms[scan(ms, "3")]                            # one scan's integrations
ms[Ti(t0 .. t1)]                             # a time window
ms[baseline_id = At(7)]                      # the store's own dimension name
correlated(ms)[baseline(ms, ("BR", "LA"))]   # one baseline, on a layer
```

A dimension is named either by its type or, as the store names it, by keyword;
the two select identically. A dimension no part of the Measurement Set carries
is an error rather than a selection that quietly keeps everything.

A selector naming a single value drops that dimension, as it does on an array.
A coordinate stated over the dropped dimension alone is left holding one value,
which is the answer to what was selected: after `ms[baseline(ms, ("BR", "LA"))]`
the Measurement Set has no `baseline_id` dimension and `baseline_antenna1_name`
holds `"BR"`. Keep the dimension by selecting a one-element collection instead,
`BaselineID(At([id]))`.
"""
function select(ps::ProcessingSet; kwargs...)
    asked = values(kwargs)
    unknown = setdiff(keys(asked), keys(_SELECT_BY))
    isempty(unknown) || throw(
        ArgumentError(
            "`select` has no keyword $(join(string.('`', unknown, '`'), ", ")); " *
                "it accepts " * join(string.(keys(_SELECT_BY)), ", ")
        )
    )
    return filter(ps) do ms
        all(k -> _states(_SELECT_BY[k](ms), asked[k]), keys(asked))
    end
end

# Whether what a Measurement Set states — one name, as `spectralwindow` gives,
# or many, as `fields` does — includes any of the names asked for.
_states(have, want) = any(w -> _includes(have, w), _wanted(want))

_includes(have::AbstractString, want) = have == want
_includes(have::AbstractVector{Intent}, want) = any(h -> _covers(Intent(want), h), have)
_includes(have, want) = want in have

# A wanted intent that names no subintent asks only what the scan was for, so
# it matches whichever subintent the store states. One that names a subintent
# asks for both.
_covers(want::Intent, have::Intent) =
    want.intent == have.intent &&
    (want.subintent === nothing || want.subintent == have.subintent)

_wanted(want::AbstractString) = (want,)
_wanted(want::Union{Intent, ScanIntents.T}) = (want,)
_wanted(want::Union{AbstractVector, Tuple}) = want
_wanted(want) = throw(
    ArgumentError(
        "a `select` keyword takes a name or a collection of names; got $(repr(want))"
    )
)

# The values of an index coordinate, or an empty vector when the tree has no
# such dimension.
function _lookup_values(tree, name::Symbol)
    d = dims(tree, dimtype(name))
    d === nothing && return []
    return collect(val(d))
end

function Base.show(io::IO, ::MIME"text/plain", ps::ProcessingSet)
    println(io, "ProcessingSet with ", length(ps), length(ps) == 1 ? " Measurement Set" : " Measurement Sets")
    for (i, (name, ms)) in enumerate(pairs(ps))
        if i > 10
            println(io, "  ⋮ ", length(ps) - 10, " more")
            break
        end
        println(io, "  ", name, ": ", _one_line(ms))
    end
    return
end

Base.show(io::IO, ps::ProcessingSet) = print(io, "ProcessingSet(", length(ps), " sets)")

function _one_line(ms::MeasurementSet)
    type = get(metadata(ms), :type, "?")
    shape = join(("$(msv4name(d))=$(length(d))" for d in dims(ms)), ", ")
    return string(type, " (", shape, ")")
end
