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
    scans(ms) -> Vector{String}

The distinct scan names in `ms`, in the order they first occur.
"""
scans(ms::MeasurementSet) = haskey(ms, :scan_name) ? unique(parent(ms[:scan_name])) : String[]

"""
    fields(ms) -> Vector{String}

The distinct field names in `ms`, in the order they first occur.
"""
fields(ms::MeasurementSet) = haskey(ms, :field_name) ? unique(parent(ms[:field_name])) : String[]

"""
    intents(ms) -> Vector{String}

What the observation was for, from `observation_info`. Values look like
`"OBSERVE_TARGET#ON_SOURCE"` or `"CALIBRATE_BANDPASS#ON_SOURCE"`.
"""
function intents(ms::MeasurementSet)
    info = get(metadata(ms), :observation_info, nothing)
    info isa AbstractDict || return String[]
    return String.(get(info, :intents, String[]))
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
store says nothing about holds `nothing` rather than an empty curve, and a
Measurement Set with no gain curve dataset at all gives an empty result.

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
    curve = Array{Union{Nothing, AbstractGainCurve}}(
        nothing, length(stations), length(receptors)
    )
    sensitivity = fill(NaN, length(stations), length(receptors))
    for (i, a) in enumerate(stations), (j, r) in enumerate(receptors)
        s = sensitivities[A(At(a)), R(At(r))]
        k = collect(coefficients[A(At(a)), R(At(r))])
        # The dataset spans every antenna and receptor it names, so one it
        # states no curve for is the padding the writer left: no coefficients
        # and no sensitivity.
        all(iszero, k) && isnan(s) && continue
        curve[i, j] = gain_curve_from_msv4(types[A(At(a))], k)
        sensitivity[i, j] = s
    end

    over = (antenna, receptor)
    return DimStack((;
        curve = DimArray(curve, over),
        sensitivity = DimArray(sensitivity, over; metadata = metadata(sensitivities)),
    ))
end

function _no_gain_curves()
    over = (dimtype(:antenna_name)(String[]), dimtype(:receptor_label)(String[]))
    return DimStack((;
        curve = DimArray(Array{Union{Nothing, AbstractGainCurve}}(nothing, 0, 0), over),
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
