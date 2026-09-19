# FITS-IDI (AIPS Memo 114) to MSv4.
#
# The conversion is a regridding. FITS-IDI stores one row per time-baseline
# pair, with rows simply absent where a baseline did not observe, and packs
# every spectral window into that row's data matrix. MSv4 stores a dense
# rectangular `(time, baseline_id, frequency, polarization)` array per spectral
# window, padded where FITS-IDI had no row. So each partition here reads a
# contiguous range of rows, allocates the grid its distinct times and baselines
# imply, and scatters the rows into it.

# Speed of light, for UU/VV/WW: FITS-IDI writes them in seconds of light
# travel, MSv4 in metres.
#
# The sign is passed through, because both formats put the baseline the same way
# round: the first antenna's position minus the second's. Neither says so
# unambiguously — MSv4 states the direction and the difference as opposites of
# each other — so the test suite measures it from antenna positions on a file of
# each format rather than taking either at its word.
const _C_LIGHT = 299792458.0

# Julian date of the Unix epoch. `DATE` is the whole-day part of a timestamp,
# so subtracting this before scaling to seconds keeps the fractional `TIME`
# from being rounded away.
const _UNIX_EPOCH_JD = 2440587.5

# FITS-IDI tables with no MSv4 home at all. `INTERFEROMETER_MODEL` is listed by
# MSv4 as future schema work; the rest describe calibration machinery MSv4 does
# not model. A file carrying them loses that content here, which the caller has
# to be told rather than left to discover.
const _UNMAPPED_TABLES = [
    "INTERFEROMETER_MODEL", "MODEL_COMPS", "BASELINE", "CALIBRATION", "BANDPASS",
]

# `EQUINOX` names a mean equinox; MSv4 names an astropy frame. The mapping is
# literal, and an unrecognized equinox is an error rather than a default: the
# difference between the candidates is tens of milliarcseconds and would not be
# visible in the output.
#
# AIPS Memo 114r spells the two values `'J2000'` and `'1950.0B'`; `"B1950"` is
# not in the standard but is written in practice. The VLBA writes no `EQUINOX`
# column at all, but a numeric `EPOCH` of 2000.0 or 1950.0, which the memo says
# must be read as an equinox rather than as an epoch.
const _IDI_FRAMES = Dict(
    "J2000" => "fk5", "2000.0" => "fk5",
    "1950.0B" => "fk4", "B1950" => "fk4", "1950.0" => "fk4",
)

# The MJD zero point, which is what the reference implementation writes when a
# release date is unknown.
const _UNKNOWN_RELEASE_DATE = "1858-11-17T00:00:00.000"

# What MSv4 calls the spectral window a FITS-IDI band is. `NO_BAND` counts
# spectral windows — AIPS calls them bands, or IFs — and the file says nothing
# to tell one from another but its position on that axis, so the name carries
# the index in MSv4's own words rather than in the file's.
spectral_window_name(window::Integer) = "spw_$(window - 1)"

# The reference frequency of each spectral window: the file's `REF_FREQ` plus
# the window's own offset from the `FREQUENCY` table. It is what a store
# carries as the frequency axis's `reference_frequency`, and the physical
# quantity that identifies a spectral window wherever one must be labelled.
reference_frequencies(ref_freq, band_freq) = ref_freq .+ band_freq

_iso8601(t::DateTime) = Dates.format(t, "yyyy-mm-ddTHH:MM:SS.sss")

"""
    idi_frame(equinox) -> String

The astropy frame name for a FITS-IDI `EQUINOX`.
"""
idi_frame(equinox::AbstractString) = get(_IDI_FRAMES, uppercase(_trim(equinox))) do
    throw(ArgumentError(
        "FITS-IDI EQUINOX $(repr(equinox)) names no MSv4 frame; known: " *
            join(sort(collect(keys(_IDI_FRAMES))), ", ")
    ))
end

# The VLBA's numeric `EPOCH` column, which the standard says carries equinox
# values. Trailing zeros are dropped so 2000.0 and 2000.00 both land.
idi_frame(equinox::Real) = idi_frame(string(Float64(equinox)))

"""
    source_frames(source_hdu) -> Vector{String}

The astropy frame of each `SOURCE` row, from its `EQUINOX` column or, where the
VLBA writes that instead, its numeric `EPOCH`.
"""
function source_frames(source_hdu)
    data = source_hdu.data
    for name in (:EQUINOX, :EPOCH)
        hascol(data, name) && return map(idi_frame, column(data, name))
    end
    throw(ArgumentError(
        "FITS-IDI SOURCE table has neither an EQUINOX nor an EPOCH column, so " *
            "the frame of its directions is unstated"
    ))
end

"""
    unix_seconds(date, time) -> Float64

POSIX seconds from the two halves of a FITS-IDI timestamp: `date`, the Julian
date of the start of the day, and `time`, the fraction of a day since it.

The format splits them so that the whole-day part stays exactly representable.
Shifting the epoch before scaling to seconds keeps that: the day contributes
about 1.8e9 seconds, whose `Float64` spacing is 240 ns, rather than the 2.5e5
seconds a Julian date would contribute at a spacing of 0.5 ms.
"""
unix_seconds(date, time) = (date - _UNIX_EPOCH_JD) * 86400 + time * 86400

"""
    weight_scaling(weightyp, mode) -> Bool

Whether the on-disk `WEIGHT` column needs the radiometer conversion to become
an MSv4 weight.

`WEIGHTYP` (AIPS Memo 114, Table 14) declares what the column holds.
`"NORMAL"` means it is already a true weight, `1/σ²`. `"CORRELAT"` and
`"CORRTIME"` mean it is a validity fraction in `[0, 1]` — the fraction of the
integration that correlated — which is not a variance at all, and becomes one
only as `f · 2Δν · τ · η²`. An absent `WEIGHTYP` means `"CORRELAT"`, which the
memo states and which matters because that is the case needing conversion.

`mode` is `:auto` to follow the keyword, `:radiometer` to convert regardless,
or `:validity` to write the raw column through unchanged.
"""
function weight_scaling(weightyp::AbstractString, mode::Symbol)
    mode in (:auto, :radiometer, :validity) || throw(ArgumentError(
        "weight_mode must be :auto (follow WEIGHTYP), :radiometer (always " *
            "convert) or :validity (never convert); got $(repr(mode))"
    ))
    mode === :validity && return false
    mode === :radiometer && return true
    return uppercase(_trim(weightyp)) != "NORMAL"
end

"""
    _integration_times(column, keyword, nrow) -> Vector{Float64}

The integration time of every row, in seconds.

`INTTIM` is optional in FITS-IDI, but MSv4's time axis must state an
`integration_time`, so a file omitting it can only be converted with the caller
supplying one. Inferring it from the spacing of the timestamps would present a
guess as a measurement.
"""
function _integration_times(column, keyword, nrow::Int)
    if keyword !== nothing
        column === nothing || @warn "the file states INTTIM; `integration_time` \
                                     overrides it." integration_time = keyword
        return fill(Float64(keyword), nrow)
    end
    column === nothing && throw(ArgumentError(
        "UV_DATA has no INTTIM column, which FITS-IDI permits but MSv4's time " *
            "axis does not: pass `integration_time` to state it in seconds"
    ))
    return column
end

# Scan segmentation. AIPS Memo 114r defines no scan column, so scans are
# normally inferred from source changes and time gaps; a correlator writing one
# anyway states what an inference can only guess, and is believed.
_scan_ranges(scan, source_id, jd, gap) = scan === nothing ?
    scan_ranges(source_id, jd; gap) : scan_ranges(scan, jd; gap = Inf)

# The partition row ranges, given the scan segmentation the file's timestamps
# imply. `()` puts the whole file in one partition per spectral window, which
# only a small file can afford.
function _partition_ranges(partition, scans, source_id, jd, nrow)
    isempty(partition) && return [1:nrow]
    :scan in partition && return scans
    :source in partition && return scan_ranges(source_id, jd; gap = Inf)
    throw(ArgumentError(
        "partition must be (), (:scan,) or (:source,); got $(repr(partition))"
    ))
end

"""
    antenna_bases(antenna_hdu) -> Dict{Int, Tuple{String, String}}

Each antenna's receptor pair, from the `ANTENNA` table's `POLTYA`/`POLTYB`,
keyed by `ANTENNA_NO`.

A multi-band file repeats an antenna once per frequency setup; the receptor
pair does not vary between them, so the first row for each antenna wins.
"""
function antenna_bases(antenna_hdu)
    data = antenna_hdu.data
    no = column(data, :ANTENNA_NO)
    poltya = column(data, :POLTYA)
    poltyb = column(data, :POLTYB)
    bases = Dict{Int, Tuple{String, String}}()
    for r in eachindex(no)
        get!(bases, Int(no[r]), (_trim(poltya[r]), _trim(poltyb[r])))
    end
    return bases
end

"""
    receptor_angles(antenna_hdu, nband) -> Dict{Int, Matrix{Float64}}

Each antenna's two feed position angles in radians, as a `(receptor, band)`
matrix keyed by `ANTENNA_NO`.

`POLAA` and `POLAB` are one angle per band, in degrees; MSv4's
`ANTENNA_RECEPTOR_ANGLE` carries no band axis, so a partition takes its own
band's column.
"""
function receptor_angles(antenna_hdu, nband::Int)
    data = antenna_hdu.data
    no = column(data, :ANTENNA_NO)
    polaa, polab = column(data, :POLAA), column(data, :POLAB)
    angles = Dict{Int, Matrix{Float64}}()
    for r in axes(no, 1)
        haskey(angles, Int(no[r])) && continue
        a, b = rowvalues(polaa, r), rowvalues(polab, r)
        length(a) == nband && length(b) == nband || throw(ArgumentError(
            "FITS-IDI ANTENNA row $r gives $(length(a)) POLAA and $(length(b)) " *
                "POLAB angles, the UV_DATA data matrix $nband bands"
        ))
        angles[Int(no[r])] = deg2rad.(permutedims([a b]))
    end
    return angles
end

"""
    antenna_positions(array_geometry_hdu) -> Dict{Int, Vector{Float64}}

Each antenna's geocentric position in metres, keyed by `NOSTA`.

`STABXYZ` is measured from the array center the `ARRAYX`/`ARRAYY`/`ARRAYZ`
keywords give, which VLBI files leave at the geocenter and connected-element
files do not.
"""
function antenna_positions(ag_hdu)
    frame = card_string(ag_hdu, "FRAME", "GEOCENTRIC")
    frame == "GEOCENTRIC" || throw(ArgumentError(
        "FITS-IDI ARRAY_GEOMETRY FRAME is $(repr(frame)); MSv4 positions are " *
            "geocentric and this converter cannot transform from another frame"
    ))
    center = Float64[card(ag_hdu, k, 0.0) for k in ("ARRAYX", "ARRAYY", "ARRAYZ")]
    data = ag_hdu.data
    nosta = column(data, :NOSTA)
    stabxyz = column(data, :STABXYZ)
    return Dict{Int, Vector{Float64}}(
        Int(nosta[r]) => center .+ Float64.(rowvalues(stabxyz, r)) for r in axes(nosta, 1)
    )
end

"""
    antenna_mounts(array_geometry_hdu) -> Dict{Int, XRadio.AbstractMount}

Each antenna's mount, keyed by `NOSTA`, carrying its `STAXOF` axis offset.
"""
function antenna_mounts(ag_hdu)
    data = ag_hdu.data
    nosta = column(data, :NOSTA)
    mntsta = column(data, :MNTSTA)
    staxof = column(data, :STAXOF)
    mounts = Dict{Int, XRadio.AbstractMount}()
    for r in axes(nosta, 1)
        offset = Float64.(rowvalues(staxof, r))
        length(offset) == 3 || throw(ArgumentError(
            "FITS-IDI ARRAY_GEOMETRY row $r gives $(length(offset)) STAXOF " *
                "values; the standard defines three, one per station axis"
        ))
        mounts[Int(nosta[r])] = idi_mount(mntsta[r], NTuple{3, Float64}(offset))
    end
    return mounts
end

"""
    antenna_diameters(array_geometry_hdu) -> Dict{Int, Float64} or nothing

Each antenna's dish diameter in metres, keyed by `NOSTA`, or `nothing` where
the file omits the optional `DIAMETER` column.
"""
function antenna_diameters(ag_hdu)
    data = ag_hdu.data
    hascol(data, :DIAMETER) || return nothing
    nosta = column(data, :NOSTA)
    diameter = column(data, :DIAMETER)
    return Dict{Int, Float64}(
        Int(nosta[r]) => Float64(diameter[r]) for r in axes(nosta, 1)
    )
end

# The polarization labels for one baseline, and the array-wide labels for a
# whole partition. MSv4's `polarization` coordinate is one list for the entire
# dataset, which presumes every antenna shares a receptor basis; where they do
# not, the majority basis names the axis and the optional `polarization_mixed`
# coordinate carries what each baseline actually measures.
function _partition_polarization(products, bases, codes)
    perbaseline = map(codes) do code
        a1, a2 = baseline_antennas(code)
        for a in (a1, a2)
            haskey(bases, a) || throw(ArgumentError(
                "FITS-IDI BASELINE $code names antenna $a, which the ANTENNA " *
                    "table does not describe"
            ))
        end
        pol_labels(products, bases[a1], bases[a2])
    end
    counts = Dict{Tuple{String, String}, Int}()
    for code in codes, a in baseline_antennas(code)
        counts[bases[a]] = get(counts, bases[a], 0) + 1
    end
    majority = argmax(b -> counts[b], keys(counts))
    labels = pol_labels(products, majority, majority)
    any(!=(labels), perbaseline) || return labels, nothing
    @warn "the array does not share one receptor basis, so no array-wide \
           polarization label describes every baseline. `polarization` takes the \
           majority basis; `polarization_mixed` carries what each baseline \
           measures." polarization = labels maxlog = 1
    return labels, perbaseline
end

# One row of the `FLAG` table. FITS-IDI records flags as ranges over the
# observation rather than as an array beside the data: a row asserts that some
# selection of source, antennas, frequency setup, time, bands, channels and
# polarizations is bad, and a reader expands it onto whatever grid it builds.
# `pols` is in MSv4's polarization order; everything else is as the file states
# it, with the wildcards resolved.
struct FlagRange
    source_id::Int          # 0 covers every source
    ant1::Int               # NOSTA, 0 a wildcard
    ant2::Int
    freqid::Int             # 0 covers every frequency setup
    t0::Float64             # POSIX seconds
    t1::Float64
    bands::Vector{Bool}
    chan_lo::Int            # 1-based, 0 meaning the first channel
    chan_hi::Int            # 1-based, 0 meaning the last
    pols::Vector{Bool}
end

"""
    reference_date(hdu) -> DateTime

The `RDATE` a FITS-IDI table's times are measured from.
"""
function reference_date(hdu)
    text = card_string(hdu, "RDATE")
    isempty(text) && throw(ArgumentError(
        "FITS-IDI $(card_string(hdu, "EXTNAME")) table has no RDATE keyword, " *
            "which its times are measured from"
    ))
    return try
        DateTime(text)
    catch
        throw(ArgumentError(
            "FITS-IDI RDATE $(repr(text)) is not a date; the standard writes it " *
                "as yyyy-mm-dd"
        ))
    end
end

"""
    read_flag_ranges(flag_hdu, perm, layout, rdate) -> Vector{FlagRange}

The `FLAG` table's rows, with their wildcards resolved and their times turned
into POSIX seconds from `rdate`.

`perm[p]` is the on-disk `STOKES` index of MSv4 polarization `p`, so `PFLAGS`
comes back in the order a converted partition's polarization axis uses.
"""
function read_flag_ranges(flag_hdu, perm, layout, rdate::DateTime)
    data = flag_hdu.data
    for name in (:SOURCE_ID, :ANTS, :FREQID, :TIMERANG, :BANDS, :CHANS, :PFLAGS)
        hascol(data, name) || throw(ArgumentError(
            "FITS-IDI FLAG table has no $name column, which the standard requires; " *
                "without it a row does not say what it flags"
        ))
    end
    epoch = datetime2unix(rdate)
    source = column(data, :SOURCE_ID)
    ants, freqid = column(data, :ANTS), column(data, :FREQID)
    timerang, bands = column(data, :TIMERANG), column(data, :BANDS)
    chans, pflags = column(data, :CHANS), column(data, :PFLAGS)

    ranges = Vector{FlagRange}(undef, length(source))
    for r in axes(source, 1)
        a, t = rowvalues(ants, r), rowvalues(timerang, r)
        b, c, p = rowvalues(bands, r), rowvalues(chans, r), rowvalues(pflags, r)
        length(a) == 2 && length(t) == 2 && length(c) == 2 || throw(ArgumentError(
            "FITS-IDI FLAG row $r gives $(length(a)) ANTS, $(length(t)) TIMERANG " *
                "and $(length(c)) CHANS values; the standard defines two of each"
        ))
        length(b) == layout.no_band || throw(ArgumentError(
            "FITS-IDI FLAG row $r gives $(length(b)) BANDS values, the UV_DATA " *
                "data matrix $(layout.no_band) bands"
        ))
        length(p) == layout.no_stkd || throw(ArgumentError(
            "FITS-IDI FLAG row $r gives $(length(p)) PFLAGS values, the UV_DATA " *
                "data matrix $(layout.no_stkd) polarizations"
        ))
        # A literal `(0, 0)` would be an empty interval at the reference date;
        # AIPS writes it for a flag that covers the whole observation.
        t0, t1 = iszero(t[1]) && iszero(t[2]) ? (-Inf, Inf) :
            (epoch + t[1] * 86400, epoch + t[2] * 86400)
        ranges[r] = FlagRange(
            Int(source[r]), Int(a[1]), Int(a[2]), Int(freqid[r]), t0, t1,
            [!iszero(x) for x in b], Int(c[1]), Int(c[2]),
            [!iszero(p[perm[i]]) for i in eachindex(perm)],
        )
    end
    return ranges
end

# Whether a `FLAG` row's antenna selection covers a baseline. One antenna flags
# every baseline touching it, two flag that baseline whichever way round the
# file writes them, and neither flags them all.
function _flags_baseline(e::FlagRange, a::Integer, b::Integer)
    e.ant1 == 0 && e.ant2 == 0 && return true
    e.ant2 == 0 && return e.ant1 == a || e.ant1 == b
    e.ant1 == 0 && return e.ant2 == a || e.ant2 == b
    return (e.ant1 == a && e.ant2 == b) || (e.ant1 == b && e.ant2 == a)
end

# Expand onto one partition's grid the flag ranges that touch it. A row names a
# selection rather than a cell, so it is resolved to index vectors once and
# written as a slice.
function _apply_flags!(flag, ranges, codes, times, source_at, band, freqid, nchan)
    isempty(ranges) && return flag
    antennas = map(baseline_antennas, codes)
    # A correlator writes one row per antenna dropout, so a real table has
    # thousands of them and all but a few cover times outside any one scan.
    # `times` is sorted, so its ends bound the partition.
    span = extrema(times)
    for e in ranges
        (e.t1 >= span[1] && e.t0 <= span[2]) || continue
        (e.freqid == 0 || e.freqid == freqid) || continue
        e.bands[band] || continue
        ks = [
            k for k in eachindex(times) if e.t0 <= times[k] <= e.t1 &&
                (e.source_id == 0 || e.source_id == source_at[k])
        ]
        isempty(ks) && continue
        bs = [b for b in eachindex(antennas) if _flags_baseline(e, antennas[b]...)]
        isempty(bs) && continue
        ps = findall(e.pols)
        isempty(ps) && continue
        cs = max(e.chan_lo, 1):(e.chan_hi == 0 ? nchan : min(e.chan_hi, nchan))
        isempty(cs) && continue
        flag[ps, cs, bs, ks] .= true
    end
    return flag
end

# ── Calibration and environment tables ───────────────────────────────────────
#
# FITS-IDI samples these on the instrument's own clock rather than on the
# correlated data's time axis, and MSv4 gives each of them a time dimension of
# its own for that reason. Every one of these tables is at most a few thousand
# rows, so each is read once in full and every partition selects from it the
# rows whose interval overlaps its own time span and whose antenna it
# describes.

# One time-sampled calibration table: the columns they all share, already
# restricted to the file's frequency setup, and `values` holding the rest.
struct CalTable{V}
    time::Vector{Float64}       # POSIX seconds, centre of the interval
    interval::Vector{Float64}   # seconds
    antenna::Vector{Int}        # NOSTA; 0 stands for every antenna
    values::V
end

# The first of `names` the table carries. Revision 1 of the WEATHER and
# PHASE-CAL tables spells several columns differently from revision 2 and the
# standard names both spellings, so a reader that knows only one rejects files
# that are in the standard.
function _cal_column(hdu, names::Symbol...)
    data = hdu.data
    for name in names
        hascol(data, name) && return collect(column(data, name))
    end
    throw(ArgumentError(
        "FITS-IDI $(card_string(hdu, "EXTNAME")) table has no " *
            join(names, " or ") * " column, which the standard requires"
    ))
end

# The first of `names` the table carries as one value a row, or `nothing` where
# it carries none of them. A column of `leng` values a row comes back as an
# `(nrow, leng)` matrix, and `leng` may be zero: the VLBA declares `WVR_H2O`
# and `IONOS_ELECTRON` that way in tracks that measured neither, so a column
# being present is not the same as its holding anything.
function _cal_optional(hdu, names::Symbol...)
    data = hdu.data
    for name in names
        hascol(data, name) || continue
        col = collect(column(data, name))
        col isa AbstractVector && return col
    end
    return nothing
end

# The same, where the standard requires the column.
function _cal_scalar(hdu, names::Symbol...)
    col = _cal_optional(hdu, names...)
    col === nothing && throw(ArgumentError(
        "FITS-IDI $(card_string(hdu, "EXTNAME")) table has no " *
            join(names, " or ") * " column of one value a row, which the " *
            "standard requires"
    ))
    return col
end

# A column FITS-IDI writes once per receptor, as `NAME_1` and `NAME_2`, read as
# a `(value, receptor, row)` array. A table of one polarization leaves the
# second receptor unmeasured, which is `NaN` rather than a repeat of the first.
# Integer-coded columns come back as floats so that the same absence is
# expressible in them.
function _receptor_column(hdu, name::AbstractString, npol::Int, len::Int, keep)
    cols = [_cal_column(hdu, Symbol(name, "_", p), Symbol(name, " ", p)) for p in 1:npol]
    # The file's own precision, not `Float64`. An MSv4 array these feed is
    # `Float64` and the dataset builders convert to it; a gain curve reaches a
    # caller through `fitsidi_gain_curves` without passing through a store, and
    # FITS-IDI writes these columns in single precision.
    T = float(mapreduce(eltype, promote_type, cols))
    out = fill(convert(T, NaN), len, 2, length(keep))
    for (p, col) in pairs(cols)
        for (i, r) in enumerate(keep)
            v = rowvalues(col, r)
            length(v) == len || throw(ArgumentError(
                "FITS-IDI $(card_string(hdu, "EXTNAME")) row $r gives $(length(v)) " *
                    "$(name)_$p values, expected $len"
            ))
            out[:, p, i] = v
        end
    end
    return out
end

# How many receptors the table describes. `NO_POL` states it; a file that
# omits the keyword still has the second receptor's columns exactly when it
# describes two.
function _no_pol(hdu, probe::AbstractString)
    n = Int(card(hdu, "NO_POL", hascol(hdu.data, Symbol(probe, "_2")) ? 2 : 1))
    n in (1, 2) || throw(ArgumentError(
        "FITS-IDI $(card_string(hdu, "EXTNAME")) table has NO_POL = $n; a feed " *
            "has one or two receptors"
    ))
    return n
end

# The number of entries a per-band column packs for each band. FITS-IDI writes
# these columns as `(n, nband)` matrices with `n` varying fastest, and states
# `n` in a keyword the memo names two ways; the column's own width states it
# unambiguously.
function _per_band_width(hdu, name::Symbol, nband::Int)
    col = _cal_column(hdu, name)
    width, rem = divrem(length(rowvalues(col, 1)), nband)
    (rem == 0 && width >= 1) || throw(ArgumentError(
        "FITS-IDI $(card_string(hdu, "EXTNAME")) column $name holds " *
            "$(length(rowvalues(col, 1))) values per row, which is not a whole " *
            "number of entries for each of $nband bands"
    ))
    return width
end

# The rows of a calibration table another frequency setup does not claim, and
# their times as POSIX seconds. FITS-IDI measures them in days from the table's
# own `RDATE`.
function _cal_rows(hdu, freqid::Int)
    time = Float64.(_cal_scalar(hdu, :TIME))
    interval = Float64.(_cal_scalar(hdu, :TIME_INTERVAL))
    antenna = Int.(_cal_scalar(hdu, :ANTENNA_NO, :ANTENNA_NUMBER))
    keep = if hascol(hdu.data, :FREQID)
        f = Int.(_cal_scalar(hdu, :FREQID))
        findall(x -> x == 0 || x == freqid, f)
    else
        collect(eachindex(time))
    end
    epoch = datetime2unix(reference_date(hdu))
    return keep, epoch .+ time[keep] .* 86400, interval[keep] .* 86400, antenna[keep]
end

# Where a table's rows land on the `(antenna, time)` grid one partition wants,
# as `(row, antenna index, time index)` triples. A row is taken when its
# interval overlaps the partition's time span, so a scan is described by the
# measurements that cover it rather than by the whole track's. `0` in the
# antenna column stands for every antenna, which WEATHER writes for conditions
# measured once for the site.
function _cal_grid(tab::CalTable, numbers, span)
    index = Dict(a => i for (i, a) in enumerate(numbers))
    keep = [
        r for r in eachindex(tab.time)
            if tab.time[r] + tab.interval[r] / 2 >= span[1] &&
                tab.time[r] - tab.interval[r] / 2 <= span[2] &&
                (tab.antenna[r] == 0 || haskey(index, tab.antenna[r]))
    ]
    isempty(keep) && return nothing
    times = sort!(unique(tab.time[keep]))
    at = Dict(t => k for (k, t) in enumerate(times))
    places = NTuple{3, Int}[]
    for r in keep
        k = at[tab.time[r]]
        if tab.antenna[r] == 0
            append!(places, ((r, i, k) for i in eachindex(numbers)))
        else
            push!(places, (r, index[tab.antenna[r]], k))
        end
    end
    return (; times, places)
end

"""
    read_system_calibration(hdus, layout, freqid) -> CalTable or nothing

The `SYSTEM_TEMPERATURE` table's system and antenna temperatures, as
`(receptor, band, row)` arrays, or `nothing` where the file omits the table.
"""
function read_system_calibration(hdus, layout, freqid::Int)
    hdu = table(hdus, "SYSTEM_TEMPERATURE")
    hdu === nothing && return nothing
    keep, time, interval, antenna = _cal_rows(hdu, freqid)
    npol, nband = _no_pol(hdu, "TSYS"), layout.no_band
    # MSv4 states these as `Float64`, so the conversion happens here, where the
    # arrays are being shaped for a store.
    receptor(name) = Float64.(
        permutedims(_receptor_column(hdu, name, npol, nband, keep), (2, 1, 3))
    )
    return CalTable(
        time, interval, antenna, (; tsys = receptor("TSYS"), tant = receptor("TANT"))
    )
end

"""
    read_phase_calibration(hdus, layout, freqid) -> CalTable or nothing

The `PHASE-CAL` table's injected tones and cable delay, or `nothing` where the
file omits the table.

`cable` is one delay per row; `tone`, `real` and `imag` are
`(tone, receptor, band, row)`. `PC_RATE` and `STATE` have no MSv4 field and are
not read.
"""
function read_phase_calibration(hdus, layout, freqid::Int)
    hdu = table(hdus, "PHASE-CAL")
    hdu === nothing && return nothing
    keep, time, interval, antenna = _cal_rows(hdu, freqid)
    npol, nband = _no_pol(hdu, "PC_FREQ"), layout.no_band
    ntone = _per_band_width(hdu, :PC_FREQ_1, nband)
    # MSv4 states these as `Float64`, so the conversion happens here, where the
    # arrays are being shaped for a store.
    receptor(name) = Float64.(
        permutedims(
            reshape(
                _receptor_column(hdu, name, npol, ntone * nband, keep),
                ntone, nband, 2, length(keep)
            ),
            (1, 3, 2, 4)
        )
    )
    cable = Float64.(_cal_scalar(hdu, :CABLE_CAL))[keep]
    return CalTable(
        time, interval, antenna,
        (;
            cable, tone = receptor("PC_FREQ"),
            real = receptor("PC_REAL"), imag = receptor("PC_IMAG"),
        )
    )
end

"""
    read_weather(hdus, freqid) -> CalTable or nothing

The `WEATHER` table, in MSv4's units — kelvin, pascals and radians rather than
FITS-IDI's centigrade, millibars and degrees — or `nothing` where the file
omits the table.

Each field is `nothing` where the file omits its column; MSv4 makes them all
optional too. `PRECIPITATION` has no MSv4 field and is not read.
"""
function read_weather(hdus, freqid::Int)
    hdu = table(hdus, "WEATHER")
    hdu === nothing && return nothing
    keep, time, interval, antenna = _cal_rows(hdu, freqid)
    take(f, names...) = let col = _cal_optional(hdu, names...)
        col === nothing ? nothing : f.(Float64.(col[keep]))
    end
    return CalTable(
        time, interval, antenna,
        (;
            temperature = take(c -> c + 273.15, :TEMPERATURE),
            dew_point = take(c -> c + 273.15, :DEWPOINT),
            pressure = take(mbar -> mbar * 100, :PRESSURE),
            wind_speed = take(identity, :WIND_VELOCITY),
            wind_direction = take(deg2rad, :WIND_DIRECTION),
            h2o = take(identity, :WVR_H2O, :H2O_COLUMN),
            ionos_electron = take(identity, :IONOS_ELECTRON, :ELECTRON_COLUMN),
        )
    )
end

# The angle a FITS-IDI gain curve is a function of, from its `Y_TYP` code.
# AIPS Memo 114r's table of value types also admits 0, which is no angle at
# all and appears where a curve does not use the slot.
const _IDI_GAIN_VARIABLES = Dict(
    1 => :elevation, 2 => :zenith_angle, 3 => :hour_angle,
    4 => :declination, 5 => :codeclination,
)

# A curve's own variable, refusing a code that names none. FITS-IDI writes
# these angles in degrees and XRadio holds them in radians.
function _idi_gain_variable(code, what::AbstractString)
    v = isfinite(code) ? get(_IDI_GAIN_VARIABLES, Int(code), nothing) : nothing
    v === nothing && throw(ArgumentError(
        "FITS-IDI GAIN_CURVE gives $what = $code, which names no angle a gain " *
            "curve can depend on; the standard defines 1 to 5"
    ))
    return v
end

"""
    read_gain_curves(hdus, layout, freqid) -> Vector{<:NamedTuple}

Every gain curve the `GAIN_CURVE` table holds, one row per antenna, receptor
and band: `(; antenna, receptor, band, curve, sensitivity)`, where `antenna` is
a `NOSTA` number and `curve` an [`XRadio.AbstractGainCurve`](@ref).

All three shapes AIPS Memo 114r defines are read — the polynomial, the table of
measured values, and the spherical-harmonic expansion in hour angle and
co-declination — whether or not MSv4 can store them. `sensitivity` is the
curve's `SENS` in kelvin per jansky, which is what turns a dimensionless curve
into a flux scale and which both formats state beside the curve rather than
within it.
"""
# The precision a value read from the table is carried in. Keeping the file's
# own is the default, since nothing `fitsidi_gain_curves` returns passes through
# a store.
_at_precision(::Nothing, x) = x
_at_precision(T::Type, x) = T.(x)

function read_gain_curves(hdus, layout, freqid::Int, T = nothing)
    hdu = table(hdus, "GAIN_CURVE")
    hdu === nothing && return NamedTuple[]
    data = hdu.data
    nband = layout.no_band
    antenna = Int.(_cal_scalar(hdu, :ANTENNA_NO))
    keep = if hascol(data, :FREQID)
        f = Int.(_cal_scalar(hdu, :FREQID))
        findall(r -> f[r] == 0 || f[r] == freqid, eachindex(f))
    else
        collect(eachindex(antenna))
    end
    npol = _no_pol(hdu, "TYPE")
    ntab = _per_band_width(hdu, :GAIN_1, nband)

    type = _receptor_column(hdu, "TYPE", npol, nband, keep)
    nterm = _receptor_column(hdu, "NTERM", npol, nband, keep)
    xtyp = _receptor_column(hdu, "X_TYP", npol, nband, keep)
    ytyp = _receptor_column(hdu, "Y_TYP", npol, nband, keep)
    xval = _receptor_column(hdu, "X_VAL", npol, nband, keep)
    sens = _receptor_column(hdu, "SENS", npol, nband, keep)
    shape(name) = reshape(
        _receptor_column(hdu, name, npol, ntab * nband, keep),
        ntab, nband, 2, length(keep)
    )
    yval, gain = shape("Y_VAL"), shape("GAIN")

    rows = NamedTuple[]
    for (i, r) in enumerate(keep), p in 1:npol, band in 1:nband
        isnan(type[band, p, i]) && continue
        n = nterm[band, p, i]
        (isfinite(n) && 1 <= n <= ntab) || throw(ArgumentError(
            "FITS-IDI GAIN_CURVE row $r gives NTERM_$p = $n for band $band; a " *
                "curve has between one and the table's $ntab entries"
        ))
        n = Int(n)
        coefficients = _at_precision(T, gain[1:n, band, p, i])
        # Codes 1, 2 and 3 are a table of measured values, a polynomial, and a
        # spherical-harmonic expansion.
        curve = if type[band, p, i] == 1
            TabulatedGain(
                _idi_gain_variable(ytyp[band, p, i], "Y_TYP_$p"), :power,
                deg2rad.(_at_precision(T, yval[1:n, band, p, i])), coefficients
            )
        elseif type[band, p, i] == 2
            GainPolynomial(
                _idi_gain_variable(ytyp[band, p, i], "Y_TYP_$p"), :power,
                XRadio.to_radian_angle(coefficients)
            )
        elseif type[band, p, i] == 3
            # Memo 114r lists ten coefficients and lets a file state fewer; a
            # term the expansion stops before contributes nothing. More than
            # ten is not an expansion of this shape, which the constructor
            # says.
            SphericalHarmonicGain(
                :power, append!(coefficients, zeros(max(10 - n, 0)))
            )
        else
            throw(ArgumentError(
                "FITS-IDI GAIN_CURVE row $r gives TYPE_$p = $(type[band, p, i]) " *
                    "for band $band; the standard defines 1 (tabulated), " *
                    "2 (polynomial) and 3 (spherical harmonic)"
            ))
        end
        # A curve against hour angle holds at one declination, which the
        # standard puts in `X_VAL` and names in `X_TYP`.
        if isfinite(xtyp[band, p, i]) && Int(xtyp[band, p, i]) == 4
            curve = AtDeclination(curve, deg2rad(_at_precision(T, xval[band, p, i])))
        end
        push!(
            rows,
            (;
                antenna = antenna[r], receptor = p, band,
                curve, sensitivity = _at_precision(T, sens[band, p, i]),
            )
        )
    end
    return rows
end

# The gain curves a conforming store can hold, keyed the way a partition asks
# for them, or `nothing` where the file has none it can. A curve MSv4 cannot
# express costs the whole dataset rather than part of it: an antenna dataset
# describing some of the array's gain curves and silently not others would be
# read as though the rest had none.
function msv4_gain_curves(rows, path::AbstractString)
    isempty(rows) && return nothing
    curves = Dict{Tuple{Int, Int, Int}, NamedTuple}()
    for row in rows
        written = XRadio.msv4_gain_curve(row.curve)
        if written === nothing
            @warn "MSv4 models a gain curve as a polynomial in elevation or in \
                   zenith angle, and this file holds one of another shape. No \
                   `gain_curve_xds` is written; `fitsidi_gain_curves` reads \
                   them all." file = path curve = typeof(row.curve)
            return nothing
        end
        curves[(row.antenna, row.receptor, row.band)] =
            (; written.type, written.coefficients, row.sensitivity)
    end
    return curves
end

# A measurement FITS-IDI writes as `NaN` is one it does not have, so an array
# of nothing else says nothing and is left out rather than written.
_all_nan(x) = all(isnan, x)

# The identity coordinates every antenna-keyed sub-dataset repeats: what the
# antenna is called, where it stands, and what its receptors carry. MSv4's
# system calibration dataset is the one that leaves out the mount and the
# telescope, which is why they are separable here.
function _antenna_identity(s, numbers; mount::Bool = true)
    names = [s.antenna_name[a] for a in numbers]
    identity = (;
        antenna_name = names,
        station_name = copy(names),
        receptor_label = ["pol_0", "pol_1"],
        polarization_type = [s.bases[a][i] for i in 1:2, a in numbers],
    )
    mount || return identity
    return (;
        identity...,
        mount = [XRadio.msv4_mount(s.mounts[a]) for a in numbers],
        telescope_name = fill(s.telescope_name, length(numbers)),
    )
end

# A calibration dataset's own time axis is the instrument's clock, not the
# correlated data's. The standard fixes what a `time` coordinate carries and
# leaves `time_phase_cal` open, so that one is stated here.
const _CAL_CLOCK = (; units = "s", scale = "utc", format = "unix", type = "time")

# The system and antenna temperatures covering one partition. FITS-IDI states
# one temperature per band, and a partition is one band.
function _system_calibration_dataset(s, numbers, span, band::Int)
    s.system_calibration === nothing && return nothing
    grid = _cal_grid(s.system_calibration, numbers, span)
    grid === nothing && return nothing
    vals = s.system_calibration.values

    nant, ntime = length(numbers), length(grid.times)
    tsys, tant = fill(NaN, 2, ntime, nant), fill(NaN, 2, ntime, nant)
    for (r, a, k) in grid.places
        tsys[:, k, a] = vals.tsys[:, band, r]
        tant[:, k, a] = vals.tant[:, band, r]
    end
    # The memo writes an unmeasured temperature as `NaN`, so a partition whose
    # every entry is one has nothing to say and the array is left out.
    _all_nan(tsys) && _all_nan(tant) && return nothing

    return XRadio.subdataset(
        :system_calibration;
        _antenna_identity(s, numbers; mount = false)...,
        time_system_cal = grid.times,
        tsys = _all_nan(tsys) ? nothing : tsys,
        tant = _all_nan(tant) ? nothing : tant,
    )
end

# The injected phase-calibration tones covering one partition.
function _phase_calibration_dataset(s, numbers, span, band::Int)
    s.phase_calibration === nothing && return nothing
    grid = _cal_grid(s.phase_calibration, numbers, span)
    grid === nothing && return nothing
    vals = s.phase_calibration.values

    nant, ntime = length(numbers), length(grid.times)
    ntone = size(vals.tone, 1)
    tone_freq = fill(NaN, ntone, 2, ntime, nant)
    tones = fill(ComplexF64(NaN, NaN), ntone, 2, ntime, nant)
    cable, interval = fill(NaN, ntime, nant), fill(NaN, ntime, nant)
    for (r, a, k) in grid.places
        tone_freq[:, :, k, a] = vals.tone[:, :, band, r]
        tones[:, :, k, a] =
            complex.(vals.real[:, :, band, r], vals.imag[:, :, band, r])
        cable[k, a] = vals.cable[r]
        interval[k, a] = s.phase_calibration.interval[r]
    end

    return XRadio.subdataset(
        :phase_calibration;
        _antenna_identity(s, numbers)...,
        time_phase_cal = grid.times,
        tone_label = ["tone_$(i - 1)" for i in 1:ntone],
        phase_cal = tones,
        phase_cal_tone_frequency = tone_freq,
        phase_cal_cable_cal = cable,
        phase_cal_interval = interval,
        array_metadata = (; time_phase_cal = _CAL_CLOCK),
    )
end

# The weather covering one partition. MSv4 keys it by station rather than by
# antenna, and a FITS-IDI row whose `ANTENNA_NO` is zero describes the whole
# site, so it is spread across every station the partition names.
function _weather_dataset(s, numbers, span)
    s.weather === nothing && return nothing
    grid = _cal_grid(s.weather, numbers, span)
    grid === nothing && return nothing
    vals = s.weather.values

    nant, ntime = length(numbers), length(grid.times)
    measured = Dict{Symbol, Any}()
    for key in keys(vals)
        col = vals[key]
        col === nothing && continue
        field = fill(NaN, ntime, nant)
        for (r, a, k) in grid.places
            field[k, a] = col[r]
        end
        _all_nan(field) || (measured[key] = field)
    end
    # Only the station positions would be left when the file measured nothing,
    # and positions are already in the antenna dataset.
    isempty(measured) && return nothing

    return XRadio.subdataset(
        :weather;
        station_name = [s.antenna_name[a] for a in numbers],
        cartesian_pos_label = ["x", "y", "z"],
        time_weather = grid.times,
        station_position = [s.positions[a][i] for i in 1:3, a in numbers],
        pairs(measured)...,
    )
end

# The gain curves of one partition's antennas for its spectral window. The
# file states a curve for each, and MSv4's gain curve dataset has no frequency
# axis at all, so a partition takes its own window's — as it does for the
# receptor angles.
function _gain_curve_dataset(s, numbers, band::Int)
    s.gain_curve === nothing && return nothing
    curves = s.gain_curve
    present = [(a, p) for a in numbers for p in 1:2 if haskey(curves, (a, p, band))]
    isempty(present) && return nothing

    types = unique(curves[(a, p, band)].type for (a, p) in present)
    length(types) == 1 || throw(ArgumentError(
        "FITS-IDI GAIN_CURVE gives band $band curves of more than one kind " *
            "($(join(types, ", "))); MSv4 states one `gain_curve_type` for " *
            "each antenna"
    ))
    nterm = maximum(length(curves[(a, p, band)].coefficients) for (a, p) in present)

    nant = length(numbers)
    # A polynomial that stops early has no further terms, which is a zero
    # coefficient rather than an absent one. These are against the angle in
    # degrees, which is what casacore and AIPS write and MSv4 inherits.
    coefficients = zeros(2, nterm, nant)
    sensitivity = fill(NaN, 2, nant)
    for (i, a) in enumerate(numbers), p in 1:2
        curve = get(curves, (a, p, band), nothing)
        curve === nothing && continue
        coefficients[p, eachindex(curve.coefficients), i] = curve.coefficients
        sensitivity[p, i] = curve.sensitivity
    end

    return XRadio.subdataset(
        :gain_curve;
        _antenna_identity(s, numbers)...,
        poly_term = collect(Int64, 0:(nterm - 1)),
        gain_curve_type = fill(only(types), nant),
        gain_curve = coefficients,
        gain_curve_sensitivity = sensitivity,
        # FITS-IDI states one gain curve for a file, so it holds for as long as
        # the file's observation lasts.
        gain_curve_interval = fill(s.observation_span, nant),
        metadata = (; measured_date = s.gain_curve_date),
    )
end

# The antennas of one partition, in antenna-number order. A scan's baselines
# name the antennas that actually observed together, which is fewer than the
# array whenever a station was down.
function _antenna_dataset(s, numbers, band::Int)
    identity = _antenna_identity(s, numbers)

    unnamed = [
        n for (n, a) in zip(identity.antenna_name, numbers)
            if s.mounts[a] isa XRadio.MountOrbiting || s.mounts[a] isa XRadio.MountOther
    ]
    isempty(unnamed) || @warn "MSv4 has no `mount` name for these antennas' \
                               mounts; writing \"OTHER\", which does not say \
                               which they are." antennas = unnamed maxlog = 1

    return XRadio.subdataset(
        :antenna;
        identity...,
        cartesian_pos_label = ["x", "y", "z"],
        antenna_position = [s.positions[a][i] for i in 1:3, a in numbers],
        antenna_receptor_angle =
            [s.receptor_angles[a][i, band] for i in 1:2, a in numbers],
        antenna_dish_diameter = s.diameters === nothing ? nothing :
            [s.diameters[a] for a in numbers],
    )
end

# The fields of one partition. FITS-IDI states one direction per source, which
# is both where the source is and what the visibilities are phased to: the
# format has no notion of pointing somewhere other than at the source.
function _field_and_source_dataset(s, sources, frame::AbstractString)
    names = [s.source_names[i] for i in sources]
    direction = permutedims(
        [deg2rad.(s.source_ra[sources]) deg2rad.(s.source_dec[sources])]
    )
    return XRadio.subdataset(
        :field_and_source;
        field_name = names,
        sky_dir_label = ["ra", "dec"],
        source_name = copy(names),
        field_phase_center_direction = direction,
        source_direction = copy(direction),
        array_metadata = (;
            field_phase_center_direction = (; frame),
            source_direction = (; frame),
        ),
    )
end

"""
    fitsidi_gain_curves(path) -> Vector{<:NamedTuple}

See [`XRadio.fitsidi_gain_curves`](@ref).
"""
function XRadio.fitsidi_gain_curves(path::AbstractString)
    hdus = fits(path)
    uv = require_table(hdus, "UV_DATA")
    fq = require_table(hdus, "FREQUENCY")
    freqid = Int(column(fq.data, :FREQID)[1])
    agd = require_table(hdus, "ARRAY_GEOMETRY").data
    named = Dict{Int, String}(
        Int(n) => _trim(a)
            for (n, a) in zip(column(agd, :NOSTA), column(agd, :ANNAME))
    )
    layout = uv_matrix_layout(uv)
    rows = read_gain_curves(hdus, layout, freqid)

    numbers = sort!(unique(row.antenna for row in rows))
    antenna = AntennaName(Categorical([get(named, a, string(a)) for a in numbers]))
    receptor = ReceptorLabel(Categorical(["pol_0", "pol_1"]))

    # MSv4 has no dimension for a spectral window — a Measurement Set is one —
    # and states which window a frequency axis describes in that axis's own
    # attributes. So the axis here is `frequency`, one entry a window, carrying
    # the reference frequency the store carries for the same window.
    freqs = reference_frequencies(
        Float64(card(uv, "REF_FREQ", 0.0)),
        Float64.(rowvalues(column(fq.data, :BANDFREQ), 1)),
    )
    allunique(freqs) || throw(ArgumentError(
        "FITS-IDI FREQUENCY gives $(length(freqs)) spectral windows sharing " *
            "reference frequencies ($(join(freqs, ", ")) Hz), so no frequency " *
            "tells them apart"
    ))
    frequency = Frequency(
        Sampled(
            freqs;
            metadata = Dict{Symbol, Any}(
                :type => "spectral_coord", :units => "Hz", :observer => "TOPO"
            )
        )
    )
    at = Dict(a => i for (i, a) in enumerate(numbers))

    # A file need not describe every antenna, receptor and spectral window, and
    # a combination it does not describe has no curve rather than an empty one.
    curve = Array{Union{Nothing, XRadio.AbstractGainCurve}}(
        nothing, length(numbers), 2, layout.no_band
    )
    # The file's own precision, as the curves themselves keep: nothing here
    # passes through a store, whose `Float64` this would otherwise impose.
    S = isempty(rows) ? Float64 : float(mapreduce(r -> typeof(r.sensitivity), promote_type, rows))
    sensitivity = fill(convert(S, NaN), length(numbers), 2, layout.no_band)
    for row in rows
        curve[at[row.antenna], row.receptor, row.band] = row.curve
        sensitivity[at[row.antenna], row.receptor, row.band] = row.sensitivity
    end

    dims = (antenna, receptor, frequency)
    return DimStack((;
        curve = DimArray(curve, dims),
        sensitivity = DimArray(
            sensitivity, dims;
            metadata = Dict{Symbol, Any}(:type => "quantity", :units => "K/Jy")
        ),
    ))
end

# MSv4 constrains `units` to a single literal for every array but `VISIBILITY` and
# `SPECTRUM`, where it is a free string defaulting to `Jy`. Correlator output is
# rarely in janskys, so the value comes from the `FLUX` column's `TUNIT`, spelled
# as the file spells it. `JY` becomes `Jy` so a calibrated store agrees with the
# rest of MSv4 on the name of a jansky.
function _visibility_units(uv, path)
    u = column_units(uv, "FLUX")
    if u === nothing
        @warn "FITS-IDI FLUX column carries no TUNIT; writing `uncalib` as \
               `VISIBILITY.units`. Data calibrated elsewhere needs its units \
               stated in the file." file = path
        return "uncalib"
    end
    return uppercase(u) == "JY" ? "Jy" : u
end

"""
    fitsidi2msv4(path, dest; kwargs...) -> dest

See [`XRadio.fitsidi2msv4`](@ref).
"""
function XRadio.fitsidi2msv4(
        path::AbstractString, dest::AbstractString;
        partition = (:scan,),
        gap::Real = 60.0,
        scan_intents = [XRadio.Intent(XRadio.ScanIntents.UNSPECIFIED)],
        spectral_window_intents::AbstractVector{<:AbstractString} = ["UNSPECIFIED"],
        release_date::Union{Nothing, AbstractString} = nothing,
        weight_mode::Symbol = :auto,
        weight_efficiency::Real = 1.0,
        integration_time::Union{Nothing, Real} = nothing,
        element_type::Type = Float32,
        name::AbstractString = first(splitext(basename(path))),
        zarr_format::Integer = 2,
        compressor = XRadio.Zarr.BloscCompressor(cname = "lz4", clevel = 5),
        chunks = Dict(:time => 128, :frequency => 64),
        overwrite::Bool = false,
        validate::Bool = true,
    )
    element_type <: AbstractFloat || throw(ArgumentError(
        "element_type must be a floating-point type; got $element_type"
    ))
    if release_date === nothing
        @warn "FITS-IDI carries no release date; writing the MJD zero point. \
               Pass `release_date` to state the real one." file = path
        release_date = _UNKNOWN_RELEASE_DATE
    end

    hdus = fits(path)
    primary = hdus[1]
    # An antenna number is only unique within its array, and the 256*a1 + a2
    # BASELINE packing carries no array. Two ARRAY_GEOMETRY tables would put two
    # numberings through the same map.
    narray = count(h -> card_string(h, "EXTNAME") == "ARRAY_GEOMETRY", hdus)
    narray == 1 || throw(ArgumentError(
        "FITS-IDI file defines $narray arrays; this converter handles one"
    ))
    ag = require_table(hdus, "ARRAY_GEOMETRY")
    an = require_table(hdus, "ANTENNA")
    fq = require_table(hdus, "FREQUENCY")
    sr = require_table(hdus, "SOURCE")
    uv = require_table(hdus, "UV_DATA")

    dropped = filter(t -> table(hdus, t) !== nothing, _UNMAPPED_TABLES)
    isempty(dropped) || @warn "MSv4 has no home for these FITS-IDI tables; \
                               their content is not carried into the store." dropped

    layout = uv_matrix_layout(uv)
    weightyp = card_string(uv, "WEIGHTYP", "CORRELAT")
    radiometer = weight_scaling(weightyp, weight_mode)

    # `VIS_SCAL` lets a near-real-time correlator write accumulated sums and
    # leave the normalization to the reader. Absent, it is 1. Where weights are
    # true weights it squares into them, since it divides into an amplitude.
    vis_scal = Float64(card(uv, "VIS_SCAL", 1.0))
    vis_scal > 0 || throw(ArgumentError(
        "UV_DATA VIS_SCAL is $vis_scal; it divides every visibility amplitude " *
            "and must be positive"
    ))
    weight_scal = uppercase(_trim(weightyp)) == "NORMAL" ? vis_scal^2 : 1.0
    _, products, perm = pol_axis(
        layout.no_stkd,
        Int(card(uv, "STK_1", round(Int, Float64(card(uv, "CRVAL2", -1.0))))),
        round(Int, Float64(card(uv, "CDELT2", -1.0))),
    )
    products = products[perm]

    # Frequencies. Only the first FREQID is used; a file with several carries
    # configurations this converter would silently mix into one axis.
    nfreqid = _nrows(fq)
    nfreqid == 1 || throw(ArgumentError(
        "FITS-IDI FREQUENCY table has $nfreqid frequency setups; this converter " *
            "handles one"
    ))
    fqd = fq.data
    ref_freq = Float64(card(uv, "REF_FREQ", 0.0))
    band_freq = Float64.(rowvalues(column(fqd, :BANDFREQ), 1))
    reference_frequency = reference_frequencies(ref_freq, band_freq)
    ch_width = Float64.(rowvalues(column(fqd, :CH_WIDTH), 1))
    length(band_freq) == layout.no_band || throw(ArgumentError(
        "FITS-IDI FREQUENCY gives $(length(band_freq)) bands, the UV_DATA " *
            "data matrix $(layout.no_band)"
    ))
    freqid = Int(column(fqd, :FREQID)[1])

    # Flags. FITS-IDI keeps them apart from the data as ranges to apply, so
    # nothing in `UV_DATA` says a datum is bad and the table has to be expanded
    # onto every partition's grid.
    fl = table(hdus, "FLAG")
    flag_ranges = fl === nothing ? FlagRange[] :
        read_flag_ranges(fl, perm, layout, reference_date(fl))

    # Calibration and environment measurements, each on its own clock and each
    # absent from plenty of files.
    system_calibration = read_system_calibration(hdus, layout, freqid)
    phase_calibration = read_phase_calibration(hdus, layout, freqid)
    weather = read_weather(hdus, freqid)
    # MSv4 states `GAIN_CURVE` as `Float64`. Building the curves in it keeps the
    # store the file's own numbers: the coefficients are held here against
    # radians and written back against the degrees the file gave, and in single
    # precision that scaling out and back costs a part in 1e7.
    gain_curve = msv4_gain_curves(
        read_gain_curves(hdus, layout, freqid, Float64), path
    )
    if gain_curve !== nothing
        @warn "FITS-IDI GAIN_CURVE does not say when the gain curves were \
               measured; writing the MJD zero point as `measured_date`." file = path
    end

    # Sources, by SOURCE_ID, which is what UV_DATA's SOURCE column holds.
    srd = sr.data
    source_ids = Int.(collect(column(srd, :SOURCE_ID)))
    source_names = _trim.(collect(column(srd, :SOURCE)))
    source_ra = Float64.(collect(column(srd, :RAEPO)))
    source_dec = Float64.(collect(column(srd, :DECEPO)))
    frames = source_frames(sr)
    # `UVW` takes the phase center's frame, and MSv4 admits only three there.
    # A B1950 observation has to be precessed before it can be written: naming
    # its `UVW` `"fk5"` would misplace every baseline by the 0.7 degrees
    # between the equinoxes without saying so.
    unwritable = unique(f for f in frames if !(f in XRadio.UVW_FRAMES))
    isempty(unwritable) || throw(ArgumentError(
        "FITS-IDI SOURCE gives directions in $(join(unwritable, ", ")); MSv4's " *
            "`UVW` frame is one of $(join(XRadio.UVW_FRAMES, ", ")), so no " *
            "conforming store can record them. Precess the observation to J2000 " *
            "before converting."
    ))
    source_of = Dict(source_ids[i] => i for i in eachindex(source_ids))

    bases = antenna_bases(an)
    receptor_angle = receptor_angles(an, layout.no_band)
    positions = antenna_positions(ag)
    mounts = antenna_mounts(ag)
    diameters = antenna_diameters(ag)
    agd = ag.data
    antenna_name = Dict{Int, String}(
        Int(n) => _trim(a)
            for (n, a) in zip(column(agd, :NOSTA), column(agd, :ANNAME))
    )
    # FITS-IDI names the array, not each antenna's home telescope; a VLBI
    # station belongs to the array that scheduled it.
    telescope_name = let t = card_string(ag, "ARRNAM")
        isempty(t) ? card_string(primary, "TELESCOP") : t
    end

    rows = read_uv_rows(uv)
    nrow = length(rows.jd)
    nrow > 0 || throw(ArgumentError("FITS-IDI UV_DATA table is empty"))
    times = unix_seconds.(rows.jd, rows.day)
    jd = rows.jd .+ rows.day
    inttim = _integration_times(rows.integration_time, integration_time, nrow)
    scans = _scan_ranges(rows.scan, rows.source_id, jd, gap)
    # A file's own scan numbers name the scans when it has them; otherwise the
    # inferred segments are numbered in the order they occur.
    scan_of = Vector{String}(undef, nrow)
    for (i, r) in enumerate(scans)
        scan_of[r] .= rows.scan === nothing ? string(i) : string(rows.scan[first(r)])
    end
    ranges = _partition_ranges(partition, scans, rows.source_id, jd, nrow)

    provenance = "FITS-IDI $(basename(path))"
    observer = let o = card_string(primary, "OBSERVER")
        isempty(o) ? card_string(uv, "OBSERVER") : o
    end
    visibility_units = _visibility_units(uv, path)
    shared = (;
        products, bases, antenna_name, source_names, frames, source_of, scan_of,
        source_ra, source_dec, positions, mounts, diameters, telescope_name,
        receptor_angles = receptor_angle, flag_ranges, freqid,
        system_calibration, phase_calibration, weather, gain_curve,
        gain_curve_date = _UNKNOWN_RELEASE_DATE,
        observation_span = maximum(times) - minimum(times),
        ref_freq, band_freq, reference_frequency,
        ch_width, radiometer, vis_scal, weight_scal, inttim,
        efficiency = Float64(weight_efficiency), perm, layout,
        # Stored as MSv4 states them, but validated on the way in: an intent
        # outside the ASDM enumerations cannot reach the store.
        scan_intents = [string(XRadio.Intent(i)) for i in scan_intents],
        spectral_window_intents = collect(String, spectral_window_intents),
        creation_date = _iso8601(now(UTC)),
        visibility_units,
        observation_info = Dict{Symbol, Any}(
            :observer => [observer],
            :release_date => String(release_date),
            :project_UID => card_string(primary, "OBSCODE", card_string(uv, "OBSCODE")),
        ),
        processor_info = Dict{Symbol, Any}(
            :type => "CORRELATOR",
            :sub_type => card_string(primary, "CORRELAT", "UNKNOWN"),
        ),
        provenance,
    )

    root = XRadio._create_store(
        dest, Dict{Symbol, Any}(:type => "processing_set");
        zarr_format, overwrite
    )
    index = 0
    width = max(2, ndigits(length(ranges) * layout.no_band - 1))
    for range in ranges
        block = read_uv_block(uv, range; T = element_type)
        for band in 1:layout.no_band
            ms = _build_partition(
                shared, rows, times, block, range, band, element_type
            )
            validate && XRadio.validate(ms)
            XRadio._write_tree(
                root, string(name, "_", lpad(index, width, '0')), ms;
                zarr_format, compressor, chunks
            )
            index += 1
        end
    end
    XRadio._consolidate(root)
    return dest
end

# One MSv4: the rows of `range` for one band, regridded onto the dense grid
# their distinct times and baselines imply.
function _build_partition(s, rows, times, block, range, band, ::Type{F}) where {F}
    C = Complex{F}
    layout = s.layout
    npol, nchan = layout.no_stkd, layout.no_chan

    grid_times = sort!(unique(times[range]))
    codes = sort!(unique(rows.baseline[range]))
    ntime, nbase = length(grid_times), length(codes)
    time_of = Dict(t => k for (k, t) in enumerate(grid_times))
    base_of = Dict(c => b for (b, c) in enumerate(codes))

    # MSv4 states dimensions slowest-varying first — (time, baseline_id,
    # frequency, polarization) — and a Julia array runs the other way.
    vis = fill(C(NaN, NaN), npol, nchan, nbase, ntime)
    weight = zeros(F, npol, nchan, nbase, ntime)
    flag = trues(npol, nchan, nbase, ntime)
    # `UVW` is `Float64` whatever the science arrays are. A VLBI baseline is
    # ~1e7 m and the phase it fixes turns over in a wavelength, which at 43 GHz
    # is 7 mm; `Float32` resolves that baseline to about a metre.
    uvw = fill(NaN, 3, nbase, ntime)
    field_name = Vector{String}(undef, ntime)
    scan_name = Vector{String}(undef, ntime)
    # The `FLAG` table selects by SOURCE_ID, which is a property of a timestamp
    # rather than of the grid.
    source_at = Vector{Int}(undef, ntime)

    # `WEIGHTYP = 'CORRELAT'` weights are validity fractions rather than
    # variances; `f · 2Δν · τ · η²` is the thermal weight they imply. The
    # band-constant part is hoisted, leaving the per-row integration time.
    wfactor = s.radiometer ?
        F(2 * abs(s.ch_width[band]) * s.efficiency^2) : zero(F)
    visnorm, wscal = F(inv(s.vis_scal)), F(s.weight_scal)
    percell = layout.ncomplex == 3

    for (i, r) in enumerate(range)
        k = time_of[times[r]]
        b = base_of[rows.baseline[r]]
        for axis in 1:3
            uvw[axis, b, k] = rows.uvw[r, axis] * _C_LIGHT
        end
        field_name[k] = s.source_names[s.source_of[Int(rows.source_id[r])]]
        scan_name[k] = s.scan_of[r]
        source_at[k] = Int(rows.source_id[r])
        scale = wscal * (iszero(wfactor) ? one(F) : wfactor * F(s.inttim[r]))
        for c in 1:nchan, (p, sp) in enumerate(s.perm)
            re = block.flux[1, sp, c, band, i]
            im = block.flux[2, sp, c, band, i]
            w = percell ? block.flux[3, sp, c, band, i] :
                block.weight === nothing ? one(F) : block.weight[sp, band, i]
            # Conjugated: FITS-IDI stores V = ⟨E_a1 · conj(E_a2)⟩ (AIPS Memo
            # 114r §2.1) and MSv4, following casacore, stores the other one
            # (casacore `FitsIDItoMS.cc`: "FITS-IDI convention is conjugate of
            # AIPS and CASA convention"). `UVW` and the baseline antenna order
            # are shared by both and are carried across untouched.
            vis[p, c, b, k] = C(re * visnorm, -im * visnorm)
            good = w > 0 && isfinite(re) && isfinite(im)
            weight[p, c, b, k] = good ? w * scale : zero(F)
            flag[p, c, b, k] = !good
        end
    end

    _apply_flags!(flag, s.flag_ranges, codes, grid_times, source_at, band, s.freqid, nchan)

    labels, mixed = _partition_polarization(s.products, s.bases, codes)
    sources = sort!(unique(s.source_of[Int(rows.source_id[r])] for r in range))
    partition_frames = unique(s.frames[i] for i in sources)
    length(partition_frames) == 1 || throw(ArgumentError(
        "a partition's sources are in different frames " *
            "($(join(partition_frames, ", "))); MSv4 gives one frame per dataset"
    ))
    frame = only(partition_frames)

    freqs = channel_frequencies(s.ref_freq, s.band_freq[band], s.ch_width[band], nchan)
    integration = s.inttim[first(range)]

    time = Ti(
        Sampled(
            grid_times;
            metadata = Dict{Symbol, Any}(
                :type => "time", :units => "s", :scale => "utc", :format => "unix",
                :integration_time => XRadio.Measure(
                    integration, Dict{Symbol, Any}(:units => "s", :type => "quantity")
                ),
            )
        )
    )
    frequency = Frequency(
        Sampled(
            collect(freqs);
            metadata = Dict{Symbol, Any}(
                :type => "spectral_coord", :units => "Hz", :observer => "TOPO",
                :spectral_window_name => spectral_window_name(band),

                :spectral_window_intents => s.spectral_window_intents,
                :reference_frequency => XRadio.Measure(
                    s.reference_frequency[band],
                    Dict{Symbol, Any}(
                        :units => "Hz", :type => "spectral_coord", :observer => "TOPO"
                    )
                ),
                :channel_width => XRadio.Measure(
                    abs(s.ch_width[band]),
                    Dict{Symbol, Any}(:units => "Hz", :type => "quantity")
                ),
            )
        )
    )
    polarization = Polarization(Categorical(labels))
    baseline = BaselineID(Sampled(collect(Int64, 1:nbase)))
    uvwlabel = UVWLabel(Categorical(["u", "v", "w"]))

    ms = MeasurementSet(;
        metadata = Dict{Symbol, Any}(
            :type => "visibility",
            :schema_version => string(XRadio.MSV4_SCHEMA_VERSION),
            :creation_date => s.creation_date,
            :creator => Dict{Symbol, Any}(
                :software_name => "XRadio.jl", :version => string(pkgversion(XRadio))
            ),
            :observation_info => s.observation_info,
            :processor_info => s.processor_info,
            :data_groups => Dict{Symbol, Any}(
                :base => Dict{Symbol, Any}(
                    :correlated_data => "VISIBILITY", :flag => "FLAG",
                    :weight => "WEIGHT", :uvw => "UVW",
                    :field_and_source => "field_and_source_base_xds",
                    :description => s.provenance,
                    :date => s.creation_date,
                ),
            ),
        )
    )
    grid = (polarization, frequency, baseline, time)
    ms[:visibility] = DimArray(vis, grid; metadata = Dict{Symbol, Any}(:units => s.visibility_units))
    ms[:flag] = DimArray(flag, grid)
    ms[:weight] = DimArray(weight, grid)
    ms[:uvw] = DimArray(
        uvw, (uvwlabel, baseline, time);
        metadata = Dict{Symbol, Any}(:type => "uvw", :units => "m", :frame => frame)
    )
    ms[:field_name] = DimArray(field_name, (time,))
    ms[:scan_name] = DimArray(
        scan_name, (time,);
        metadata = Dict{Symbol, Any}(:scan_intents => s.scan_intents)
    )
    antennas = map(baseline_antennas, codes)
    ms[:baseline_antenna1_name] = DimArray([s.antenna_name[a[1]] for a in antennas], (baseline,))
    ms[:baseline_antenna2_name] = DimArray([s.antenna_name[a[2]] for a in antennas], (baseline,))
    if mixed !== nothing
        ms[:polarization_mixed] = DimArray(
            reduce(hcat, mixed), (polarization, baseline)
        )
    end

    numbers = sort!(unique(Iterators.flatten(antennas)))
    branches(ms)[:antenna] = _antenna_dataset(s, numbers, band)
    branches(ms)[:field_and_source_base] = _field_and_source_dataset(s, sources, frame)

    # The calibration and environment sub-datasets, each left out where the
    # file has no such table or where nothing it holds covers this partition.
    span = (first(grid_times), last(grid_times))
    for (key, sub) in (
            :system_calibration => _system_calibration_dataset(s, numbers, span, band),
            :phase_calibration => _phase_calibration_dataset(s, numbers, span, band),
            :weather => _weather_dataset(s, numbers, span),
            :gain_curve => _gain_curve_dataset(s, numbers, band),
        )
        sub === nothing || (branches(ms)[key] = sub)
    end
    return ms
end
