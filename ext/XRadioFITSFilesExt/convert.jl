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

# The antennas of one partition, in antenna-number order. A scan's baselines
# name the antennas that actually observed together, which is fewer than the
# array whenever a station was down.
function _antenna_dataset(s, codes, band::Int)
    numbers = sort!(unique(Iterators.flatten(map(baseline_antennas, codes))))
    names = [s.antenna_name[a] for a in numbers]
    mounts = [s.mounts[a] for a in numbers]

    unnamed = [n for (n, m) in zip(names, mounts) if m isa XRadio.MountOrbiting ||
        m isa XRadio.MountOther]
    isempty(unnamed) || @warn "MSv4 has no `mount` name for these antennas' \
                               mounts; writing \"OTHER\", which does not say \
                               which they are." antennas = unnamed maxlog = 1

    antenna = AntennaName(Categorical(names))
    receptor = ReceptorLabel(Categorical(["pol_0", "pol_1"]))
    cart = CartesianPosLabel(Categorical(["x", "y", "z"]))

    ds = Dataset(; metadata = Dict{Symbol, Any}(:type => "antenna"))
    ds[:antenna_position] = DimArray(
        [s.positions[a][i] for i in 1:3, a in numbers], (cart, antenna);
        metadata = Dict{Symbol, Any}(
            :type => "location", :units => "m", :frame => "ITRS",
            :coordinate_system => "geocentric", :origin_object_name => "earth",
        )
    )
    # FITS-IDI names an antenna once. A VLBI antenna does not move between
    # stations, so the one name is both.
    ds[:station_name] = DimArray(copy(names), (antenna,))
    ds[:mount] = DimArray([XRadio.msv4_mount(m) for m in mounts], (antenna,))
    ds[:telescope_name] = DimArray(fill(s.telescope_name, length(numbers)), (antenna,))
    ds[:polarization_type] = DimArray(
        [s.bases[a][i] for i in 1:2, a in numbers], (receptor, antenna)
    )
    ds[:antenna_receptor_angle] = DimArray(
        [s.receptor_angles[a][i, band] for i in 1:2, a in numbers], (receptor, antenna);
        metadata = Dict{Symbol, Any}(:type => "quantity", :units => "rad")
    )
    s.diameters === nothing || (
        ds[:antenna_dish_diameter] = DimArray(
            [s.diameters[a] for a in numbers], (antenna,);
            metadata = Dict{Symbol, Any}(:type => "quantity", :units => "m")
        )
    )
    return ds
end

# The fields of one partition. FITS-IDI states one direction per source, which
# is both where the source is and what the visibilities are phased to: the
# format has no notion of pointing somewhere other than at the source.
function _field_and_source_dataset(s, sources, frame::AbstractString)
    names = [s.source_names[i] for i in sources]
    direction = permutedims([deg2rad.(s.source_ra[sources]) deg2rad.(s.source_dec[sources])])
    attrs = Dict{Symbol, Any}(:type => "sky_coord", :units => "rad", :frame => frame)

    field = FieldName(Categorical(names))
    skylabel = SkyDirLabel(Categorical(["ra", "dec"]))

    ds = Dataset(; metadata = Dict{Symbol, Any}(:type => "field_and_source"))
    ds[:source_name] = DimArray(copy(names), (field,))
    ds[:field_phase_center_direction] =
        DimArray(direction, (skylabel, field); metadata = attrs)
    ds[:source_direction] = DimArray(copy(direction), (skylabel, field); metadata = copy(attrs))
    return ds
end

"""
    fitsidi2msv4(path, dest; kwargs...) -> dest

See [`XRadio.fitsidi2msv4`](@ref).
"""
function XRadio.fitsidi2msv4(
        path::AbstractString, dest::AbstractString;
        partition = (:scan,),
        gap::Real = 60.0,
        scan_intents::AbstractVector{<:AbstractString} = ["UNSPECIFIED"],
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
    shared = (;
        products, bases, antenna_name, source_names, frames, source_of, scan_of,
        source_ra, source_dec, positions, mounts, diameters, telescope_name,
        receptor_angles = receptor_angle, flag_ranges, freqid,
        ref_freq, band_freq, ch_width, radiometer, vis_scal, weight_scal, inttim,
        efficiency = Float64(weight_efficiency), perm, layout,
        scan_intents = collect(String, scan_intents),
        spectral_window_intents = collect(String, spectral_window_intents),
        creation_date = _iso8601(now(UTC)),
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
            vis[p, c, b, k] = C(re * visnorm, im * visnorm)
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
                :spectral_window_name => "band_$(band - 1)",
                :spectral_window_intents => s.spectral_window_intents,
                :reference_frequency => XRadio.Measure(
                    s.ref_freq + s.band_freq[band],
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
    baseline = BaselineId(Sampled(collect(Int64, 1:nbase)))
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
    ms[:visibility] = DimArray(vis, grid; metadata = Dict{Symbol, Any}(:units => "Jy"))
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

    branches(ms)[:antenna] = _antenna_dataset(s, codes, band)
    branches(ms)[:field_and_source_base] = _field_and_source_dataset(s, sources, frame)
    return ms
end
