# A minimal FITS-IDI file built from known values.
#
# The converter builds MSv4 datasets straight from a file's tables with no
# model in between, so a fixture file is the only thing a conversion test can
# assert against. This writer emits the five tables the conversion requires —
# ARRAY_GEOMETRY, ANTENNA, FREQUENCY, SOURCE and UV_DATA — plus FLAG and the
# four calibration and environment tables on request, and returns the dense
# arrays it wrote so a test compares against the values themselves rather than
# against a second derivation of them.
#
# It is a fixture generator, not a FITS-IDI writer: one FREQID, one array, no
# scaled columns, and every band the same width. `UV_DATA`'s column order and
# element types match the VLBA validation file, so a fixture exercises the same
# decode path as real data.

using Dates
using FITSFiles: HDU, Bintable, Card

# A column holding `n` values per row is written as a vector of vectors, except
# that FITSFiles cannot serialize a length-1 vector column — a single band has
# to be written as a scalar column instead. The reader's `rowvalues` accepts
# both shapes.
_idi_column(rows::AbstractVector{<:AbstractVector}) =
    length(first(rows)) == 1 ? [only(r) for r in rows] : rows

# The source table states its mean equinox as an `EQUINOX` string, or — as the
# VLBA writes it — as a numeric `EPOCH` the standard says to read as an equinox.
_equinox_column(::Val{:EPOCH}, equinox, n) =
    (; EPOCH = fill(equinox == "J2000" ? 2000.0 : 1950.0, n))
_equinox_column(::Val{:EQUINOX}, equinox, n) = (; EQUINOX = fill(rpad(equinox, 8), n))
_equinox_column(which::Symbol, equinox, n) = _equinox_column(Val(which), equinox, n)

# `UU`/`VV`/`WW` may carry a projection suffix, and `INTTIM` and a correlator's
# own `SCAN` column are both absent from plenty of real files, so which columns
# exist is itself a thing to test against.
function _uv_columns(
        suffix, uu, vv, ww, date, time, baseline, filter, source, freqid,
        inttim, scan, weight, flux; source_key = :SOURCE
    )
    cols = Pair{Symbol, Any}[
        Symbol("UU", suffix) => uu,
        Symbol("VV", suffix) => vv,
        Symbol("WW", suffix) => ww,
        :DATE => date, :TIME => time, :BASELINE => baseline, :FILTER => filter,
    ]
    for (name, values) in (source_key => source, :FREQID => freqid,
            :INTTIM => inttim, :SCAN => scan, :WEIGHT => weight)
        values === nothing || push!(cols, name => values)
    end
    push!(cols, :FLUX => flux)
    return NamedTuple(cols)
end

# `TUNITn` names the units of column `n`, so which cards exist depends on which
# columns do.
const _IDI_UV_UNITS = Dict(
    :DATE => "DAYS", :TIME => "DAYS", :INTTIM => "SECONDS", :FLUX => "UNCALIB"
)
# `flux_units` of `nothing` omits `FLUX`'s `TUNIT` entirely, which FITS permits.
function _tunit_cards(columns::NamedTuple, flux_units)
    cards = Card[]
    for (n, name) in enumerate(keys(columns))
        units = name === :FLUX ? flux_units :
            any(p -> startswith(String(name), p), ("UU", "VV", "WW")) ? "SECONDS" :
            get(_IDI_UV_UNITS, name, nothing)
        units === nothing || push!(cards, Card("TUNIT$n", units))
    end
    return cards
end

# The keywords AIPS Memo 114 requires on every FITS-IDI table.
function _idi_common_cards(extname::AbstractString, tabrev::Integer, g)
    return Card[
        Card("EXTNAME", extname),
        Card("EXTVER", Int32(1)),
        Card("TABREV", Int32(tabrev)),
        Card("NO_STKD", Int32(g.no_stkd)),
        Card("STK_1", Int32(g.stk_1)),
        Card("NO_BAND", Int32(g.nband)),
        Card("NO_CHAN", Int32(g.nchan)),
        Card("REF_FREQ", Float64(g.ref_freq)),
        Card("CHAN_BW", Float64(g.ch_width)),
        Card("REF_PIXL", 1.0),
        Card("OBSCODE", g.obscode),
        Card("RDATE", g.rdate),
    ]
end

"""
    synthetic_visibility_value(pol, chan, band, baseline, time) -> ComplexF32

The default fixture visibility. Every axis contributes a distinct decade, so a
transposed or mis-strided read cannot produce the expected value by accident.
"""
synthetic_visibility_value(s, c, band, b, k) =
    ComplexF32(1000k + 100b + 10band + c, s)

"""
    synthetic_weight_value(pol, band, baseline, time) -> Float32

The default fixture weight. FITS-IDI stores one weight per (polarization,
band), shared across every channel, so this takes no channel index.
"""
synthetic_weight_value(s, band, b, k) = Float32(10band + s)

"""
    synthetic_uvw_value(axis, baseline, time) -> Float32

The default fixture uvw coordinate, in **seconds of light travel**, which is
what FITS-IDI writes and what the conversion to metres has to multiply by `c`.
"""
synthetic_uvw_value(axis, b, k) = Float32((100b + k) * 1.0e-5 * axis)

"""
    synthetic_fitsidi(path; kw...) -> NamedTuple

Write a small FITS-IDI file of known content to `path` and return what was
written.

The returned named tuple carries the dense arrays behind the file —
`visibility` as `(polarization, channel, band, baseline, time)`, `weight` as
`(polarization, band, baseline, time)`, `uvw` as `(axis, baseline, time)` in
seconds — with `NaN` wherever `present` says a row was omitted, which is the
padding a dense MSv4 grid has to reproduce. It also carries the coordinates and
metadata a conversion has to recover: antenna names and positions, source
directions, channel frequencies, and `datetimes`, the timestamps in a form
independent of the file's own `DATE`/`TIME` split.

# Keywords

Array: `antenna_names`, `positions` (3 x n, metres, geocentric), `mntsta`,
`staxof` (3 x n, metres), `diameters`, `poltya`, `poltyb`, `array_name`,
`array_center`.

Sky: `source_names`, `ra` and `dec` (degrees, one per source), `equinox`.

Spectrum: `ref_freq`, `band_freq` (one offset per band), `ch_width`, `nchan`,
`sideband`, and the `STOKES` axis as `no_stkd`, `stk_1`, `cdelt2`.

Observation: `rdate`, `scans` — each a `(; source, start, ntime)` giving a
source index, a start in seconds from `rdate` midnight, and a number of
integrations — plus `dt`, `inttim`, `baselines` (`NOSTA` pairs, defaulting to
every pair including autocorrelations) and `present`, a predicate on
`(baseline_index, time_index)` selecting which rows are written.

Values: `visibility`, `weight` and `uvw`, functions of the indices above, and
`flux_type`, the float type the `FLUX` column is written in — `Float64` for a
`D`-format column.

Shape of the file itself: `uvw_suffix` names the projection on the `UU`/`VV`/`WW`
columns, `equinox_column` is `:EQUINOX` or the VLBA's numeric `:EPOCH`,
`inttim_column`, `scan_column`, `source_column`, `freqid_column` and
`weight_column` add or drop those columns, `source_key` spells the source
column `:SOURCE`, `:SOURCE_ID` or `Symbol("SOURCE ID")`, `narrays` repeats the
`ARRAY_GEOMETRY` table, and `vis_scal` writes a `VIS_SCAL` normalization the
reader has to divide out.

Provenance: `telescope`, `obscode`, `observer`, `uv_observer`, `correlator`,
`weightyp`, and `flux_units`, which is the `FLUX` column's `TUNIT` and may be
`nothing` to omit it.

`flags` adds a FLAG table, one entry per row, each a named tuple that may set
`source_id`, `ants`, `timerange`, `bands`, `chans`, `pflags`, `reason` and
`severity`.

`system_temperatures`, `phase_cals`, `weathers` and `gain_curves` add the four
optional calibration and environment tables the same way. Their times are
seconds from `rdate` midnight, as `scans` are, and `antenna` is a `NOSTA`
number — zero, in a weather row, for the whole site.

- a system temperature row sets `time`, `interval`, `antenna` and the
  `(receptor, band)` matrices `tsys` and `tant`;
- a phase-cal row sets `time`, `interval`, `antenna`, `cable`, and the
  `(tone, receptor, band)` arrays `freq` and `cal`;
- a weather row sets `time`, `interval`, `antenna` and any of `temperature`
  (centigrade), `pressure` (millibars), `dewpoint`, `wind_velocity`,
  `wind_direction` (degrees), `wvr_h2o` and `ionos_electron`;
- a gain-curve row sets `antenna`, the `(receptor, band)` integer matrices
  `type`, `nterm`, `x_typ` and `y_typ`, the `(term, receptor, band)` array
  `gain` and the `(receptor, band)` matrix `sens`. `gain_curve_terms` is the
  table's `NO_TABS`. A tabulated curve adds the `(term, receptor, band)` array
  `y_val` of the angles it is tabulated against, and a curve against hour
  angle the `(receptor, band)` matrix `x_val` of the declination it holds at;
  the slots a row does not use are `NaN`, as the standard asks.
"""
function synthetic_fitsidi(
        path;
        antenna_names = ["AA", "BB", "CC"],
        positions = [
            -2.11e6 -1.32e6 1.45e6
            -3.70e6 -4.82e6 -4.46e6
            4.55e6 3.55e6 4.30e6
        ],
        mntsta = zeros(Int32, length(antenna_names)),
        staxof = zeros(Float32, 3, length(antenna_names)),
        diameters = fill(25.0f0, length(antenna_names)),
        poltya = fill("R", length(antenna_names)),
        poltyb = fill("L", length(antenna_names)),
        array_name = "TESTARR",
        array_center = (0.0, 0.0, 0.0),
        frame = "GEOCENTRIC",
        source_names = ["SRC0", "SRC1"],
        ra = [187.7, 83.6],
        dec = [12.4, 22.0],
        equinox = "J2000",
        ref_freq = 1.0e9,
        band_freq = [0.0, 1.28e8],
        ch_width = 5.0e5,
        nchan = 4,
        sideband = 1,
        no_stkd = 4,
        stk_1 = -1,
        cdelt2 = -1,
        rdate = "2026-01-01",
        scans = [(; source = 1, start = 0.0, ntime = 3),
                 (; source = 2, start = 600.0, ntime = 2)],
        dt = 1.0,
        inttim = 1.0,
        baselines = [(a, b) for a in 1:length(antenna_names)
                         for b in a:length(antenna_names)],
        present = (b, k) -> true,
        visibility = synthetic_visibility_value,
        flux_type = Float32,
        weight = synthetic_weight_value,
        uvw = synthetic_uvw_value,
        telescope = "TEST",
        obscode = "TEST01",
        observer = "Nobody",
        uv_observer = observer,
        correlator = "TESTCORR",
        weightyp = "CORRELAT",
        flux_units = "UNCALIB",
        vis_scal = nothing,
        uvw_suffix = "---SIN",
        equinox_column = :EQUINOX,
        inttim_column = true,
        scan_column = false,
        source_column = true,
        source_key = :SOURCE,
        freqid_column = true,
        weight_column = true,
        narrays = 1,
        flags = NamedTuple[],
        system_temperatures = NamedTuple[],
        phase_cals = NamedTuple[],
        weathers = NamedTuple[],
        gain_curves = NamedTuple[],
        gain_curve_terms = 4,
    )
    nant = length(antenna_names)
    nband = length(band_freq)
    nbaseline = length(baselines)
    size(positions) == (3, nant) || throw(DimensionMismatch(
        "positions must be 3 x $nant, got $(size(positions))"
    ))
    all(1 <= s.source <= length(source_names) for s in scans) || throw(ArgumentError(
        "every scan's source must index into source_names"
    ))

    g = (; no_stkd, stk_1, nband, nchan, ref_freq, ch_width, obscode, rdate)

    # The timeline, flattened across scans. FITS-IDI splits a timestamp into a
    # whole-day `DATE` and a fractional `TIME` so the day stays exactly
    # representable; `jd0` is the constant day and `seconds` the offset into it.
    midnight = DateTime(rdate)
    jd0 = datetime2julian(midnight)
    seconds = Float64[]
    time_source = Int32[]
    time_scan = Int[]
    for (i, s) in enumerate(scans), j in 1:s.ntime
        push!(seconds, s.start + (j - 1) * dt)
        push!(time_source, Int32(s.source))
        push!(time_scan, i)
    end
    issorted(seconds) || throw(ArgumentError("scans must be given in time order"))
    ntime = length(seconds)
    datetimes = midnight .+ Millisecond.(round.(Int, seconds .* 1000))

    baseline_codes = Int32[256 * a + b for (a, b) in baselines]
    mask = BitMatrix(undef, nbaseline, ntime)
    for k in 1:ntime, b in 1:nbaseline
        mask[b, k] = present(b, k)
    end

    # The dense grids are built first and the file written from them, so the
    # returned arrays are exactly the values on disk rather than a second
    # rounding of the same formula.
    vis = fill(Complex{flux_type}(NaN, NaN), no_stkd, nchan, nband, nbaseline, ntime)
    wts = fill(NaN32, no_stkd, nband, nbaseline, ntime)
    uvws = fill(NaN32, 3, nbaseline, ntime)
    for k in 1:ntime, b in 1:nbaseline
        mask[b, k] || continue
        for band in 1:nband, c in 1:nchan, s in 1:no_stkd
            vis[s, c, band, b, k] = visibility(s, c, band, b, k)
        end
        for band in 1:nband, s in 1:no_stkd
            wts[s, band, b, k] = weight(s, band, b, k)
        end
        for axis in 1:3
            uvws[axis, b, k] = uvw(axis, b, k)
        end
    end

    # UV_DATA rows, time-major then baseline, which is the SORT = 'T*' order.
    nrow = count(mask)
    date_col = fill(jd0, nrow)
    time_col = Vector{Float64}(undef, nrow)
    bl_col = Vector{Int32}(undef, nrow)
    filter_col = zeros(Int32, nrow)
    source_col = Vector{Int32}(undef, nrow)
    freqid_col = ones(Int32, nrow)
    inttim_col = fill(Float32(inttim), nrow)
    uu_col = Vector{Float32}(undef, nrow)
    vv_col = Vector{Float32}(undef, nrow)
    ww_col = Vector{Float32}(undef, nrow)
    flux_rows = Vector{Vector{flux_type}}(undef, nrow)
    weight_rows = Vector{Vector{Float32}}(undef, nrow)
    row_time = Vector{Int}(undef, nrow)

    r = 0
    for k in 1:ntime, b in 1:nbaseline
        mask[b, k] || continue
        r += 1
        row_time[r] = k
        time_col[r] = seconds[k] / 86400
        bl_col[r] = baseline_codes[b]
        source_col[r] = time_source[k]
        uu_col[r], vv_col[r], ww_col[r] = uvws[1, b, k], uvws[2, b, k], uvws[3, b, k]
        # FLUX runs COMPLEX fastest, then STOKES, then FREQ, then BAND.
        flux = Vector{flux_type}(undef, 2 * no_stkd * nchan * nband)
        for band in 1:nband, c in 1:nchan, s in 1:no_stkd
            off = ((band - 1) * nchan * no_stkd + (c - 1) * no_stkd + (s - 1)) * 2
            v = vis[s, c, band, b, k]
            flux[off + 1] = real(v)
            flux[off + 2] = imag(v)
        end
        flux_rows[r] = flux
        # WEIGHT runs STOKES fastest, then BAND.
        weight_rows[r] = [wts[s, band, b, k] for band in 1:nband for s in 1:no_stkd]
    end

    freqs = [ref_freq + band_freq[band] + (c - 1) * ch_width
        for c in 1:nchan, band in 1:nband]

    primary = HDU(
        missing, Card[
            Card("SIMPLE", true),
            Card("BITPIX", Int32(8)),
            Card("NAXIS", Int32(0)),
            Card("EXTEND", true),
            Card("GROUPS", false),
            Card("TELESCOP", telescope),
            Card("INSTRUME", telescope),
            Card("OBSERVER", observer),
            Card("OBSCODE", obscode),
            Card("CORRELAT", correlator),
            Card("DATE-OBS", rdate),
        ]
    )

    array_geometry = HDU(
        Bintable,
        (
            ANNAME = rpad.(antenna_names, 8),
            STABXYZ = [Float64.(positions[:, i]) for i in 1:nant],
            NOSTA = Int32.(1:nant),
            MNTSTA = Int32.(mntsta),
            STAXOF = [Float32.(staxof[:, i]) for i in 1:nant],
            DIAMETER = Float32.(diameters),
        ),
        vcat(
            _idi_common_cards("ARRAY_GEOMETRY", 1, g), Card[
                Card("ARRAYX", Float64(array_center[1])),
                Card("ARRAYY", Float64(array_center[2])),
                Card("ARRAYZ", Float64(array_center[3])),
                Card("ARRNAM", array_name),
                Card("FRAME", frame),
                Card("NUMORB", Int32(0)),
                Card("FREQ", Float64(ref_freq)),
                Card("TIMSYS", "UTC"),
                Card("GSTIA0", 0.0),
                Card("DEGPDY", 360.9856),
                Card("UT1UTC", 0.0),
                Card("IATUTC", 37.0),
                Card("POLARX", 0.0),
                Card("POLARY", 0.0),
                Card("TUNIT2", "METERS"),
                Card("TUNIT5", "METERS"),
                Card("TUNIT6", "METERS"),
            ]
        )
    )

    antenna = HDU(
        Bintable,
        (
            TIME = zeros(Float64, nant),
            TIME_INTERVAL = fill(1.0f0, nant),
            ANNAME = rpad.(antenna_names, 8),
            ANTENNA_NO = Int32.(1:nant),
            ARRAY = ones(Int32, nant),
            FREQID = ones(Int32, nant),
            NO_LEVELS = fill(Int32(2), nant),
            POLTYA = collect(poltya),
            POLAA = _idi_column([zeros(Float32, nband) for _ in 1:nant]),
            POLTYB = collect(poltyb),
            POLAB = _idi_column([zeros(Float32, nband) for _ in 1:nant]),
        ),
        vcat(
            _idi_common_cards("ANTENNA", 1, g),
            Card[Card("NOPCAL", Int32(0)), Card("POLTYPE", "APPROX")]
        )
    )

    frequency = HDU(
        Bintable,
        (
            FREQID = Int32[1],
            BANDFREQ = _idi_column([Float64.(band_freq)]),
            CH_WIDTH = _idi_column([fill(Float64(ch_width), nband)]),
            TOTAL_BANDWIDTH = _idi_column([fill(Float64(ch_width * nchan), nband)]),
            SIDEBAND = _idi_column([fill(Int32(sideband), nband)]),
            BB_CHAN = _idi_column([zeros(Int32, nband)]),
        ),
        vcat(
            _idi_common_cards("FREQUENCY", 1, g),
            Card[Card("TUNIT2", "HZ"), Card("TUNIT3", "HZ"), Card("TUNIT4", "HZ")]
        )
    )

    nsource = length(source_names)
    source = HDU(
        Bintable,
        merge(
        (
            SOURCE_ID = Int32.(1:nsource),
            SOURCE = rpad.(source_names, 16),
            QUAL = zeros(Int32, nsource),
            CALCODE = fill(rpad("", 4), nsource),
            FREQID = ones(Int32, nsource),
            RAEPO = Float64.(ra),
            DECEPO = Float64.(dec),
        ),
        _equinox_column(equinox_column, equinox, nsource),
        (
            # DiFX copies the mean place into the apparent-place columns rather
            # than precessing it, so the fixture does the same.
            RAAPP = Float64.(ra),
            DECAPP = Float64.(dec),
        )),
        vcat(
            _idi_common_cards("SOURCE", 1, g), Card[
                Card("TUNIT6", "DEGREES"), Card("TUNIT7", "DEGREES"),
                Card("TUNIT9", "DEGREES"), Card("TUNIT10", "DEGREES"),
            ]
        )
    )

    uv_columns = _uv_columns(
        uvw_suffix, uu_col, vv_col, ww_col, date_col, time_col, bl_col,
        filter_col,
        source_column ? source_col : nothing,
        freqid_column ? freqid_col : nothing,
        inttim_column ? inttim_col : nothing,
        scan_column ? Int32.(time_scan[row_time]) : nothing,
        weight_column ? weight_rows : nothing, flux_rows; source_key
    )

    uv_data = HDU(
        Bintable,
        uv_columns,
        vcat(
            _idi_common_cards("UV_DATA", 1, g), Card[
                Card("NMATRIX", Int32(1)),
                Card("TMATX$(findfirst(==(:FLUX), keys(uv_columns)))", true),
                Card("MAXIS", Int32(6)),
                Card("MAXIS1", Int32(2)),
                Card("CTYPE1", "COMPLEX"),
                Card("CRVAL1", 1.0), Card("CDELT1", 1.0), Card("CRPIX1", 1.0),
                Card("MAXIS2", Int32(no_stkd)),
                Card("CTYPE2", "STOKES"),
                Card("CRVAL2", Float64(stk_1)), Card("CDELT2", Float64(cdelt2)),
                Card("CRPIX2", 1.0),
                Card("MAXIS3", Int32(nchan)),
                Card("CTYPE3", "FREQ"),
                Card("CRVAL3", Float64(ref_freq)), Card("CDELT3", Float64(ch_width)),
                Card("CRPIX3", 1.0),
                Card("MAXIS4", Int32(nband)),
                Card("CTYPE4", "BAND"),
                Card("CRVAL4", 1.0), Card("CDELT4", 1.0), Card("CRPIX4", 1.0),
                Card("MAXIS5", Int32(1)),
                Card("CTYPE5", "RA"),
                Card("CRVAL5", 0.0), Card("CDELT5", 1.0), Card("CRPIX5", 1.0),
                Card("MAXIS6", Int32(1)),
                Card("CTYPE6", "DEC"),
                Card("CRVAL6", 0.0), Card("CDELT6", 1.0), Card("CRPIX6", 1.0),
                Card("TELESCOP", telescope),
                Card("INSTRUME", telescope),
                Card("OBSERVER", uv_observer),
                Card("DATE-OBS", rdate),
                Card("EQUINOX", equinox),
                Card("WEIGHTYP", weightyp),
                (vis_scal === nothing ? Card[] :
                    Card[Card("VIS_SCAL", Float64(vis_scal))])...,
                Card("SORT", "T*"),
            ],
            _tunit_cards(uv_columns, flux_units)
        )
    )

    hdus = HDU[primary]
    append!(hdus, fill(array_geometry, narrays))
    append!(hdus, HDU[antenna, frequency, source, uv_data])
    isempty(flags) || push!(hdus, _synthetic_flag_hdu(flags, g, nchan, nband))
    isempty(system_temperatures) ||
        push!(hdus, _synthetic_syscal_hdu(system_temperatures, g))
    isempty(phase_cals) || push!(hdus, _synthetic_phasecal_hdu(phase_cals, g))
    isempty(weathers) || push!(hdus, _synthetic_weather_hdu(weathers, g))
    isempty(gain_curves) ||
        push!(hdus, _synthetic_gaincurve_hdu(gain_curves, g, gain_curve_terms))
    write(path, hdus)

    return (;
        path,
        antenna_names, nosta = collect(1:nant), positions, mntsta, staxof,
        diameters, poltya, poltyb, array_name,
        source_names, ra, dec, equinox,
        ref_freq, band_freq, ch_width, nchan, nband, no_stkd, stk_1, cdelt2,
        frequencies = freqs,
        rdate, jd0, seconds, datetimes, time_source, time_scan,
        inttim, nscan = length(scans),
        baselines, baseline_codes, present = mask,
        visibility = vis, weight = wts, uvw = uvws,
        telescope, obscode, observer, uv_observer, weightyp,
        flags, system_temperatures, phase_cals, weathers, gain_curves,
    )
end

# FLAG rows are ranges to apply rather than a per-datum array: an antenna pair
# with 0 as a wildcard, a time range in days from RDATE, a channel range, a
# per-band selector and four polarization switches.
function _synthetic_flag_hdu(flags, g, nchan::Integer, nband::Integer)
    n = length(flags)
    get_range(f, key, default) = Int32.(collect(get(f, key, default)))
    data = (
        SOURCE_ID = Int32[get(f, :source_id, 0) for f in flags],
        ARRAY = ones(Int32, n),
        ANTS = [get_range(f, :ants, (0, 0)) for f in flags],
        FREQID = ones(Int32, n),
        TIMERANG = [Float64.(collect(get(f, :timerange, (0.0, 1.0)))) for f in flags],
        BANDS = _idi_column([get_range(f, :bands, ntuple(_ -> 1, nband)) for f in flags]),
        CHANS = [get_range(f, :chans, (1, nchan)) for f in flags],
        PFLAGS = [get_range(f, :pflags, (1, 1, 1, 1)) for f in flags],
        REASON = [rpad(get(f, :reason, "TEST"), 24) for f in flags],
        SEVERITY = Int32[get(f, :severity, 0) for f in flags],
    )
    cards = vcat(
        _idi_common_cards("FLAG", 2, g),
        Card[Card("TUNIT5", "DAYS")]
    )
    return HDU(Bintable, data, cards)
end

# The calibration and environment tables, each optional and each written only
# when a test asks for it. FITS-IDI writes their times in days from `RDATE` and
# splits every per-receptor quantity into a `_1` and a `_2` column.
#
# A per-receptor column of `n` values a band is written with `n` varying
# fastest, so `(n, band)` arrays are flattened in column-major order, which is
# what `vec` already does.
_receptor_rows(::Type{T}, rows, key, p) where {T} = _idi_column(
    [T.(vec(selectdim(r[key], ndims(r[key]) - 1, p))) for r in rows]
)

# SYSTEM_TEMPERATURE rows: `tsys` and `tant` are `(receptor, band)`.
function _synthetic_syscal_hdu(rows, g)
    n = length(rows)
    data = (
        TIME = Float64[r.time / 86400 for r in rows],
        TIME_INTERVAL = Float32[get(r, :interval, 1.0) / 86400 for r in rows],
        SOURCE_ID = Int32[get(r, :source_id, 1) for r in rows],
        ANTENNA_NO = Int32[r.antenna for r in rows],
        ARRAY = ones(Int32, n),
        FREQID = ones(Int32, n),
        TSYS_1 = _receptor_rows(Float32, rows, :tsys, 1),
        TANT_1 = _receptor_rows(Float32, rows, :tant, 1),
        TSYS_2 = _receptor_rows(Float32, rows, :tsys, 2),
        TANT_2 = _receptor_rows(Float32, rows, :tant, 2),
    )
    cards = vcat(
        _idi_common_cards("SYSTEM_TEMPERATURE", 1, g),
        Card[Card("NO_POL", Int32(2)), Card("TUNIT1", "DAYS"), Card("TUNIT2", "DAYS")]
    )
    return HDU(Bintable, data, cards)
end

# PHASE-CAL rows: `freq` and `cal` are `(tone, receptor, band)`.
function _synthetic_phasecal_hdu(rows, g)
    n = length(rows)
    part(::Type{T}, f, key, p) where {T} = _idi_column(
        [T.(vec(f.(selectdim(r[key], 2, p)))) for r in rows]
    )
    data = (
        TIME = Float64[r.time / 86400 for r in rows],
        TIME_INTERVAL = Float32[get(r, :interval, 1.0) / 86400 for r in rows],
        SOURCE_ID = Int32[get(r, :source_id, 1) for r in rows],
        ANTENNA_NO = Int32[r.antenna for r in rows],
        ARRAY = ones(Int32, n),
        FREQID = ones(Int32, n),
        CABLE_CAL = Float64[r.cable for r in rows],
        PC_FREQ_1 = part(Float64, identity, :freq, 1),
        PC_REAL_1 = part(Float32, real, :cal, 1),
        PC_IMAG_1 = part(Float32, imag, :cal, 1),
        PC_FREQ_2 = part(Float64, identity, :freq, 2),
        PC_REAL_2 = part(Float32, real, :cal, 2),
        PC_IMAG_2 = part(Float32, imag, :cal, 2),
    )
    cards = vcat(
        _idi_common_cards("PHASE-CAL", 2, g),
        Card[Card("NO_POL", Int32(2)), Card("TUNIT1", "DAYS"), Card("TUNIT2", "DAYS")]
    )
    return HDU(Bintable, data, cards)
end

# WEATHER rows. `antenna` zero means the whole site, and the three atmospheric
# columns are optional in the standard, so a row that omits them is written
# without them.
function _synthetic_weather_hdu(rows, g)
    n = length(rows)
    scalar(key) = Float32[get(r, key, NaN32) for r in rows]
    data = (;
        TIME = Float64[r.time / 86400 for r in rows],
        TIME_INTERVAL = Float32[get(r, :interval, 1.0) / 86400 for r in rows],
        ANTENNA_NO = Int32[r.antenna for r in rows],
        TEMPERATURE = scalar(:temperature),
        PRESSURE = scalar(:pressure),
        DEWPOINT = scalar(:dewpoint),
        WIND_VELOCITY = scalar(:wind_velocity),
        WIND_DIRECTION = scalar(:wind_direction),
        (any(r -> haskey(r, :wvr_h2o), rows) ? (; WVR_H2O = scalar(:wvr_h2o)) : (;))...,
        (any(r -> haskey(r, :ionos_electron), rows) ?
            (; IONOS_ELECTRON = scalar(:ionos_electron)) : (;))...,
    )
    cards = vcat(
        _idi_common_cards("WEATHER", 2, g),
        Card[
            Card("TUNIT1", "DAYS"), Card("TUNIT2", "DAYS"),
            Card("TUNIT4", "CENTIGRADE"), Card("TUNIT5", "MILLIBARS"),
            Card("TUNIT6", "CENTIGRADE"), Card("TUNIT7", "M/SEC"),
            Card("TUNIT8", "DEGREES"),
        ]
    )
    return HDU(Bintable, data, cards)
end

# `X_VAL` and `Y_VAL` are the values a tabulated curve is tabulated against and
# the declination an hour-angle curve holds at, so a row that sets neither
# writes the `NaN` the standard asks for in the slots it does not use.
_gaincurve_values(::Type{T}, rows, key, p, len) where {T} = _idi_column(
    [
        haskey(r, key) ? T.(vec(selectdim(r[key], ndims(r[key]) - 1, p))) :
            fill(T(NaN), len) for r in rows
    ]
)

# GAIN_CURVE rows: `gain` is `(term, receptor, band)` and `sens`, `type`,
# `nterm` and `y_typ` are `(receptor, band)`. A curve has no time.
function _synthetic_gaincurve_hdu(rows, g, ntab::Integer)
    n = length(rows)
    ints(key, p) = _idi_column([Int32.(vec(r[key][p, :])) for r in rows])
    data = (
        ANTENNA_NO = Int32[r.antenna for r in rows],
        ARRAY = ones(Int32, n),
        FREQID = ones(Int32, n),
        TYPE_1 = ints(:type, 1), NTERM_1 = ints(:nterm, 1),
        X_TYP_1 = ints(:x_typ, 1), Y_TYP_1 = ints(:y_typ, 1),
        X_VAL_1 = _gaincurve_values(Float32, rows, :x_val, 1, g.nband),
        Y_VAL_1 = _gaincurve_values(Float32, rows, :y_val, 1, ntab * g.nband),
        GAIN_1 = _receptor_rows(Float32, rows, :gain, 1),
        SENS_1 = _idi_column([Float32.(r.sens[1, :]) for r in rows]),
        TYPE_2 = ints(:type, 2), NTERM_2 = ints(:nterm, 2),
        X_TYP_2 = ints(:x_typ, 2), Y_TYP_2 = ints(:y_typ, 2),
        X_VAL_2 = _gaincurve_values(Float32, rows, :x_val, 2, g.nband),
        Y_VAL_2 = _gaincurve_values(Float32, rows, :y_val, 2, ntab * g.nband),
        GAIN_2 = _receptor_rows(Float32, rows, :gain, 2),
        SENS_2 = _idi_column([Float32.(r.sens[2, :]) for r in rows]),
    )
    cards = vcat(
        _idi_common_cards("GAIN_CURVE", 1, g),
        Card[Card("NO_POL", Int32(2)), Card("NO_TABS", Int32(ntab))]
    )
    return HDU(Bintable, data, cards)
end
