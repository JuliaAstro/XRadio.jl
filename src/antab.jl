"""
    read_antab!(ps, path; ifs = nothing) -> ps

Read the system temperatures and gain curves of the ANTAB file at `path` into
every Measurement Set of `ps`, as the standard `system_calibration` (`TSYS`)
and `gain_curve` datasets. [`system_temperatures`](@ref) and
[`gaincurves`](@ref) read them back.

ANTAB is the a-priori calibration format of AIPS's `ANTAB` task. A `TSYS`
block's `INDEX` names the AIPS IF each column applies to (`'R1:32'` is one value
for right circular polarization in IFs 1 to 32, `'L1|R1'` one value for both
receptors of IF 1, `'X'` an unused column). An IF is a FITS-IDI band, which
[`fitsidi2msv4`](@ref) writes as one spectral window and records as the window's
`fitsidi_band`; for a store without that attribute, `ifs` maps each spectral
window name to its IF, such as `Dict("spw_0" => 1, "spw_1" => 2)`, or to the IF
of each of its channels, such as `Dict("spw_0" => 1:32)` for a UVFITS file whose
32 one-channel IFs were read as the channels of one window. A window given one
IF per channel takes a system temperature per channel, over `frequency`. The
IFs the file covers and the IFs of `ps`'s windows must be the same.

Stations are matched to antennas by name, and polarizations to receptors by the
antenna's `polarization_type`. Each Measurement Set takes the rows that fall
within its own time span, on the file's own clock (`TimeSystemCal`); an antenna
the file says nothing about is `NaN`. ANTAB rows give a day of the year and no
year, which is taken as the one that puts each row nearest the set's first
sample. A `GAIN` record holds where its `FREQ` range, if it states one, contains
the window's reference frequency; `ELEV` curves are polynomials in elevation and
`ALTAZ` curves in zenith angle, both in degrees, with `DPFU` as the sensitivity
in kelvin per jansky (right, then left, polarization when it gives two). ANTAB
does not say when a curve was measured, so `measured_date` is the MJD zero
point.

A Measurement Set that already has either dataset is an error rather than being
overwritten; `delete!(branches(ms), :system_calibration)` removes one.
"""
function read_antab!(ps::ProcessingSet, path::AbstractString; ifs = nothing)
    antab = _parse_antab(path)
    members = collect(values(ps))
    isempty(members) && throw(ArgumentError("the processing set holds no Measurement Sets"))
    for ms in members, key in (:system_calibration, :gain_curve)
        haskey(branches(ms), key) && throw(ArgumentError(
            "spectral window `$(spectralwindow(ms))` already has a `$key` dataset, and " *
                "ANTAB data are written only where there is none; " *
                "`delete!(branches(ms), :$key)` removes it"
        ))
    end

    windows = [_antab_if(ms, ifs) for ms in members]
    _check_antab_ifs(antab, Set(Iterators.flatten(windows)))
    spans = [_antab_span(ms) for ms in members]
    t0 = minimum(first, spans)
    interval = maximum(last, spans) - t0
    clock = Dict(s => _antab_unix.(b.doy, b.seconds, t0) for (s, b) in antab.tsys)

    named(a) = haskey(antab.tsys, a) || haskey(antab.gains, a)
    any(named, (a for ms in members for a in antennas(ms))) || throw(ArgumentError(
        "$path names stations $(join(sort!(collect(union(keys(antab.tsys), keys(antab.gains)))), ", ")), " *
            "and no antenna of the processing set has any of those names"
    ))

    written = map(members, windows, spans) do ms, k, span
        roster = _antab_identity(ms)
        (;
            system_calibration = _antab_tsys_dataset(antab, clock, roster, k, span, ms),
            gain_curve = _antab_gain_dataset(antab, roster, _reference_frequency(ms), interval),
        )
    end
    for (ms, sub) in zip(members, written), (key, xds) in pairs(sub)
        isnothing(xds) || (ms[key] = xds)
    end
    any(sub -> !isnothing(sub.gain_curve), written) && @warn(
        "ANTAB does not say when the gain curves were measured; writing the MJD " *
            "zero point as `measured_date`.", file = path
    )
    return ps
end

function _antab_if(ms::MeasurementSet, ifs)
    window = spectralwindow(ms)
    if isnothing(ifs)
        band = get(metadata(lookup(dims(ms, dimtype(:frequency)))), :fitsidi_band, nothing)
        isnothing(band) && throw(ArgumentError(
            "spectral window `$window` does not record its FITS-IDI band (`fitsidi_band`), " *
                "so nothing says which ANTAB IF it is; pass `ifs = Dict(\"$window\" => k, ...)`"
        ))
        return Int(band)
    end
    haskey(ifs, window) || throw(ArgumentError("`ifs` gives no IF for spectral window `$window`"))
    k = ifs[window]
    k isa Integer && return Int(k)
    nchan = length(frequencies(ms))
    length(k) == nchan || throw(ArgumentError(
        "`ifs` gives $(length(k)) IFs for spectral window `$window`, which has $nchan channels"
    ))
    return collect(Int, k)
end

function _check_antab_ifs(antab, windows)
    covered = Set{Int}()
    for block in values(antab.tsys), column in block.columns, (_, range) in column
        union!(covered, range)
    end
    isempty(covered) && return nothing
    covered == windows || throw(ArgumentError(
        "$(antab.path) covers IFs $(sort!(collect(covered))) and the processing set's " *
            "spectral windows are IFs $(sort!(collect(windows)))"
    ))
    return nothing
end

function _antab_span(ms::MeasurementSet)
    meta = metadata(lookup(dims(ms, dimtype(:time))))
    scale, format = get(meta, :scale, nothing), get(meta, :format, nothing)
    (scale, format) == ("utc", "unix") || throw(ArgumentError(
        "ANTAB times are UTC, and spectral window `$(spectralwindow(ms))` states its " *
            "time as $scale $format rather than utc unix"
    ))
    lo, hi = extrema(times(ms))
    half = haskey(meta, :integration_time) ? meta[:integration_time].value / 2 : 0.0
    return (lo - half, hi + half)
end

function _antab_unix(doy, seconds, t0)
    y = Dates.year(unix2datetime(t0))
    candidates = (datetime2unix(DateTime(yy)) + 86400 * (doy - 1) + seconds for yy in (y - 1):(y + 1))
    return argmin(t -> abs(t - t0), candidates)
end

_reference_frequency(ms::MeasurementSet) =
    metadata(lookup(dims(ms, dimtype(:frequency))))[:reference_frequency].value

# What an antenna-keyed sub-dataset repeats of the antenna dataset, with
# `polarization_type` laid out (receptor, antenna) as the builder takes it.
function _antab_identity(ms::MeasurementSet)
    haskey(branches(ms), :antenna) || throw(ArgumentError(
        "spectral window `$(spectralwindow(ms))` has no antenna dataset to match ANTAB stations against"
    ))
    ax = branches(ms)[:antenna]
    A, R = dimtype(:antenna_name), dimtype(:receptor_label)
    names, receptors = collect(lookup(ax, A)), collect(lookup(ax, R))
    poltype = ax[:polarization_type]
    return (;
        antenna_name = names,
        station_name = collect(ax[:station_name]),
        mount = collect(ax[:mount]),
        telescope_name = collect(ax[:telescope_name]),
        receptor_label = receptors,
        polarization_type = [poltype[R(j), A(i)] for j in eachindex(receptors), i in eachindex(names)],
    )
end

function _antab_receptor(roster, i, letter, station)
    j = findall(==(string(letter)), view(roster.polarization_type, :, i))
    length(j) == 1 || throw(ArgumentError(
        "ANTAB station $station gives polarization $letter, and the antenna's receptors are " *
            join(roster.polarization_type[:, i], ", ")
    ))
    return only(j)
end

function _antab_tsys_dataset(antab, clock, roster, k, span, ms)
    rows = [
        (station, i, r) for (i, station) in pairs(roster.antenna_name) if haskey(antab.tsys, station)
            for r in eachindex(clock[station]) if span[1] <= clock[station][r] <= span[2]
    ]
    isempty(rows) && return nothing
    grid = sort!(unique(clock[s][r] for (s, _, r) in rows))
    at = Dict(t => n for (n, t) in pairs(grid))
    channel_ifs = k isa Integer ? [k] : k
    tsys = fill(NaN, length(channel_ifs), length(roster.receptor_label), length(grid), length(roster.antenna_name))
    for (station, i, r) in rows
        block = antab.tsys[station]
        for (c, column) in pairs(block.columns), (letter, range) in column, (n, kn) in pairs(channel_ifs)
            kn in range || continue
            j = _antab_receptor(roster, i, letter, station)
            tsys[n, j, at[clock[station][r]], i] = block.values[r, c]
        end
    end
    all(isnan, tsys) && return nothing
    (; antenna_name, station_name, receptor_label, polarization_type) = roster
    k isa Integer && return subdataset(
        :system_calibration;
        antenna_name, station_name, receptor_label, polarization_type,
        time_system_cal = grid, tsys = dropdims(tsys; dims = 1),
    )
    F = dims(ms, dimtype(:frequency))
    return subdataset(
        :system_calibration;
        antenna_name, station_name, receptor_label, polarization_type,
        time_system_cal = grid, frequency = collect(F), tsys,
        array_metadata = (; frequency = metadata(lookup(F))),
    )
end

function _antab_gain_dataset(antab, roster, reference_frequency, interval)
    curves = map(roster.antenna_name) do station
        held = filter(get(antab.gains, station, _AntabGain[])) do g
            isnothing(g.freq) || g.freq[1] <= reference_frequency <= g.freq[2]
        end
        length(held) > 1 && throw(ArgumentError(
            "$(antab.path) has $(length(held)) GAIN records for $station that hold at " *
                "$(reference_frequency / 1.0e6) MHz"
        ))
        isempty(held) ? nothing : only(held)
    end
    all(isnothing, curves) && return nothing

    nreceptor, nant = length(roster.receptor_label), length(roster.antenna_name)
    nterm = maximum(length(c.poly) for c in curves if !isnothing(c))
    coefficients = zeros(nreceptor, nterm, nant)
    sensitivity = fill(NaN, nreceptor, nant)
    types = fill("POWER(EL)", nant)
    for (i, curve) in pairs(curves)
        isnothing(curve) && continue
        types[i] = curve.kind
        for j in 1:nreceptor
            coefficients[j, eachindex(curve.poly), i] = curve.poly
            sensitivity[j, i] = curve.dpfu[_dpfu_slot(curve, roster.polarization_type[j, i])]
        end
    end
    return subdataset(
        :gain_curve; roster...,
        poly_term = collect(Int64, 0:(nterm - 1)),
        gain_curve_type = types,
        gain_curve = coefficients,
        gain_curve_sensitivity = sensitivity,
        gain_curve_interval = fill(interval, nant),
        metadata = (; measured_date = "1858-11-17T00:00:00.000"),
    )
end

function _dpfu_slot(curve, polarization)
    length(curve.dpfu) == 1 && return 1
    polarization in ("R", "X") && return 1
    polarization in ("L", "Y") && return 2
    throw(ArgumentError("a DPFU pair is right then left polarization, and this receptor is $polarization"))
end

# The text of an ANTAB file, before anything is matched to a Measurement Set.
# Times stay day of year and seconds of that day, since the file states no year.

# One TSYS column: the receptors and IFs its value applies to, so `'R1:32'` is
# `[('R', 1:32)]`, `'L1|R1'` is `[('L', 1:1), ('R', 1:1)]` and `'X'` is empty.
const _AntabColumn = Vector{Tuple{Char, StepRange{Int, Int}}}

struct _AntabGain
    kind::String
    dpfu::Vector{Float64}
    poly::Vector{Float64}
    freq::Union{Nothing, Tuple{Float64, Float64}}
end

struct _AntabTsys
    columns::Vector{_AntabColumn}
    doy::Vector{Int}
    seconds::Vector{Float64}
    values::Matrix{Float64}
end

struct _Antab
    path::String
    gains::Dict{String, Vector{_AntabGain}}
    tsys::Dict{String, _AntabTsys}
end

_antab_error(path, lineno, msg) = throw(ArgumentError("$path, line $lineno: $msg"))

_strip_antab_comment(line) = first(split(line, '!'; limit = 2))

function _parse_antab(path::AbstractString)
    lines = String.(rstrip.(_strip_antab_comment.(readlines(path))))
    gains = Dict{String, Vector{_AntabGain}}()
    tsys = Dict{String, _AntabTsys}()
    i = firstindex(lines)
    while i <= lastindex(lines)
        line = strip(lines[i])
        keyword = isempty(line) ? "" : uppercase(first(split(line)))
        if keyword == "GAIN"
            station, gain = _parse_antab_gain(path, i, line)
            push!(get!(gains, station, _AntabGain[]), gain)
            i += 1
        elseif keyword == "TSYS"
            station, block, next = _parse_antab_tsys(path, lines, i)
            haskey(tsys, station) && _antab_error(path, i, "a second TSYS block for $station")
            tsys[station] = block
            i = next
        else
            i += 1
        end
    end
    for station in keys(tsys)
        haskey(gains, station) || throw(ArgumentError(
            "$path has a TSYS block for $station and no GAIN record to scale it by"
        ))
    end
    return _Antab(String(path), gains, tsys)
end

const _ANTAB_GAIN_KINDS = Dict("ELEV" => "POWER(EL)", "ALTAZ" => "POWER(ZA)")

function _parse_antab_gain(path, lineno, line)
    tokens = split(replace(line, '/' => ' '))
    length(tokens) >= 3 || _antab_error(path, lineno, "a GAIN record names a station and a curve type")
    station, kind = String(tokens[2]), uppercase(tokens[3])
    haskey(_ANTAB_GAIN_KINDS, kind) || _antab_error(path, lineno,
        "GAIN $station is of type $kind; only ELEV and ALTAZ, polynomials in elevation " *
            "and in zenith angle, are read"
    )
    rest = join(tokens[4:end], ' ')
    value(key) = match(Regex("$key\\s*=\\s*([-+0-9.eEdD,\\s]+?)(?=\\s+[A-Za-z]+\\s*=|\\s*\$)", "i"), rest)
    dpfu, poly, freq = value("DPFU"), value("POLY"), value("FREQ")
    isnothing(dpfu) && _antab_error(path, lineno, "GAIN $station has no DPFU")
    isnothing(poly) && _antab_error(path, lineno, "GAIN $station has no POLY")
    dpfus, coefficients = _antab_numbers(dpfu[1]), _antab_numbers(poly[1])
    length(dpfus) in (1, 2) || _antab_error(path, lineno,
        "GAIN $station gives $(length(dpfus)) DPFU values; ANTAB states one, or one per receptor"
    )
    isempty(coefficients) && _antab_error(path, lineno, "GAIN $station has an empty POLY")
    range = nothing
    if !isnothing(freq)
        mhz = _antab_numbers(freq[1])
        length(mhz) == 2 || _antab_error(path, lineno, "GAIN $station's FREQ is a lower and an upper frequency")
        range = (1.0e6 * mhz[1], 1.0e6 * mhz[2])
    end
    return station, _AntabGain(_ANTAB_GAIN_KINDS[kind], dpfus, coefficients, range)
end

_antab_numbers(s) = [parse(Float64, replace(x, r"[dD]" => "e")) for x in split(s, [',', ' ', '\t']; keepempty = false)]

# A TSYS block is a header, which may run over several lines and ends at the
# first `/`, then rows of `DOY TIME value...` up to the next `/`.
function _parse_antab_tsys(path, lines, start)
    station = String(split(strip(lines[start]))[2])
    header = IOBuffer()
    i = start
    rest = join(split(strip(lines[i]))[3:end], ' ')
    while true
        slash = findfirst('/', rest)
        if !isnothing(slash)
            print(header, ' ', rest[1:prevind(rest, slash)])
            break
        end
        print(header, ' ', rest)
        i += 1
        i <= lastindex(lines) || _antab_error(path, start, "TSYS $station has no `/` closing its header")
        rest = strip(lines[i])
    end
    text = String(take!(header))
    occursin(r"INDEX2\s*="i, text) && _antab_error(path, start, "TSYS $station uses INDEX2, which is not read")
    index = match(r"INDEX\s*=\s*((?:'[^']*'\s*,?\s*)+)"i, text)
    isnothing(index) && _antab_error(path, start, "TSYS $station has no INDEX")
    columns = [_parse_antab_index(path, start, m[1]) for m in eachmatch(r"'([^']*)'", index[1])]

    doy, seconds, rows = Int[], Float64[], Vector{Float64}[]
    i += 1
    while i <= lastindex(lines)
        line = strip(lines[i])
        if startswith(line, '/')
            values = Matrix{Float64}(undef, length(rows), length(columns))
            for (r, row) in pairs(rows)
                values[r, :] = row
            end
            return station, _AntabTsys(columns, doy, seconds, values), i + 1
        end
        if !isempty(line)
            tokens = split(line)
            length(tokens) == 2 + length(columns) || _antab_error(path, i,
                "TSYS $station's INDEX names $(length(columns)) columns, and this row has $(length(tokens) - 2) values"
            )
            push!(doy, parse(Int, tokens[1]))
            push!(seconds, _antab_seconds(path, i, tokens[2]))
            push!(rows, parse.(Float64, tokens[3:end]))
        end
        i += 1
    end
    return _antab_error(path, start, "TSYS $station has no `/` closing its rows")
end

# `R1`, `L1:32`, `L10:16:2`, receptors sharing one value (`L1|R1`), or `X`
# for a column that applies to nothing.
function _parse_antab_index(path, lineno, entry)
    column = Tuple{Char, StepRange{Int, Int}}[]
    uppercase(strip(entry)) == "X" && return column
    for part in split(entry, '|')
        m = match(r"^\s*([RLXY])(\d+)(?::(\d+)(?::(\d+))?)?\s*$"i, part)
        isnothing(m) && _antab_error(path, lineno, "unrecognized INDEX entry '$entry'")
        lo = parse(Int, m[2])
        hi = isnothing(m[3]) ? lo : parse(Int, m[3])
        step = isnothing(m[4]) ? 1 : parse(Int, m[4])
        push!(column, (uppercase(only(m[1])), lo:step:hi))
    end
    return column
end

# `HH:MM:SS.s`, or `HH:MM.m` with decimal minutes.
function _antab_seconds(path, lineno, s)
    parts = split(s, ':')
    if length(parts) == 3
        return 3600 * parse(Int, parts[1]) + 60 * parse(Int, parts[2]) + parse(Float64, parts[3])
    elseif length(parts) == 2
        return 3600 * parse(Int, parts[1]) + 60 * parse(Float64, parts[2])
    end
    return _antab_error(path, lineno, "unrecognized time '$s'")
end
