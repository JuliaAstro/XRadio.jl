# Reading FITS-IDI (AIPS Memo 114).
#
# There is no intermediate model here. `fitsidi2msv4` builds MSv4 datasets
# straight from the file's tables, so this file holds only what that needs:
# access primitives for cards and columns, the arithmetic the format states
# rather than stores, and one bulk reader for `UV_DATA`'s per-row columns.
#
# The two access patterns exist because the file has two scales. Every table
# but `UV_DATA` is at most a few thousand rows and is read column by column
# through FITSFiles. `UV_DATA` is the observation itself — hundreds of
# thousands of rows of tens of kilobytes each — and its per-row columns all sit
# in the leading bytes of a row, so those are read in a single strided pass.
# Column-at-a-time access there would re-traverse the whole file once per
# column.

# ── Cards and columns ────────────────────────────────────────────────────────

# FITS pads strings two different ways and neither is visible in output: card
# values are space-padded to their field width, table string columns are
# NUL-padded to theirs. Untrimmed, both reach MSv4 coordinate labels and every
# name comparison in between.
_trim(s::AbstractString) = String(strip(s, ('\0', ' ')))
_trim(s) = _trim(string(s))

# The lazy descriptor behind an HDU. `hdu.data` materializes on some FITSFiles
# versions, which for `UV_DATA` means reading the entire observation; the
# descriptor carries the byte layout without touching the data block.
_descriptor(hdu) = getfield(hdu, :data)

_nrows(hdu) = Int(_descriptor(hdu).format.shape[2])
_rowlength(hdu) = Int(_descriptor(hdu).format.shape[1])

# The named HDU, by trimmed `EXTNAME`, or `nothing` when the file omits it.
# FITS-IDI makes several tables optional.
function table(hdus, name::AbstractString)
    i = findfirst(hdus) do h
        haskey(h.cards, "EXTNAME") && _trim(h.cards["EXTNAME"]) == name
    end
    return i === nothing ? nothing : hdus[i]
end

function require_table(hdus, name::AbstractString)
    t = table(hdus, name)
    t === nothing && throw(ArgumentError("FITS-IDI file has no $name table"))
    return t
end

card(hdu, key, default) = haskey(hdu.cards, key) ? hdu.cards[key] : default
card_string(hdu, key, default = "") = _trim(card(hdu, key, default))

# A `LazyStructuredData` answers `keys` but not `haskey`, so optional columns
# have to be guarded against the key list.
hascol(data, name::Symbol) = name in keys(data)
column(data, name::Symbol) = getproperty(data, name)
column(data, name::Symbol, default) = hascol(data, name) ? column(data, name) : default

# A column of `leng` elements per row comes back as an `(nrows, leng)` matrix,
# but a `leng == 1` column comes back as a plain vector. `rowvalues` gives row
# `r` as a vector either way.
rowvalues(m::AbstractMatrix, r::Integer) = vec(m[r, :])
rowvalues(v::AbstractVector, r::Integer) = [v[r]]

# ── Arithmetic the format states rather than stores ──────────────────────────

"""
    IDI_MOUNTS

Maps a FITS-IDI `MNTSTA` mount code to its name. An antenna whose code is
absent from this table cannot be described in MSv4, whose `mount` field is a
name rather than a number.
"""
const IDI_MOUNTS = Dict(
    0 => "ALT-AZ",
    1 => "EQUATORIAL",
    2 => "ORBITING",
    3 => "X-Y",
    4 => "ALT-AZ+NASMYTH-R",
    5 => "ALT-AZ+NASMYTH-L",
    6 => "OTHER",
)

"""
    IDI_POL_PRODUCTS

Maps an AIPS polarization code, as it appears on the FITS-IDI `STOKES` axis, to
the correlation product it names — but *generically*, by receptor index rather
than by polarization: `"PP"` correlates both antennas' first receptor, `"PQ"`
the first with the second, and so on.

The circular block (-1 to -4) and the linear block (-5 to -8) collapse to the
same four products, because which polarization a receptor actually carries is a
property of the antenna, not of the correlator axis: FITS-IDI states it
per-antenna in `POLTYA`/`POLTYB`. Naming the axis `"RR"` globally is only
correct when every antenna shares a basis, and a VLBI array need not — mixing
circular- and linear-feed stations is ordinary. [`pol_labels`](@ref) resolves a
product to actual polarizations for one baseline, which is the level at which
the question has an answer.
"""
const IDI_POL_PRODUCTS = Dict(
    -1 => "PP", -2 => "QQ", -3 => "PQ", -4 => "QP",
    -5 => "PP", -6 => "QQ", -7 => "PQ", -8 => "QP",
)

# Parallel hands first, then the two cross hands. AIPS interleaves them as
# PP, QQ, PQ, QP; ordering by this rank turns one into the other, and a file
# carrying only a subset keeps the same relative order.
const _IDI_PRODUCT_RANK = Dict("PP" => 1, "PQ" => 2, "QP" => 3, "QQ" => 4)

"""
    pol_axis(no_stkd, stk_1, cdelt) -> (codes, products, permutation)

Decode a FITS-IDI `STOKES` axis into AIPS codes, generic correlation products,
and the permutation that puts the parallel hands first.

The axis is described arithmetically rather than enumerated: code `k` is
`stk_1 + (k - 1) * cdelt` for `k` in `1:no_stkd`. Products come back in the
file's own order; `products[permutation]` is the canonical order.

Errors on a code outside the AIPS correlation vocabulary. Codes 1 to 4 name
Stokes parameters I, Q, U and V, which are combinations rather than
correlations and have no receptor pattern, so a file carrying them is rejected
here rather than mislabelled.

    pol_axis(4, -1, -1) == ([-1, -2, -3, -4], ["PP", "QQ", "PQ", "QP"], [1, 3, 4, 2])
"""
function pol_axis(no_stkd::Integer, stk_1::Integer, cdelt::Integer)
    codes = [Int(stk_1) + (k - 1) * Int(cdelt) for k in 1:no_stkd]
    products = map(codes) do c
        haskey(IDI_POL_PRODUCTS, c) || throw(ArgumentError(
            "FITS-IDI STOKES axis: $c is not an AIPS correlation code (-1 to -8)"
        ))
        IDI_POL_PRODUCTS[c]
    end
    allunique(products) || throw(ArgumentError(
        "FITS-IDI STOKES axis repeats a correlation product: $products"
    ))
    perm = sortperm(products; by = p -> _IDI_PRODUCT_RANK[p])
    return codes, products, perm
end

"""
    pol_labels(products, basis1, basis2) -> Vector{String}

Resolve generic correlation products into the polarizations one baseline
actually correlates, given each antenna's receptor pair from `POLTYA`/`POLTYB`.

For two circular-feed antennas this returns the familiar `["RR", "LL", "RL",
"LR"]`. For a baseline between a circular-feed and a linear-feed station it
returns something like `["RX", "LY", "RY", "LX"]`, which is what that baseline
measures and which no array-wide axis label can express. The first letter is
`basis1`'s receptor, the second `basis2`'s.

    pol_labels(["PP", "QQ", "PQ", "QP"], ("R", "L"), ("X", "Y")) ==
        ["RX", "LY", "RY", "LX"]
"""
function pol_labels(products, basis1::NTuple{2, <:AbstractString},
        basis2::NTuple{2, <:AbstractString})
    return map(products) do p
        basis1[p[1] == 'P' ? 1 : 2] * basis2[p[2] == 'P' ? 1 : 2]
    end
end

"""
    baseline_antennas(code) -> (antenna1, antenna2)

Decode a FITS-IDI `BASELINE` value into its two `NOSTA` antenna numbers.

The column packs them as `256 * antenna1 + antenna2`, which only round-trips
for numbers below 256; arrays larger than that use a separate convention this
does not implement, so a code implying an out-of-range antenna is an error.

    baseline_antennas(258) == (1, 2)
"""
function baseline_antennas(code::Integer)
    a1, a2 = divrem(Int(code), 256)
    (1 <= a1 <= 255 && 1 <= a2 <= 255) || throw(ArgumentError(
        "FITS-IDI BASELINE $code does not decode to two antennas in 1:255"
    ))
    return a1, a2
end

"""
    channel_frequencies(ref_freq, band_freq, ch_width, no_chan) -> Vector{Float64}

The sky frequency of every channel of one band, in hertz.

`band_freq` is the band's offset from `ref_freq`, both from the `FREQUENCY`
table. A lower-sideband band carries a negative `ch_width`, so its channels run
down from the band edge.

    channel_frequencies(1e9, 0.0, 1e6, 4) == [1.0e9, 1.001e9, 1.002e9, 1.003e9]
"""
channel_frequencies(ref_freq, band_freq, ch_width, no_chan::Integer) =
    ref_freq .+ band_freq .+ (0:(no_chan - 1)) .* ch_width

"""
    scan_ranges(source_id, jd; gap = 60.0) -> Vector{UnitRange{Int}}

Segment time-ordered rows into scans.

FITS-IDI records no scan number, so a scan has to be inferred: it ends where
the source changes, or where more than `gap` seconds pass between consecutive
timestamps. `jd` is in days and must already be in time order, which the `SORT`
keyword asserts and [`read_uv_rows`](@ref) does not enforce.
"""
function scan_ranges(source_id, jd; gap::Real = 60.0)
    n = length(jd)
    length(source_id) == n || throw(DimensionMismatch(
        "source_id and jd must match: $(length(source_id)) vs $n"
    ))
    n == 0 && return UnitRange{Int}[]
    gap_days = gap / 86400
    out = UnitRange{Int}[]
    start = 1
    for i in 2:n
        if source_id[i] != source_id[start] || (jd[i] - jd[i - 1]) > gap_days
            push!(out, start:(i - 1))
            start = i
        end
    end
    push!(out, start:n)
    return out
end

# ── UV_DATA's per-row columns ────────────────────────────────────────────────

# The strided reader decodes raw bytes, so it applies neither the TSCAL/TZERO
# scaling nor the multi-element layout that column access would. Both failures
# are silent, so reject a column carrying either.
function _scalar(field)
    scale = something(field.scale, 1)
    zero = something(field.zero, 0)
    (scale == 1 && zero == 0) || throw(ArgumentError(
        "UV_DATA column $(rstrip(String(field.name))) is scaled " *
            "(TSCAL $scale, TZERO $zero), which this reader does not apply"
    ))
    field.leng == 1 || throw(ArgumentError(
        "UV_DATA column $(rstrip(String(field.name))) holds $(field.leng) " *
            "values per row; this reader decodes scalar columns"
    ))
    return field
end

# One strided pass collecting every row's leading `prefix` bytes into a single
# buffer — tens of megabytes for a file of tens of gigabytes. Decoding then
# runs over memory, one type-stable loop per column.
function _read_row_prefixes(desc, prefix::Int, n::Int, rowlen::Int)
    buf = Vector{UInt8}(undef, n * prefix)
    io = FITSFiles.open_lazy_source(desc)
    try
        for r in 1:n
            seek(io, desc.begpos + rowlen * (r - 1))
            unsafe_read(io, pointer(buf, (r - 1) * prefix + 1), prefix)
        end
    finally
        close(io)
    end
    return buf
end

# FITS is big-endian; `S` is the type the file stores, `T` the type wanted.
function _decode_column!(
        out::Vector{T}, ::Type{S}, buf::Vector{UInt8}, prefix::Int, off::Int
    ) where {T, S}
    for r in eachindex(out)
        out[r] = T(ntoh(unsafe_load(Ptr{S}(pointer(buf, (r - 1) * prefix + off)))))
    end
    return out
end

function _row_column(::Type{T}, field, buf::Vector{UInt8}, prefix::Int, n::Int) where {T}
    out = Vector{T}(undef, n)
    off = first(field.slice)
    S = field.type
    S === Float64 && return _decode_column!(out, Float64, buf, prefix, off)
    S === Float32 && return _decode_column!(out, Float32, buf, prefix, off)
    S === Int64 && return _decode_column!(out, Int64, buf, prefix, off)
    S === Int32 && return _decode_column!(out, Int32, buf, prefix, off)
    S === Int16 && return _decode_column!(out, Int16, buf, prefix, off)
    throw(ArgumentError(
        "UV_DATA column $(rstrip(String(field.name))) has type $S, which is " *
            "not a scalar this reader decodes"
    ))
end

"""
    read_uv_rows(uv_hdu) -> NamedTuple

Read the per-row columns of a `UV_DATA` table — the timestamps, baselines and
uvw coordinates, but not the visibilities.

Returns `(; jd, day, baseline, source_id, freqid, integration_time, uvw)`, one
entry per row. `jd` and `day` are the timestamp FITS-IDI splits across `DATE`
and `TIME` so that the whole-day part stays exactly representable; adding them
costs precision, so callers do it once. `uvw` is `(nrows, 3)` in **seconds of
light travel**, as the file writes it, and `integration_time` is in seconds.
`baseline` is the packed code that [`baseline_antennas`](@ref) unpacks.

One pass over the file, reading only the leading bytes of each row.
"""
function read_uv_rows(uv_hdu)
    desc = _descriptor(uv_hdu)
    n = _nrows(uv_hdu)
    rowlen = _rowlength(uv_hdu)

    byname = Dict(Symbol(rstrip(String(f.name))) => f for f in desc.fields)
    function field(name::Symbol)
        haskey(byname, name) ||
            throw(ArgumentError("UV_DATA has no $name column"))
        return _scalar(byname[name])
    end

    # UU/VV/WW carry a suffix naming the projection the correlator used
    # ("UU---SIN"); the standard also allows the bare names.
    function uvw_field(axis::Char)
        base = string(axis, axis)
        for name in keys(byname)
            startswith(String(name), base) && return _scalar(byname[name])
        end
        throw(ArgumentError("UV_DATA has no $base column"))
    end

    wanted = [
        field(:DATE), field(:TIME), field(:BASELINE), field(:SOURCE),
        field(:FREQID), field(:INTTIM),
        uvw_field('U'), uvw_field('V'), uvw_field('W'),
    ]
    prefix = maximum(f -> last(f.slice), wanted)
    buf = _read_row_prefixes(desc, prefix, n, rowlen)
    col(::Type{T}, f) where {T} = _row_column(T, f, buf, prefix, n)

    uvw = Matrix{Float64}(undef, n, 3)
    for (k, f) in enumerate(wanted[7:9])
        uvw[:, k] = col(Float64, f)
    end

    return (;
        jd = col(Float64, wanted[1]),
        day = col(Float64, wanted[2]),
        baseline = col(Int32, wanted[3]),
        source_id = col(Int32, wanted[4]),
        freqid = col(Int32, wanted[5]),
        integration_time = col(Float64, wanted[6]),
        uvw,
    )
end
