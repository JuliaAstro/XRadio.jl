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
# thousands of rows of tens of kilobytes each — so it is decoded straight from
# the memory-mapped data unit, its per-row columns in one strided pass over the
# leading bytes of each row.

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

# The memory-mapped data unit; it becomes invalid if the file is truncated.
_bytes(hdu) = getfield(_descriptor(hdu), :bytes)

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

# The `TUNIT` of the column named `name`, or `nothing` where the column has none.
# FITS makes `TUNIT` optional, so an absent one means the units are unstated
# rather than dimensionless.
function column_units(hdu, name::AbstractString)
    i = 1
    while haskey(hdu.cards, "TTYPE$i")
        card_string(hdu, "TTYPE$i") == name && return let u = card_string(hdu, "TUNIT$i")
            isempty(u) ? nothing : u
        end
        i += 1
    end
    return nothing
end

column(data, name::Symbol) = getproperty(data, name)
column(data, name::Symbol, default) = haskey(data, name) ? column(data, name) : default

# A column of `leng` elements per row comes back as an `(nrows, leng)` matrix,
# but a `leng == 1` column comes back as a plain vector. `rowvalues` gives row
# `r` as a vector either way.
rowvalues(m::AbstractMatrix, r::Integer) = vec(m[r, :])
rowvalues(v::AbstractVector, r::Integer) = [v[r]]

# ── Arithmetic the format states rather than stores ──────────────────────────

"""
    IDI_MOUNTS

Maps a FITS-IDI `MNTSTA` mount code to the constructor of the
[`XRadio.AbstractMount`](@ref) it names, each taking the antenna's `STAXOF`
axis offset.
"""
const IDI_MOUNTS = Dict(
    0 => MountAltAz,
    1 => MountEquatorial,
    2 => MountOrbiting,
    3 => MountXY,
    4 => MountNasmythR,
    5 => MountNasmythL,
    6 => MountBWGR,
    7 => MountBWGL,
)

"""
    idi_mount(mntsta, staxof) -> AbstractMount

The mount an AIPS `MNTSTA` code names, carrying `staxof` as its axis offset.
"""
function idi_mount(mntsta::Integer, staxof)
    haskey(IDI_MOUNTS, Int(mntsta)) || throw(ArgumentError(
        "FITS-IDI MNTSTA $mntsta names no mount; AIPS defines " *
            join(sort(collect(keys(IDI_MOUNTS))), ", ")
    ))
    return IDI_MOUNTS[Int(mntsta)](staxof)
end

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

_fieldname(field) = rstrip(String(field.name))

# The byte readers here decode raw bytes, so they apply none of the TSCAL/TZERO
# scaling that column access would. The failure is silent, so reject a column
# carrying either.
function _unscaled(field)
    scale = something(field.scale, 1)
    zero = something(field.zero, 0)
    (scale == 1 && zero == 0) || throw(ArgumentError(
        "UV_DATA column $(_fieldname(field)) is scaled " *
            "(TSCAL $scale, TZERO $zero), which this reader does not apply"
    ))
    return field
end

function _scalar(field)
    _unscaled(field)
    field.leng == 1 || throw(ArgumentError(
        "UV_DATA column $(_fieldname(field)) holds $(field.leng) " *
            "values per row; this reader decodes scalar columns"
    ))
    return field
end

# FITS is big-endian; `S` is the type the file stores and `T` the type wanted.
# `out` is one column per row, `stride` the distance between rows within `buf`,
# and `off` the 1-based byte offset in `buf` of the field in the first row.
# The caller keeps every byte read within `buf`.
function _decode!(
        out::Matrix{T}, ::Type{S}, buf::Vector{UInt8}, stride::Int, off::Int
    ) where {T, S}
    GC.@preserve buf for r in axes(out, 2)
        base = (r - 1) * stride + off
        for k in axes(out, 1)
            out[k, r] = T(ntoh(unsafe_load(Ptr{S}(pointer(buf, base + (k - 1) * sizeof(S))))))
        end
    end
    return out
end

# Decodes `field` for the rows of `out`, the first of them starting at byte
# `start + 1` of `buf`.
function _decode_field!(out::Matrix, field, buf::Vector{UInt8}, stride::Int, start::Int = 0)
    off = start + first(field.slice)
    S = field.type
    S === Float64 && return _decode!(out, Float64, buf, stride, off)
    S === Float32 && return _decode!(out, Float32, buf, stride, off)
    S === Int64 && return _decode!(out, Int64, buf, stride, off)
    S === Int32 && return _decode!(out, Int32, buf, stride, off)
    S === Int16 && return _decode!(out, Int16, buf, stride, off)
    throw(ArgumentError(
        "UV_DATA column $(_fieldname(field)) has type $S, which is not a " *
            "numeric type this reader decodes"
    ))
end

_row_column(::Type{T}, field, buf::Vector{UInt8}, rowlen::Int, n::Int) where {T} =
    vec(_decode_field!(Matrix{T}(undef, 1, n), field, buf, rowlen))

"""
    read_uv_rows(uv_hdu) -> NamedTuple

Read the per-row columns of a `UV_DATA` table — the timestamps, baselines and
uvw coordinates, but not the visibilities.

Returns `(; jd, day, baseline, source_id, freqid, integration_time, scan, uvw)`,
one entry per row. `jd` and `day` are the timestamp FITS-IDI splits across
`DATE` and `TIME` so that the whole-day part stays exactly representable; adding
them costs precision, so callers do it once. `uvw` is `(nrows, 3)` in **seconds
of light travel**, as the file writes it. `baseline` is the packed code that
[`baseline_antennas`](@ref) unpacks.

Only `DATE`, `TIME`, `BASELINE` and the `UU`/`VV`/`WW` triple are required of a
FITS-IDI file. `source_id` and `freqid` fall back to 1, which is what a file
holding one source or one frequency setup means by omitting them.
`integration_time` and `scan` are `nothing` when the file does not state them —
the first is optional in the format, and the second is not part of it at all,
so a `SCAN` column is a correlator's own addition and is read when offered.

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
    # AIPS Memo 114r makes several random parameters conditional — `SOURCE_ID`
    # only when the file holds more than one source, `FREQID` only when it holds
    # more than one setup, `INTTIM` never — and states three spellings of
    # `SOURCE_ID` to be synonymous.
    function optional_field(names::Symbol...)
        for name in names
            haskey(byname, name) && return _scalar(byname[name])
        end
        return nothing
    end

    # UU/VV/WW carry a suffix naming the projection the correlator used. The
    # ---NCP system puts `w` on the north pole rather than on the source, so its
    # coordinates are not the ones MSv4 wants and reading them as ---SIN would
    # be wrong in silence. The VLBA archive's `-L` suffix means ---SIN in
    # seconds despite naming wavelengths, which the memo tells readers to expect.
    function uvw_field(axis::Char)
        base = string(axis, axis)
        for name in keys(byname)
            full = String(name)
            startswith(full, base) || continue
            occursin("NCP", full) && throw(ArgumentError(
                "UV_DATA column $full is in the ---NCP system, whose w axis points " *
                    "at the north pole rather than at the source; this reader " *
                    "handles ---SIN"
            ))
            return _scalar(byname[name])
        end
        throw(ArgumentError("UV_DATA has no $base column"))
    end

    date, time, baseline = field(:DATE), field(:TIME), field(:BASELINE)
    source = optional_field(:SOURCE_ID, Symbol("SOURCE ID"), :SOURCE)
    freqid = optional_field(:FREQID)
    inttim = optional_field(:INTTIM)
    array = optional_field(:ARRAY)
    # No FITS-IDI table records scans — AIPS Memo 114r has no scan concept at
    # all — so they are normally inferred. A correlator that writes one anyway
    # knows better than any inference, and is believed.
    scan = optional_field(:SCAN, :SCAN_NUMBER, Symbol("SCAN NO"), :SCANNO)

    uu, vv, ww = uvw_field('U'), uvw_field('V'), uvw_field('W')
    bytes = _bytes(uv_hdu)
    col(::Type{T}, f) where {T} = _row_column(T, f, bytes, rowlen, n)
    col(::Type{T}, ::Nothing, default) where {T} = fill(T(default), n)
    col(::Type{T}, f, default) where {T} = col(T, f)

    uvw = Matrix{Float64}(undef, n, 3)
    for (k, f) in enumerate((uu, vv, ww))
        uvw[:, k] = col(Float64, f)
    end

    # A single array need not number itself, but two arrays sharing one
    # UV_DATA would collide in the 256*a1 + a2 baseline packing.
    if array !== nothing
        arrays = unique(col(Int32, array))
        length(arrays) == 1 || throw(ArgumentError(
            "UV_DATA spans arrays $(Int.(arrays)); this reader handles one array, " *
                "whose antenna numbering the BASELINE packing assumes"
        ))
    end

    return (;
        jd = col(Float64, date),
        day = col(Float64, time),
        baseline = col(Int32, baseline),
        source_id = col(Int32, source, 1),
        freqid = col(Int32, freqid, 1),
        integration_time = inttim === nothing ? nothing : col(Float64, inttim),
        scan = scan === nothing ? nothing : col(Int32, scan),
        uvw,
    )
end

# ── UV_DATA's data matrix ────────────────────────────────────────────────────

"""
    uv_matrix_layout(uv_hdu) -> (; ncomplex, no_stkd, no_chan, no_band)

The shape of one row's `FLUX` array, from the `MAXISn`/`CTYPEn` keywords.

FITS-IDI writes the regular axes of the data matrix fastest-varying first as
`COMPLEX`, `STOKES`, `FREQ`, `BAND`, optionally followed by degenerate `RA` and
`DEC` axes. `ncomplex` is 2 when a datum is a real and an imaginary part, and 3
when it carries a per-channel weight as well.
"""
function uv_matrix_layout(uv_hdu)
    maxis = Int(card(uv_hdu, "MAXIS", 0))
    maxis >= 4 || throw(ArgumentError(
        "UV_DATA data matrix has MAXIS = $maxis; the COMPLEX, STOKES, FREQ " *
            "and BAND axes are all required"
    ))
    for (n, want) in enumerate(("COMPLEX", "STOKES", "FREQ", "BAND"))
        got = card_string(uv_hdu, "CTYPE$n")
        got == want || throw(ArgumentError(
            "UV_DATA data-matrix axis $n is $(repr(got)), expected $(repr(want))"
        ))
    end
    axislength(n) = Int(card(uv_hdu, "MAXIS$n", 0))
    for n in 5:maxis
        axislength(n) == 1 || throw(ArgumentError(
            "UV_DATA data-matrix axis $n ($(card_string(uv_hdu, "CTYPE$n"))) has " *
                "length $(axislength(n)); only degenerate axes may follow BAND"
        ))
    end
    ncomplex = axislength(1)
    ncomplex in (2, 3) || throw(ArgumentError(
        "UV_DATA COMPLEX axis has length $ncomplex, expected 2 (real, imaginary) " *
            "or 3 (real, imaginary, weight)"
    ))
    return (; ncomplex, no_stkd = axislength(2), no_chan = axislength(3),
        no_band = axislength(4))
end

"""
    flux_element_type(uv_hdu) -> Type{<:AbstractFloat}

The floating-point type the `UV_DATA` `FLUX` column is stored in: `Float64` for
a `D`-format column, `Float32` for every other form the reader decodes,
including the integer ones, which carry fewer digits than a `Float32` holds.
"""
function flux_element_type(uv_hdu)
    fields = _descriptor(uv_hdu).fields
    i = findfirst(f -> _fieldname(f) == "FLUX", fields)
    i === nothing && throw(ArgumentError("UV_DATA has no FLUX column"))
    return fields[i].type === Float64 ? Float64 : Float32
end

"""
    read_uv_block(uv_hdu, rows; T = Float32) -> (; flux, weight, layout)

Read the data matrix and the `WEIGHT` column for a contiguous range of
`UV_DATA` rows.

`flux` is `(complex, polarization, channel, band, row)` and `weight` is
`(polarization, band, row)`, both in the file's own axis order and both
`Matrix{T}`-backed; `weight` is `nothing` when the `COMPLEX` axis is 3 and so
carries the weights itself. `layout` is the [`uv_matrix_layout`](@ref) the
shapes come from.

The rows are read as one contiguous block, which is nearly the whole cost of a
conversion, so ask for a whole partition at a time. A block is
`length(rows) * (2 to 3) * no_stkd * no_chan * no_band` values resident, which
for a file of any size means the caller must partition rather than ask for all
of them.
"""
function read_uv_block(uv_hdu, rows::AbstractUnitRange{<:Integer}; T::Type = Float32)
    desc = _descriptor(uv_hdu)
    n = _nrows(uv_hdu)
    (first(rows) >= 1 && last(rows) <= n) || throw(ArgumentError(
        "UV_DATA rows $rows lie outside the table's 1:$n"
    ))
    layout = uv_matrix_layout(uv_hdu)
    rowlen = _rowlength(uv_hdu)
    nrow = length(rows)

    byname = Dict(Symbol(rstrip(String(f.name))) => f for f in desc.fields)
    haskey(byname, :FLUX) || throw(ArgumentError("UV_DATA has no FLUX column"))
    flux_field = _unscaled(byname[:FLUX])
    nflux = layout.ncomplex * layout.no_stkd * layout.no_chan * layout.no_band
    flux_field.leng == nflux || throw(ArgumentError(
        "UV_DATA FLUX holds $(flux_field.leng) values per row; the data-matrix " *
            "axes give $nflux"
    ))

    bytes = _bytes(uv_hdu)
    start = rowlen * (first(rows) - 1)

    flux = reshape(
        _decode_field!(Matrix{T}(undef, nflux, nrow), flux_field, bytes, rowlen, start),
        layout.ncomplex, layout.no_stkd, layout.no_chan, layout.no_band, nrow
    )

    # AIPS Memo 114r ties the WEIGHT random parameter to the COMPLEX axis: with
    # MAXIS1 = 2 the weights have nowhere else to live and the parameter is
    # required, and with MAXIS1 = 3 they are in the data matrix and it must be
    # absent. Either mismatch means the weights are not what they appear.
    haveweight = haskey(byname, :WEIGHT)
    if haveweight == (layout.ncomplex == 3)
        throw(ArgumentError(
            haveweight ?
                "UV_DATA has a WEIGHT column and a COMPLEX axis of 3, which already " *
                "carries the weights; the standard permits only one of them" :
                "UV_DATA has no WEIGHT column and a COMPLEX axis of 2, which leaves " *
                "the weights unstated; the standard requires the column"
        ))
    end
    weight = if haveweight
        field = _unscaled(byname[:WEIGHT])
        nweight = layout.no_stkd * layout.no_band
        field.leng == nweight || throw(ArgumentError(
            "UV_DATA WEIGHT holds $(field.leng) values per row; NO_STKD x " *
                "NO_BAND gives $nweight"
        ))
        reshape(
            _decode_field!(Matrix{T}(undef, nweight, nrow), field, bytes, rowlen, start),
            layout.no_stkd, layout.no_band, nrow
        )
    end

    return (; flux, weight, layout)
end
