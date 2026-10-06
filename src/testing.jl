# Conforming MSv4 data for tests and examples.
#
# Every generator here builds conforming data with values that say nothing,
# as `xradio.testing` does for the reference implementation:
# visibilities of one, unit weights, no flags, zero UVW. A test sets the values
# it needs afterwards, with `transform` and `rebuild` as for any modified
# Measurement Set.
#
# ```julia
# using XRadio: Testing
# ps = Testing.processing_set(; nantenna = 4, nspw = 2, nscan = 3)
# ```
module Testing

using Dates: Dates, DateTime, datetime2unix
using OrderedCollections: OrderedDict
using ..XRadio: XRadio, ProcessingSet, MeasurementSet, Measure, MountAltAz,
    subdataset, set_mounts!

_quantity(v, units) = Measure(v, Dict{Symbol, Any}(:type => "quantity", :units => units))

"""
    Testing.antenna(names; positions, mounts, receptors, diameter = 25.0) -> Dataset

An antenna dataset for the antennas `names`.

`positions` is a 3 × n matrix of geocentric positions in metres; by default
antenna `i` sits at `(1, 2, 3) · 10⁴ i` m, so no two are within a site of one
another. `mounts` defaults to alt-azimuth, and `receptors`, every antenna's
receptor list of any length, to `("R", "L")`.
"""
function antenna(
        names::AbstractVector{<:AbstractString};
        positions = [1.0e4 * c * i for c in 1:3, i in eachindex(names)],
        mounts = fill(MountAltAz(), length(names)),
        receptors = ("R", "L"),
        diameter = 25.0,
    )
    n = length(names)
    xds = subdataset(
        :antenna;
        antenna_name = collect(String, names),
        cartesian_pos_label = ["x", "y", "z"],
        receptor_label = ["pol_$(r - 1)" for r in eachindex(receptors)],
        antenna_position = Matrix{Float64}(positions),
        station_name = collect(String, names),
        telescope_name = fill("SYNTHETIC", n),
        polarization_type = [String(r) for r in receptors, _ in 1:n],
        antenna_receptor_angle = zeros(length(receptors), n),
        antenna_dish_diameter = fill(Float64(diameter), n),
        metadata = (; overall_telescope_name = "SYNTHETIC", relocatable_antennas = false),
        array_metadata = (;
            antenna_position = (;
                frame = "ITRS", coordinate_system = "geocentric", origin_object_name = "earth",
            ),
        ),
    )
    return set_mounts!(xds, mounts)
end

"""
    Testing.field_and_source(; source = "SRC", field = "field_0",
                               direction = (0.5, 0.3)) -> Dataset

A field-and-source dataset for one field observing one source, whose phase
centre is `direction`, (right ascension, declination) in radians, ICRS.
"""
function field_and_source(;
        source::AbstractString = "SRC", field::AbstractString = "field_0",
        direction = (0.5, 0.3),
    )
    return subdataset(
        :field_and_source;
        field_name = [String(field)],
        sky_dir_label = ["ra", "dec"],
        source_name = [String(source)],
        field_phase_center_direction = reshape(collect(Float64, direction), 2, 1),
        array_metadata = (; field_phase_center_direction = (; frame = "icrs")),
    )
end

"""
    Testing.measurement_set(; antennas, times, frequencies, polarizations,
                              kw...) -> MeasurementSet

One Measurement Set: every baseline among `antennas` but those involving an
antenna in `drop`, sampled at `times` (Unix seconds) and `frequencies` (Hz, one
spectral window).

# Keywords
- `antennas = ["A1", "A2", "A3"]`, `drop = String[]`: antennas in the antenna
  dataset, and those among them on no baseline.
- `times`, `frequencies`, `polarizations = ["RR", "RL", "LR", "LL"]`.
- `integration_time`: defaults to the spacing of `times`, or 1 s for one time.
- `channel_width`: defaults to the spacing of `frequencies`, or 1 MHz.
- `reference_frequency = first(frequencies)`, `spectral_window = "spw_0"`.
- `scan = "1"`, `field = "field_0"`, `source = "SRC"`, `direction`.
- `eltype = ComplexF32`: the visibilities' element type; weights take its real
  type.
- `antenna_xds`: the antenna dataset, by default
  [`Testing.antenna`](@ref)`(antennas)`.
"""
function measurement_set(;
        antennas::AbstractVector{<:AbstractString} = ["A1", "A2", "A3"],
        drop = String[],
        times = [1.6148160e9, 1.6148160e9 + 30],
        frequencies = 2.3e11 .+ 2.0e6 .* (0:7),
        polarizations = ["RR", "RL", "LR", "LL"],
        integration_time = length(times) > 1 ? times[2] - times[1] : 1.0,
        channel_width = length(frequencies) > 1 ? frequencies[2] - frequencies[1] : 1.0e6,
        reference_frequency = first(frequencies),
        spectral_window::AbstractString = "spw_0",
        scan::AbstractString = "1",
        field::AbstractString = "field_0",
        source::AbstractString = "SRC",
        direction = (0.5, 0.3),
        eltype::Type{<:Complex} = ComplexF32,
        antenna_xds = antenna(antennas),
    )
    kept = [a for a in antennas if a ∉ drop]
    pairs_ = [(a, b) for (i, a) in enumerate(kept) for b in kept[(i + 1):end]]
    nb, nt, nf, np = length(pairs_), length(times), length(frequencies), length(polarizations)
    shape = (np, nf, nb, nt)
    freq_attrs = (;
        spectral_window_name = String(spectral_window),
        spectral_window_intents = ["UNSPECIFIED"],
        reference_frequency = Measure(
            Float64(reference_frequency),
            Dict{Symbol, Any}(:type => "spectral_coord", :units => "Hz", :observer => "icrs"),
        ),
        channel_width = _quantity(Float64(channel_width), "Hz"),
    )
    return XRadio.measurement_set(;
        time = collect(Float64, times),
        frequency = collect(Float64, frequencies),
        polarization = collect(String, polarizations),
        baseline_id = collect(0:(nb - 1)),
        baseline_antenna1_name = first.(pairs_),
        baseline_antenna2_name = last.(pairs_),
        uvw_label = ["u", "v", "w"],
        visibility = ones(eltype, shape),
        flag = falses(shape),
        weight = ones(real(eltype), shape),
        uvw = zeros(3, nb, nt),
        field_name = fill(String(field), nt),
        scan_name = fill(String(scan), nt),
        metadata = (;
            observation_info = Dict{Symbol, Any}(
                :observer => ["synthetic"], :release_date => "2000-01-01T00:00:00.000",
                :project_UID => "synthetic",
            ),
            processor_info = Dict{Symbol, Any}(:type => "CORRELATOR", :sub_type => "SYNTHETIC"),
        ),
        array_metadata = (;
            time = (; integration_time = _quantity(Float64(integration_time), "s")),
            frequency = freq_attrs,
            uvw = (; frame = "icrs"),
            scan_name = (; scan_intents = ["OBSERVE_TARGET#ON_SOURCE"]),
        ),
        subdatasets = (;
            antenna = antenna_xds,
            field_and_source_base = field_and_source(; source, field, direction),
        ),
    )
end

"""
    Testing.processing_set(; kw...) -> ProcessingSet

A processing set of `nscan` scans observed in `nspw` spectral windows, one
Measurement Set per scan and window, named `name_00`, `name_01` and so on,
scan by scan.

# Keywords
- `nantenna = 3` antennas named `A1`, `A2`, …, or `antennas` to name them;
  `drop`, `polarizations`, `eltype`: as for [`Testing.measurement_set`](@ref).
- `ntime = 4` samples per scan, `integration_time = 30.0` s apart, starting at
  `start = DateTime(2021, 3, 4)`; scans begin `scan_gap` seconds apart, by
  default two samples after the previous one ends.
- `nspw = 1` windows of `nchan = 8` channels `channel_width = 2e6` Hz wide,
  starting at `spw_starts`, by default `reference_frequency = 2.3e11` Hz and
  every `spw_separation = 1e8` Hz above it.
- `name = "synthetic"`.
"""
function processing_set(;
        nantenna::Integer = 3,
        antennas::AbstractVector{<:AbstractString} = ["A$i" for i in 1:nantenna],
        nscan::Integer = 1, ntime::Integer = 4, integration_time::Real = 30.0,
        start::DateTime = DateTime(2021, 3, 4),
        scan_gap::Real = (ntime + 2) * integration_time,
        nspw::Integer = 1, nchan::Integer = 8, channel_width::Real = 2.0e6,
        reference_frequency::Real = 2.3e11, spw_separation::Real = 1.0e8,
        spw_starts = [reference_frequency + (w - 1) * spw_separation for w in 1:nspw],
        name::AbstractString = "synthetic",
        kw...,
    )
    length(spw_starts) == nspw || throw(ArgumentError(
        "$(length(spw_starts)) spectral window starts for nspw = $nspw"
    ))
    t0 = datetime2unix(start)
    ant = antenna(antennas)
    width = max(2, ndigits(nscan * nspw - 1))
    sets = OrderedDict{Symbol, MeasurementSet}()
    for s in 1:nscan, w in 1:nspw
        sets[Symbol(name, "_", lpad(length(sets), width, '0'))] = measurement_set(;
            antennas, antenna_xds = ant,
            times = t0 + (s - 1) * scan_gap .+ integration_time .* (0:(ntime - 1)),
            frequencies = spw_starts[w] .+ channel_width .* (0:(nchan - 1)),
            integration_time, channel_width,
            reference_frequency = spw_starts[w],
            spectral_window = "spw_$(w - 1)",
            scan = string(s),
            kw...,
        )
    end
    return ProcessingSet(sets)
end

end
