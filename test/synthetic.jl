# A small, fully conforming Measurement Set built from nothing.
#
# Round-tripping a downloaded store only shows that whatever it contains
# survives the trip. Building one here fixes what conformance looks like
# independently, so the writer is tested against the standard rather than
# against a fixture that may itself be at fault.

using Dates
using DimensionalData
using DimensionalData.Lookups: Categorical, Sampled

# A dimension type of a caller's own, for the cases that turn on a dimension
# having no MSv4 name of its own: a store spells such a dimension `Dim{name}`,
# and this is what a reader may `set` one to after the fact.
DimensionalData.@dim SolutionInterval "Calibration solution interval"

"A `Measure` in the form the standard's `dataarray` attributes take."
measure(v, units, type = "quantity") =
    XRadio.Measure(v, Dict{Symbol, Any}(:units => units, :type => type))

"""
    synthetic_visibility(; ntime, nbaseline, nfreq, T, C)

A conforming interferometric Measurement Set of the given shape.

`T` is the element type of the real-valued arrays and `C` that of the
visibilities, so the same builder covers the standard's `float32`/`float64` and
`complex64`/`complex128` alternatives — and, under a non-strict check, anything
else of the same kind.
"""
function synthetic_visibility(;
        ntime::Int = 4, nbaseline::Int = 3, nfreq::Int = 5,
        T::Type = Float64, C::Type = ComplexF32,
    )
    antnames = ["AA", "BB", "CC"]
    npol = 2
    return XRadio.measurement_set(;
        time = collect(range(1.7e9; step = 10.0, length = ntime)),
        frequency = collect(range(2.3e11; step = 1.0e6, length = nfreq)),
        polarization = ["XX", "YY"],
        baseline_id = collect(1:nbaseline),
        baseline_antenna1_name = ["AA", "AA", "BB"][1:nbaseline],
        baseline_antenna2_name = ["BB", "CC", "CC"][1:nbaseline],
        uvw_label = ["u", "v", "w"],
        # Julia order: the reverse of the (time, baseline_id, frequency,
        # polarization) the standard states.
        visibility = rand(C, npol, nfreq, nbaseline, ntime),
        flag = falses(npol, nfreq, nbaseline, ntime),
        weight = ones(T, npol, nfreq, nbaseline, ntime),
        uvw = rand(T, 3, nbaseline, ntime),
        field_name = fill("field_0", ntime),
        scan_name = fill("1", ntime),
        metadata = (;
            creation_date = "2026-01-01T00:00:00.000",
            creator = Dict{Symbol, Any}(:software_name => "XRadio.jl", :version => "0.1.0"),
            observation_info = Dict{Symbol, Any}(
                :observer => ["Nobody"],
                :release_date => "2027-01-01T00:00:00.000",
                :project_UID => "uid://X0/X1/X2",
            ),
            processor_info = Dict{Symbol, Any}(:type => "CORRELATOR", :sub_type => "TEST"),
            data_groups = Dict{Symbol, Any}(
                :base => Dict{Symbol, Any}(
                    :correlated_data => "VISIBILITY", :flag => "FLAG",
                    :weight => "WEIGHT", :uvw => "UVW",
                    :field_and_source => "field_and_source_base_xds",
                    :description => "Correlator output",
                    :date => "2026-01-01T00:00:00.000",
                ),
            ),
        ),
        array_metadata = (;
            time = (;
                scale = "utc", format = "unix", integration_time = measure(10.0, "s"),
            ),
            frequency = (;
                observer = "icrs",
                spectral_window_name = "spw_0", spectral_window_intents = ["CONTINUUM"],
                reference_frequency = XRadio.Measure(
                    2.3e11,
                    Dict{Symbol, Any}(
                        :units => "Hz", :type => "spectral_coord", :observer => "icrs"
                    )
                ),
                channel_width = measure(1.0e6, "Hz"),
            ),
            uvw = (; frame = "icrs"),
            scan_name = (; scan_intents = ["OBSERVE_TARGET#ON_SOURCE"]),
        ),
        subdatasets = (;
            antenna = synthetic_antenna(antnames),
            field_and_source_base = synthetic_field_and_source(),
        ),
    )
end

function synthetic_antenna(antnames)
    ant = AntennaName(Categorical(antnames))
    receptor = ReceptorLabel(Categorical(["X", "Y"]))
    cart = CartesianPosLabel(Categorical(["x", "y", "z"]))
    n = length(antnames)

    ds = Dataset(;
        metadata = Dict{Symbol, Any}(
            :type => "antenna",
            :overall_telescope_name => "TEST",
            :relocatable_antennas => false,
        )
    )
    ds[:antenna_position] = DimArray(
        rand(3, n) .* 6.4e6, (cart, ant);
        metadata = Dict{Symbol, Any}(
            :type => "location", :units => "m", :frame => "ITRS",
            :coordinate_system => "geocentric", :origin_object_name => "earth",
        )
    )
    ds[:station_name] = DimArray(antnames, (ant,))
    ds[:mount] = DimArray(fill("ALT-AZ", n), (ant,))
    ds[:telescope_name] = DimArray(fill("TEST", n), (ant,))
    # (antenna_name, receptor_label) in MSv4 order, so reversed here.
    ds[:polarization_type] = DimArray(repeat(["X", "Y"], 1, n), (receptor, ant))
    return ds
end

function synthetic_field_and_source()
    field = FieldName(Categorical(["field_0"]))
    skylab = SkyDirLabel(Categorical(["ra", "dec"]))
    ds = Dataset(; metadata = Dict{Symbol, Any}(:type => "field_and_source"))
    ds[:source_name] = DimArray(["source_0"], (field,))
    ds[:field_phase_center_direction] = DimArray(
        [0.5; 0.3;;], (skylab, field);
        metadata = Dict{Symbol, Any}(:type => "sky_coord", :units => "rad", :frame => "icrs")
    )
    return ds
end

"A processing set holding `n` copies of the synthetic Measurement Set."
function synthetic_set(n::Int = 2; kw...)
    sets = OrderedCollections.OrderedDict{Symbol, MeasurementSet}(
        Symbol("test_", lpad(i - 1, 2, '0')) => synthetic_visibility(; kw...) for i in 1:n
    )
    return ProcessingSet(sets)
end
