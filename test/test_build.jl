using XRadio: subdataset

@testset "building datasets from the schema" begin
    stations = ["BR", "FD"]
    positions = [-2.11e6 -1.32e6; -3.70e6 -4.82e6; 4.55e6 3.55e6]

    @testset "a dataset the standard recognizes" begin
        ds = subdataset(
            :weather;
            station_name = stations,
            cartesian_pos_label = ["x", "y", "z"],
            time_weather = [1.77e9, 1.77e9 + 300],
            station_position = positions,
            temperature = [283.15 284.15; 283.65 284.65],
        )
        @test isempty(check(ds))
        @test metadata(ds)[:type] == "weather"

        # An index coordinate is the dimension it is named after; everything
        # else is a layer over those.
        @test lookup(ds[:temperature], StationName) == stations
        @test collect(keys(ds)) == [:station_position, :temperature]

        # Arrays are in Julia order, the reverse of the order the standard
        # states its dimensions in.
        @test map(XRadio.msv4name, dims(ds[:temperature])) == (:time_weather, :station_name)
        @test map(XRadio.msv4name, dims(ds[:station_position])) ==
            (:cartesian_pos_label, :station_name)

        # The attributes the standard fixes come from the schema, not from the
        # caller.
        @test metadata(ds[:temperature])[:units] == "K"
        @test metadata(ds[:temperature])[:type] == "quantity"
        @test metadata(ds[:station_position])[:frame] == "ITRS"
        @test metadata(ds[:station_position])[:coordinate_system] == "geocentric"
        tmeta = metadata(lookup(ds[:temperature], TimeWeather))
        @test (tmeta[:type], tmeta[:units], tmeta[:scale], tmeta[:format]) ==
            ("time_weather", "s", "utc", "unix")
    end

    @testset "a variable takes the shape it was given the dimensions for" begin
        # MSv4 lets weather be sampled on the correlated data's time axis or on
        # the weather station's own, and which one a dataset uses is which one
        # it carries.
        ds = subdataset(
            :weather;
            station_name = stations,
            cartesian_pos_label = ["x", "y", "z"],
            time = [1.77e9],
            station_position = positions,
            pressure = [90000.0 90100.0],
        )
        @test map(XRadio.msv4name, dims(ds[:pressure])) == (:time, :station_name)
        @test isempty(check(ds))

        # Six shapes are admitted for a system temperature, and giving the
        # coarse calibration axes rather than the fine ones selects among them.
        sc = subdataset(
            :system_calibration;
            antenna_name = stations, station_name = stations,
            receptor_label = ["pol_0", "pol_1"],
            polarization_type = ["R" "R"; "L" "L"],
            time_system_cal = [1.77e9],
            tsys = reshape([30.0, 31.0, 40.0, 41.0], 2, 1, 2),
        )
        @test map(XRadio.msv4name, dims(sc[:tsys])) ==
            (:receptor_label, :time_system_cal, :antenna_name)
        @test isempty(check(sc))
    end

    @testset "what is not given is not written" begin
        ds = subdataset(
            :weather;
            station_name = stations,
            cartesian_pos_label = ["x", "y", "z"],
            time_weather = [1.77e9],
            station_position = positions,
            temperature = [283.15 284.15],
            pressure = nothing,
        )
        # An optional variable a file does not carry stays absent rather than
        # being invented as zeros.
        @test !haskey(ds, :pressure)
        @test haskey(ds, :temperature)
    end

    @testset "attributes the standard leaves open" begin
        direction = deg2rad.([187.7; 12.4;;])
        ds = subdataset(
            :field_and_source;
            field_name = ["3C273"], sky_dir_label = ["ra", "dec"],
            source_name = ["3C273"],
            field_phase_center_direction = direction,
            array_metadata = (; field_phase_center_direction = (; frame = "fk5")),
        )
        @test metadata(ds[:field_phase_center_direction])[:frame] == "fk5"
        # The schema's own default stands where the caller says nothing.
        @test metadata(ds[:field_phase_center_direction])[:units] == "rad"

        # Dataset attributes the standard requires but does not fix.
        gc = subdataset(
            :gain_curve;
            antenna_name = stations, station_name = stations,
            mount = ["ALT-AZ", "ALT-AZ"], telescope_name = ["VLBA", "VLBA"],
            receptor_label = ["pol_0", "pol_1"],
            polarization_type = ["R" "R"; "L" "L"],
            poly_term = Int64[0, 1],
            gain_curve = zeros(2, 2, 2),
            gain_curve_interval = [3600.0, 3600.0],
            gain_curve_sensitivity = fill(0.093, 2, 2),
            metadata = (; measured_date = "2026-01-01T00:00:00.000"),
        )
        @test metadata(gc)[:measured_date] == "2026-01-01T00:00:00.000"
        @test isempty(check(gc))
    end

    @testset "what it refuses" begin
        # A name the dataset does not define is a typo, and saying which names
        # it does define is more use than a `KeyError`.
        @test_throws "is not part of the `weather` dataset" subdataset(
            :weather; station_name = stations, temperture = [280.0 281.0]
        )
        # A variable whose dimensions were not supplied cannot be placed.
        @test_throws "and the dimensions given are" subdataset(
            :weather; station_name = stations, temperature = [280.0 281.0]
        )
        @test_throws "is not an MSv4 schema" subdataset(:barometer)
    end
end

@testset "building a Measurement Set" begin
    times = [1.7e9, 1.7e9 + 10]
    freqs = [2.3e11, 2.3e11 + 1.0e6, 2.3e11 + 2.0e6]
    shape = (2, 3, 1, 2)
    fs = XRadio.Testing.field_and_source()
    build(; kw...) = XRadio.measurement_set(;
        time = times, frequency = freqs, polarization = ["RR", "LL"],
        baseline_id = [0], baseline_antenna1_name = ["A1"], baseline_antenna2_name = ["A2"],
        visibility = ones(ComplexF32, shape), flag = falses(shape), weight = ones(Float32, shape),
        field_name = fill("field_0", 2), scan_name = fill("1", 2),
        kw...,
    )

    @testset "the base data group names what was given" begin
        ms = build(; subdatasets = (; field_and_source_base = fs))
        g = datagroup(ms)
        @test g[:correlated_data] == "VISIBILITY"
        @test (g[:flag], g[:weight]) == ("FLAG", "WEIGHT")
        @test !haskey(g, :uvw)
        @test g[:field_and_source] == "field_and_source_base_xds"
        @test g[:description] == ""
        @test metadata(ms)[:type] == "visibility"
        @test metadata(ms)[:creator][:software_name] == "XRadio.jl"
        @test field_and_source(ms) == ms.field_and_source_base
        # Linked, as assignment links: it reads `time` from the Measurement Set.
        @test getfield(ms.field_and_source_base, :tree) !== nothing
    end

    @testset "arrays land in the standard's layout with its attributes" begin
        ms = build()
        @test dims(ms[:visibility]) == dims(ms, (Polarization, Frequency, BaselineID, Ti))
        @test metadata(ms[:visibility])[:units] == "Jy"
        @test lookup(ms, Frequency) == freqs
    end

    @testset "metadata adds to and replaces the defaults" begin
        ms = build(; metadata = (; creation_date = "2020-01-01T00:00:00.000", note = 1))
        @test metadata(ms)[:creation_date] == "2020-01-01T00:00:00.000"
        @test metadata(ms)[:note] == 1
        @test haskey(metadata(ms), :schema_version)
    end

    @testset "what cannot be built is an error" begin
        @test_throws "`visibilities` is not part of the `visibility` dataset" build(; visibilities = 1)
        @test_throws "is `:visibility` or `:spectrum`" XRadio.measurement_set(:antenna)
    end
end

@testset "test data" begin
    Testing = XRadio.Testing

    @testset "one Measurement Set per scan and spectral window" begin
        ps = Testing.processing_set(; nantenna = 4, nscan = 3, nspw = 2, ntime = 5, nchan = 6)
        @test isempty(check(ps))
        @test collect(keys(ps)) == [Symbol("synthetic_0", i) for i in 0:5]
        @test [scans(ms) for ms in values(ps)] == [["1"], ["1"], ["2"], ["2"], ["3"], ["3"]]
        @test [spectralwindow(ms) for ms in values(ps)] == repeat(["spw_0", "spw_1"], 3)
        ms = ps[3]
        @test size(ms[:visibility]) == (4, 6, 6, 5)
        @test all(==(1), ms[:visibility]) && !any(ms[:flag]) && all(==(1), ms[:weight])
        @test first(lookup(ms, Ti)) - first(lookup(ps[1], Ti)) == 7 * 30.0
        @test first(lookup(ps[2], Frequency)) - first(lookup(ps[1], Frequency)) == 1.0e8
        @test antennas(ms) == ["A1", "A2", "A3", "A4"]
        @test Dates.unix2datetime(first(lookup(ps[1], Ti))) == Dates.DateTime(2021, 3, 4)
    end

    @testset "an antenna on no baseline" begin
        ms = Testing.processing_set(; nantenna = 4, drop = ["A2"])[1]
        @test antennas(ms) == ["A1", "A2", "A3", "A4"]
        @test length(baselines(ms)) == 3
        @test "A2" ∉ vcat(collect(ms[:baseline_antenna1_name]), collect(ms[:baseline_antenna2_name]))
    end

    @testset "chosen axes and element type" begin
        ps = Testing.processing_set(;
            antennas = ["X", "Y"], polarizations = ["XX", "YY"], eltype = ComplexF64,
            spw_starts = [8.0e9], scan_gap = 1000.0, nscan = 2,
        )
        @test eltype(ps[1][:visibility]) == ComplexF64
        @test eltype(ps[1][:weight]) == Float64
        @test collect(lookup(ps[1], Polarization)) == ["XX", "YY"]
        @test first(lookup(ps[1], Frequency)) == 8.0e9
        @test first(lookup(ps[2], Ti)) - first(lookup(ps[1], Ti)) == 1000.0
        @test_throws "2 spectral window starts for nspw = 1" Testing.processing_set(; spw_starts = [1.0e9, 2.0e9])
    end

    @testset "values are set afterwards" begin
        ps = transform(Testing.processing_set(; nscan = 2)) do ms
            out = copy(ms)
            out[:visibility] = rebuild(ms[:visibility], fill(2.0f0 + 0im, size(ms[:visibility])))
            out
        end
        @test all(==(2), ps[2][:visibility])
        p = joinpath(mktempdir(), "synthetic.ps.zarr")
        write(p, ps)
        @test read(ProcessingSet, p)[2][:visibility] == ps[2][:visibility]
    end
end
