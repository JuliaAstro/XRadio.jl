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
