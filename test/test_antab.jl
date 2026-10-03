@testset "ANTAB" begin
    dir = mktempdir()
    antab(text) = (path = joinpath(dir, "$(hash(text)).antab"); write(path, text); path)

    @testset "parser: GAIN + TSYS layouts" begin
        parsed = XRadio._parse_antab(antab("""
        GAIN AA ELEV DPFU = 0.031000 POLY = 1.0 /
        GAIN MG ELEV DPFU = 0.0179, 0.0168 POLY = 0.727119, 0.00947339, -0.00008222 /
        GAIN NN ALTAZ FREQ = 220000, 225000 DPFU = 1.0 POLY = 1.0, 0.5 /

        TSYS AA  FT=1.0  TIMEOFF=0
        INDEX = 'L1|R1', 'L2|R2'
        /
        100 12:00:00     2.5    3.5
        100 12:30.5      2.7    3.7
        /
        TSYS MG timeoff= 0.0  FT = 1.0  INDEX = 'R1:32', 'L1:32:2', 'X' /
        100 12:00:00     200.0  220.0  0.0  ! a comment
        /
        """))
        aa_gain = only(parsed.gains["AA"])
        @test (aa_gain.kind, aa_gain.dpfu, aa_gain.poly, aa_gain.freq) == ("POWER(EL)", [0.031], [1.0], nothing)
        @test only(parsed.gains["MG"]).dpfu == [0.0179, 0.0168]
        @test only(parsed.gains["MG"]).poly == [0.727119, 0.00947339, -0.00008222]
        @test only(parsed.gains["NN"]).kind == "POWER(ZA)"
        @test only(parsed.gains["NN"]).freq == (2.2e11, 2.25e11)

        # A header over several lines, paired receptors, and decimal minutes.
        aa = parsed.tsys["AA"]
        @test aa.columns == [[('L', 1:1:1), ('R', 1:1:1)], [('L', 2:1:2), ('R', 2:1:2)]]
        @test aa.doy == [100, 100]
        @test aa.seconds == [43200.0, 45030.0]
        @test aa.values == [2.5 3.5; 2.7 3.7]

        # A one-line header, ranges with and without a step, and an unused column.
        mg = parsed.tsys["MG"]
        @test mg.columns == [[('R', 1:1:32)], [('L', 1:2:31)], Tuple{Char, StepRange{Int, Int}}[]]
        @test mg.values == [200.0 220.0 0.0]

        @test_throws "only ELEV and ALTAZ" XRadio._parse_antab(antab("GAIN AA EQUAT DPFU = 1.0 POLY = 1.0 /\n"))
        @test_throws "no DPFU" XRadio._parse_antab(antab("GAIN AA ELEV POLY = 1.0 /\n"))
        @test_throws "no GAIN record" XRadio._parse_antab(antab("TSYS AA INDEX = 'R1' /\n100 12:00:00 1.0\n/\n"))
        @test_throws "uses INDEX2" XRadio._parse_antab(antab(
            "GAIN AA ELEV DPFU = 1.0 POLY = 1.0 /\nTSYS AA INDEX = 'R1' INDEX2 = 'L1' /\n/\n"
        ))
        @test_throws "names 2 columns, and this row has 1" XRadio._parse_antab(antab(
            "GAIN AA ELEV DPFU = 1.0 POLY = 1.0 /\nTSYS AA INDEX = 'R1', 'L1' /\n100 12:00:00 1.0\n/\n"
        ))
        @test_throws "unrecognized INDEX entry" XRadio._parse_antab(antab(
            "GAIN AA ELEV DPFU = 1.0 POLY = 1.0 /\nTSYS AA INDEX = 'Q1' /\n/\n"
        ))
    end

    # `Testing.processing_set` starts at 2021-03-04T00:00:00, day 63, with four
    # 30 s integrations, in windows at 230 and 230.1 GHz.
    text = """
    GAIN AA ELEV DPFU = 0.031 POLY = 1.0 /
    GAIN MG ELEV DPFU = 0.0179, 0.0168 POLY = 0.7, 0.01 /
    GAIN NN ALTAZ FREQ = 220000, 230050 DPFU = 1.0 POLY = 1.0, 0.5 /
    GAIN NN ALTAZ FREQ = 230050, 240000 DPFU = 2.0 POLY = 1.0, 0.25 /
    GAIN QQ ELEV DPFU = 9.0 POLY = 1.0 /

    TSYS AA INDEX = 'L1|R1', 'L2|R2' /
    63 00:00:00     2.5    3.5
    63 00:00:30     2.7    3.7
    63 01:00:00     9.9    9.9
    /
    TSYS MG INDEX = 'R1:2', 'L1:2' /
    62 23:59:59     200.0  220.0
    63 00:01.0      180.0  210.0
    /
    """
    built() = XRadio.Testing.processing_set(; antennas = ["AA", "MG", "NN", "ZZ"], nspw = 2)
    ifs = Dict("spw_0" => 1, "spw_1" => 2)

    @testset "into a built processing set" begin
        ps = built()
        @test (@test_logs (:warn, r"does not say when the gain curves") read_antab!(ps, antab(text); ifs)) === ps

        for (name, k) in (:synthetic_00 => 1, :synthetic_01 => 2)
            ms = ps[name]
            @test isempty(XRadio.check(ms))
            tsys = system_temperatures(ms)
            # Rows inside the span, on the file's own clock; the row a second
            # before midnight is the previous day of the same year, and the
            # row an hour later is outside the span.
            @test lookup(tsys, TimeSystemCal) == 1.614816e9 .+ [-1.0, 0.0, 30.0, 60.0]
            at(a, r) = collect(tsys[AntennaName(At(a)), ReceptorLabel(At(r))])
            per_if = k == 1 ? [2.5, 2.7] : [3.5, 3.7]
            @test isequal(at("AA", "pol_0"), [NaN; per_if; NaN])
            @test isequal(at("AA", "pol_1"), [NaN; per_if; NaN])
            @test isequal(at("MG", "pol_0"), [200.0, NaN, NaN, 180.0])
            @test isequal(at("MG", "pol_1"), [220.0, NaN, NaN, 210.0])
            @test all(isnan, at("NN", "pol_0")) && all(isnan, at("ZZ", "pol_1"))
        end

        # DPFU is right then left; a GAIN record holds where its FREQ range
        # contains the window's reference frequency.
        low, high = gaincurves(ps[:synthetic_00]), gaincurves(ps[:synthetic_01])
        @test collect(low[:sensitivity][AntennaName(At("MG"))]) == [0.0179, 0.0168]
        @test collect(low[:sensitivity][AntennaName(At("AA"))]) == [0.031, 0.031]
        @test collect(low[:sensitivity][AntennaName(At("NN"))]) == [1.0, 1.0]
        @test collect(high[:sensitivity][AntennaName(At("NN"))]) == [2.0, 2.0]
        @test all(isnan, high[:sensitivity][AntennaName(At("ZZ"))])
        @test gain(low[:curve][AntennaName(At("MG")), ReceptorLabel(At("pol_0"))], deg2rad(40)) ≈ 0.7 + 0.01 * 40
        @test XRadio.gain_variable(high[:curve][AntennaName(At("NN")), ReceptorLabel(At("pol_0"))]) === :zenith_angle
        @test gain(high[:curve][AntennaName(At("NN")), ReceptorLabel(At("pol_0"))], deg2rad(40)) ≈ 1.0 + 0.25 * 40
    end

    @testset "what it refuses" begin
        path = antab(text)
        ps = built()
        other = built()
        @test_logs (:warn,) read_antab!(other, path; ifs)
        ps[:synthetic_01][:gain_curve] = branches(other[:synthetic_01])[:gain_curve]
        @test_throws "`spw_1` already has a `gain_curve` dataset" read_antab!(ps, path; ifs)
        @test !haskey(branches(ps[:synthetic_00]), :system_calibration)

        @test_throws "does not record its FITS-IDI band" read_antab!(built(), path)
        @test_throws "gives no IF for spectral window `spw_1`" read_antab!(built(), path; ifs = Dict("spw_0" => 1))
        @test_throws "covers IFs [1, 2] and the processing set's spectral windows are IFs [1, 3]" read_antab!(
            built(), path; ifs = Dict("spw_0" => 1, "spw_1" => 3)
        )
        @test_throws "no antenna of the processing set" read_antab!(
            XRadio.Testing.processing_set(; antennas = ["X1", "X2"], nspw = 2), path; ifs
        )
        linear = antab("GAIN AA ELEV DPFU = 1.0 POLY = 1.0 /\nTSYS AA INDEX = 'X1:2' /\n63 00:00:00 1.0\n/\n")
        @test_throws "gives polarization X, and the antenna's receptors are R, L" read_antab!(built(), linear; ifs)
    end

    @testset "into a FITS-IDI store, by band" begin
        # Bands out of frequency order: band 1 is the higher. The store records
        # each window's band, and the ANTAB IFs follow it.
        gt = synthetic_fitsidi(joinpath(dir, "antab.idifits"); band_freq = [1.28e8, 0.0])
        ps = open(ProcessingSet, fitsidi2msv4(gt.path, joinpath(dir, "antab.ps.zarr");
            release_date = "2030-01-01T00:00:00.000"))
        band(ms) = metadata(lookup(dims(ms, XRadio.Frequency)))[:fitsidi_band]
        @test [band(ps[k]) for k in (:antab_00, :antab_01)] == [1, 2]
        @test XRadio._reference_frequency(ps[:antab_00]) > XRadio._reference_frequency(ps[:antab_01])

        # The first scan starts at 2026-01-01T00:00:00 with 1 s integrations,
        # so the row 0.2 s before is the last day of 2025.
        path = antab("""
        GAIN AA ELEV DPFU = 0.5 POLY = 1.0 /
        GAIN BB ELEV DPFU = 0.25 POLY = 1.0 /
        TSYS AA INDEX = 'R1', 'L1', 'R2', 'L2' /
        365 23:59:59.8   11.0  12.0  21.0  22.0
        1 00:00:01       13.0  14.0  23.0  24.0
        1 00:10:00       15.0  16.0  25.0  26.0
        /
        TSYS BB INDEX = 'R1:2', 'L1:2' /
        1 00:10:00.5     30.0  40.0
        /
        """)
        @test_logs (:warn,) read_antab!(ps, path)
        @test all(isempty ∘ XRadio.check, values(ps))
        first_scan = system_temperatures(ps[:antab_00])
        @test lookup(first_scan, TimeSystemCal) ≈ datetime2unix(DateTime(2026)) .+ [-0.2, 1.0]
        @test collect(first_scan[AntennaName(At("AA")), ReceptorLabel(At("pol_0"))]) == [11.0, 13.0]
        @test collect(system_temperatures(ps[:antab_01])[AntennaName(At("AA")), ReceptorLabel(At("pol_1"))]) == [22.0, 24.0]
        second_scan = system_temperatures(ps[:antab_02])
        @test isequal(collect(second_scan[AntennaName(At("BB")), ReceptorLabel(At("pol_1"))]), [NaN, 40.0])
        @test all(isnan, second_scan[AntennaName(At("CC"))])
        @test collect(gaincurves(ps[:antab_01])[:sensitivity][AntennaName(At("BB"))]) == [0.25, 0.25]
    end
end
