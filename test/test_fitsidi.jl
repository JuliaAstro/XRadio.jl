# The FITS-IDI reader is internal to the extension; `using FITSFiles` in
# runtests.jl is what loads it.
const IDI = Base.get_extension(XRadio, :XRadioFITSFilesExt)
using .IDI: baseline_antennas, channel_frequencies, pol_axis, pol_labels,
    read_uv_rows, scan_ranges

@testset "FITS-IDI reading" begin
    @testset "polarization axis" begin
        # Products are named by receptor index, not by polarization: the
        # circular and linear code blocks describe the same four correlations
        # and differ only in what the antennas' receptors carry.
        codes, products, perm = pol_axis(4, -1, -1)
        @test codes == [-1, -2, -3, -4]
        @test products == ["PP", "QQ", "PQ", "QP"]
        @test products[perm] == ["PP", "PQ", "QP", "QQ"]
        @test pol_axis(4, -5, -1)[2:3] == (products, perm)

        # A subset keeps the canonical relative order.
        @test pol_axis(2, -1, -1) == ([-1, -2], ["PP", "QQ"], [1, 2])
        @test pol_axis(1, -1, -1) == ([-1], ["PP"], [1])

        # Stokes parameters are combinations, not correlations.
        @test_throws "not an AIPS correlation code" pol_axis(4, 1, 1)
        @test_throws "not an AIPS correlation code" pol_axis(4, -7, -1)
        @test_throws "repeats a correlation product" pol_axis(2, -1, 0)
    end

    @testset "polarization is per baseline" begin
        products = ["PP", "QQ", "PQ", "QP"]
        circular, linear = ("R", "L"), ("X", "Y")

        # Within one basis the familiar labels come back.
        @test pol_labels(products, circular, circular) == ["RR", "LL", "RL", "LR"]
        @test pol_labels(products, linear, linear) == ["XX", "YY", "XY", "YX"]
        # A mixed baseline resolves to products no array-wide axis label can
        # name, which is the reason the axis stays generic.
        @test pol_labels(products, circular, linear) == ["RX", "LY", "RY", "LX"]
        # The first letter is the first antenna's receptor, so swapping the
        # antennas transposes the cross hands.
        @test pol_labels(products, linear, circular) == ["XR", "YL", "XL", "YR"]
    end

    @testset "baseline packing" begin
        @test baseline_antennas(258) == (1, 2)
        @test baseline_antennas(515) == (2, 3)
        @test baseline_antennas(256 * 10 + 10) == (10, 10)
        # 0 is not an antenna number, so neither half may be zero.
        @test_throws "does not decode to two antennas" baseline_antennas(256)
        @test_throws "does not decode to two antennas" baseline_antennas(5)
    end

    @testset "channel frequencies" begin
        @test channel_frequencies(1.0e9, 0.0, 1.0e6, 4) ==
            [1.0e9, 1.001e9, 1.002e9, 1.003e9]
        # A lower-sideband band has a negative channel width, so its channels
        # run down from the band edge.
        @test channel_frequencies(1.0e9, 1.0e8, -1.0e6, 4) ==
            [1.1e9, 1.099e9, 1.098e9, 1.097e9]
    end

    @testset "scan segmentation" begin
        # Rows one second apart, with a source change and two longer gaps, so
        # that each of the two conditions breaks a scan on its own.
        jd = 2.4585e6 .+ [0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 100.0, 101.0] ./ 86400
        source_id = Int32[1, 1, 2, 2, 2, 2, 2, 2]

        # Source changes after row 2; gaps of 7 s and 89 s exceed the
        # threshold, so both also break a scan.
        @test scan_ranges(source_id, jd; gap = 5.0) == [1:2, 3:4, 5:6, 7:8]
        # A gap longer than every interval leaves only the source change.
        @test scan_ranges(source_id, jd; gap = 1000.0) == [1:2, 3:8]
        @test scan_ranges(Int32[], Float64[]) == UnitRange{Int}[]
        @test_throws DimensionMismatch scan_ranges(source_id[1:3], jd)
    end
end

# The validation file is 24 GB and is not fetched by CI; when it is present,
# check the reader against quantities the file states in more than one place.
const QBAND_PATH = joinpath(
    @__DIR__, "..", "testdata", "BT164", "BT164A",
    "VLBA_BT164A_bt164aQband_BIN0_SRC0_0_260108T230453.idifits",
)

if isfile(QBAND_PATH)
    @testset "FITS-IDI reading on the VLBA validation file" begin
        hdus = fits(QBAND_PATH)
        ag = IDI.require_table(hdus, "ARRAY_GEOMETRY")
        uv = IDI.require_table(hdus, "UV_DATA")

        # Every string comes out of the file padded; nothing may survive it.
        @test IDI.card_string(hdus[1], "TELESCOP") == "VLBA"
        @test IDI.card_string(uv, "EQUINOX") == "J2000"
        @test IDI.card_string(uv, "WEIGHTYP") == "CORRELAT"
        @test IDI._trim.(IDI.column(ag.data, :ANNAME)) ==
            ["BR", "FD", "HN", "KP", "LA", "MK", "NL", "OV", "PT", "SC"]

        no_stkd = Int(IDI.card(ag, "NO_STKD", 0))
        no_band = Int(IDI.card(ag, "NO_BAND", 0))
        no_chan = Int(IDI.card(ag, "NO_CHAN", 0))
        ref_freq = Float64(IDI.card(ag, "REF_FREQ", 0.0))
        @test (no_stkd, no_band, no_chan) == (4, 4, 256)

        _, products, perm = pol_axis(
            no_stkd, Int(IDI.card(ag, "STK_1", -1)),
            round(Int, Float64(IDI.card(uv, "CDELT2", -1.0))),
        )
        @test products == ["PP", "QQ", "PQ", "QP"]
        @test products[perm] == ["PP", "PQ", "QP", "QQ"]
        # The VLBA is homogeneously circular, so every baseline resolves the
        # same way and the familiar labels come back.
        an = IDI.require_table(hdus, "ANTENNA").data
        bases = collect(
            zip(
                IDI._trim.(IDI.column(an, :POLTYA)),
                IDI._trim.(IDI.column(an, :POLTYB)),
            )
        )
        @test all(==(("R", "L")), bases)
        @test pol_labels(products, bases[1], bases[2]) == ["RR", "LL", "RL", "LR"]

        # The four bands tile 512 MHz above REF_FREQ in 128 MHz steps.
        fq = IDI.require_table(hdus, "FREQUENCY").data
        bandfreq = IDI.rowvalues(IDI.column(fq, :BANDFREQ), 1)
        chwidth = IDI.rowvalues(IDI.column(fq, :CH_WIDTH), 1)
        @test bandfreq == [0.0, 1.28e8, 2.56e8, 3.84e8]
        @test channel_frequencies(ref_freq, bandfreq[1], chwidth[1], no_chan)[1] ==
            ref_freq
        @test channel_frequencies(ref_freq, bandfreq[4], chwidth[4], no_chan)[end] ≈
            ref_freq + 3.84e8 + 255 * 5.0e5

        src = IDI.require_table(hdus, "SOURCE").data
        @test IDI._trim.(IDI.column(src, :SOURCE)) == ["3C273", "M87", "3C279"]
        # DiFX writes the catalogue position into the apparent-place columns
        # rather than precessing it.
        @test IDI.column(src, :RAAPP) == IDI.column(src, :RAEPO)
        @test IDI.column(src, :DECAPP) == IDI.column(src, :DECEPO)

        rows = read_uv_rows(uv)
        @test length(rows.jd) == 752_160
        # SORT = 'T*' asserts time order; the conversion relies on it.
        @test IDI.card_string(uv, "SORT") == "T*"
        jd = rows.jd .+ rows.day
        @test issorted(jd)
        @test sort(unique(rows.source_id)) == Int32[1, 2, 3]
        @test unique(rows.integration_time) == [0.5]

        # 10 antennas including autocorrelations.
        pairs = baseline_antennas.(rows.baseline)
        @test length(unique(pairs)) == 55
        @test all(p -> 1 <= p[1] <= 10 && 1 <= p[2] <= 10, pairs)

        # The primary header records one DIFXJOB per scan, which is an
        # independent count of what the gap heuristic has to recover.
        njobs = count(c -> startswith(string(c), "HISTORY DIFXJOB"), hdus[1].cards)
        @test length(scan_ranges(rows.source_id, jd)) == njobs
    end
end
