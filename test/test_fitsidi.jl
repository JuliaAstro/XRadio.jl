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

@testset "FITS-IDI fixture" begin
    dir = mktempdir()

    # Locate a UV_DATA row in the dense grid the fixture reports, going from
    # the file back to the grid rather than replaying the writer's row order.
    function cell(gt, rows, r)
        t = rows.day[r] * 86400
        k = findfirst(s -> isapprox(s, t; atol = 1.0e-6), gt.seconds)
        b = findfirst(==(rows.baseline[r]), gt.baseline_codes)
        return b, k
    end

    @testset "tables and cards" begin
        gt = synthetic_fitsidi(joinpath(dir, "full.idifits"); observer = "Someone",
            uv_observer = "Somebody Else")
        hdus = fits(gt.path)

        @test [IDI.card_string(h, "EXTNAME") for h in hdus[2:end]] ==
            ["ARRAY_GEOMETRY", "ANTENNA", "FREQUENCY", "SOURCE", "UV_DATA"]
        @test IDI.table(hdus, "FLAG") === nothing
        @test_throws "has no WEATHER table" IDI.require_table(hdus, "WEATHER")

        # The two headers carry different observers, which is the disagreement
        # the conversion has to resolve in favour of the primary.
        @test IDI.card_string(hdus[1], "OBSERVER") == "Someone"
        @test IDI.card_string(IDI.require_table(hdus, "UV_DATA"), "OBSERVER") ==
            "Somebody Else"

        uv = IDI.require_table(hdus, "UV_DATA")
        @test IDI.card_string(hdus[1], "TELESCOP") == gt.telescope
        @test IDI.card_string(uv, "EQUINOX") == "J2000"
        @test IDI.card_string(uv, "WEIGHTYP") == gt.weightyp
        @test IDI.card_string(uv, "SORT") == "T*"
        @test Int(IDI.card(uv, "NO_STKD", 0)) == gt.no_stkd
        @test Int(IDI.card(uv, "NO_BAND", 0)) == gt.nband
        @test Int(IDI.card(uv, "NO_CHAN", 0)) == gt.nchan
        @test Float64(IDI.card(uv, "REF_FREQ", 0.0)) == gt.ref_freq
    end

    @testset "metadata tables" begin
        gt = synthetic_fitsidi(
            joinpath(dir, "meta.idifits");
            antenna_names = ["AA", "BB", "CC"], mntsta = Int32[0, 1, 3],
            staxof = Float32[0 1.5 -2.5; 0 0 0.5; 0 -1 0], diameters = Float32[25, 30, 6],
            poltya = ["R", "R", "X"], poltyb = ["L", "L", "Y"],
        )
        hdus = fits(gt.path)

        ag = IDI.require_table(hdus, "ARRAY_GEOMETRY").data
        @test IDI._trim.(IDI.column(ag, :ANNAME)) == gt.antenna_names
        @test IDI.column(ag, :NOSTA) == gt.nosta
        @test IDI.column(ag, :MNTSTA) == gt.mntsta
        @test permutedims(IDI.column(ag, :STAXOF)) == gt.staxof
        @test IDI.column(ag, :DIAMETER) == gt.diameters
        # STABXYZ comes back as (antenna, axis), the transpose of the fixture's.
        @test permutedims(IDI.column(ag, :STABXYZ)) == gt.positions

        # A heterogeneous array: the third station carries linear feeds, so no
        # single array-wide polarization label describes every baseline.
        an = IDI.require_table(hdus, "ANTENNA").data
        @test IDI._trim.(IDI.column(an, :POLTYA)) == gt.poltya
        @test IDI._trim.(IDI.column(an, :POLTYB)) == gt.poltyb
        _, products, _ = pol_axis(gt.no_stkd, gt.stk_1, gt.cdelt2)
        @test pol_labels(products, ("R", "L"), ("X", "Y")) == ["RX", "LY", "RY", "LX"]

        fq = IDI.require_table(hdus, "FREQUENCY").data
        @test IDI.rowvalues(IDI.column(fq, :BANDFREQ), 1) == gt.band_freq
        @test IDI.rowvalues(IDI.column(fq, :CH_WIDTH), 1) == fill(gt.ch_width, gt.nband)
        for band in 1:gt.nband
            @test channel_frequencies(
                gt.ref_freq, gt.band_freq[band], gt.ch_width, gt.nchan
            ) == gt.frequencies[:, band]
        end

        src = IDI.require_table(hdus, "SOURCE").data
        @test IDI._trim.(IDI.column(src, :SOURCE)) == gt.source_names
        @test IDI.column(src, :RAEPO) == gt.ra
        @test IDI.column(src, :DECEPO) == gt.dec
        @test IDI.column(src, :RAAPP) == IDI.column(src, :RAEPO)
    end

    @testset "UV_DATA values" begin
        gt = synthetic_fitsidi(joinpath(dir, "uv.idifits"))
        hdus = fits(gt.path)
        uv = IDI.require_table(hdus, "UV_DATA")
        rows = read_uv_rows(uv)

        nrow = count(gt.present)
        @test length(rows.jd) == nrow
        @test all(==(gt.jd0), rows.jd)
        @test issorted(rows.jd .+ rows.day)
        @test unique(rows.integration_time) == [gt.inttim]
        @test unique(rows.freqid) == Int32[1]

        # Scans are recovered from source changes and time gaps alone.
        @test length(scan_ranges(rows.source_id, rows.jd .+ rows.day)) == gt.nscan

        flux = collect(IDI.column(uv.data, :FLUX))
        wts = collect(IDI.column(uv.data, :WEIGHT))
        for r in 1:nrow
            b, k = cell(gt, rows, r)
            @test rows.source_id[r] == gt.time_source[k]
            @test baseline_antennas(rows.baseline[r]) == gt.baselines[b]
            @test rows.uvw[r, :] == Float64.(gt.uvw[:, b, k])

            # FLUX runs COMPLEX fastest, then STOKES, then FREQ, then BAND.
            parts = reshape(flux[r, :], 2, gt.no_stkd, gt.nchan, gt.nband)
            @test complex.(parts[1, :, :, :], parts[2, :, :, :]) ==
                gt.visibility[:, :, :, b, k]
            # WEIGHT runs STOKES fastest, then BAND, one value per channel set.
            @test reshape(wts[r, :], gt.no_stkd, gt.nband) == gt.weight[:, :, b, k]
        end
    end

    @testset "absent baselines" begin
        # A baseline that drops out for one integration is simply not written,
        # which is the sparsity a dense MSv4 grid has to pad.
        gt = synthetic_fitsidi(
            joinpath(dir, "gaps.idifits"); present = (b, k) -> !(b == 1 && k == 2)
        )
        rows = read_uv_rows(IDI.require_table(fits(gt.path), "UV_DATA"))

        @test count(gt.present) == length(gt.baselines) * length(gt.seconds) - 1
        @test length(rows.jd) == count(gt.present)
        @test all(isnan, real(gt.visibility[:, :, :, 1, 2]))
        @test all(isnan, gt.weight[:, :, 1, 2])
        @test all(isnan, gt.uvw[:, 1, 2])
        # No row carries the missing (baseline, time) pair.
        missing_time = gt.seconds[2] / 86400
        @test !any(
            r -> rows.baseline[r] == gt.baseline_codes[1] &&
                isapprox(rows.day[r], missing_time; atol = 1.0e-9),
            eachindex(rows.baseline)
        )
    end

    @testset "a single band and a single source" begin
        # One band collapses the per-band columns to scalars, and one source
        # leaves the time gap as the only thing that can break a scan.
        gt = synthetic_fitsidi(
            joinpath(dir, "single.idifits");
            band_freq = [0.0], source_names = ["ONLY"], ra = [10.0], dec = [-5.0],
            scans = [(; source = 1, start = 0.0, ntime = 2),
                (; source = 1, start = 900.0, ntime = 2)],
        )
        hdus = fits(gt.path)

        fq = IDI.require_table(hdus, "FREQUENCY").data
        @test IDI.column(fq, :BANDFREQ) isa AbstractVector
        @test IDI.rowvalues(IDI.column(fq, :BANDFREQ), 1) == [0.0]

        rows = read_uv_rows(IDI.require_table(hdus, "UV_DATA"))
        @test unique(rows.source_id) == Int32[1]
        @test length(scan_ranges(rows.source_id, rows.jd .+ rows.day)) == gt.nscan
    end

    @testset "FLAG table" begin
        gt = synthetic_fitsidi(
            joinpath(dir, "flagged.idifits");
            flags = [
                (; ants = (2, 0), timerange = (0.0, 1.0), reason = "ANTENNA DOWN"),
                (; source_id = 2, chans = (1, 2), pflags = (1, 0, 0, 0)),
            ],
        )
        fl = IDI.require_table(fits(gt.path), "FLAG").data

        # 0 is a wildcard, so (2, 0) flags every baseline touching antenna 2.
        @test IDI.rowvalues(IDI.column(fl, :ANTS), 1) == Int32[2, 0]
        @test IDI.rowvalues(IDI.column(fl, :TIMERANG), 1) == [0.0, 1.0]
        @test IDI._trim(IDI.column(fl, :REASON)[1]) == "ANTENNA DOWN"
        @test IDI.column(fl, :SOURCE_ID) == Int32[0, 2]
        @test IDI.rowvalues(IDI.column(fl, :CHANS), 2) == Int32[1, 2]
        @test IDI.rowvalues(IDI.column(fl, :PFLAGS), 2) == Int32[1, 0, 0, 0]
        @test IDI.rowvalues(IDI.column(fl, :BANDS), 1) == ones(Int32, gt.nband)
    end
end

@testset "FITS-IDI conversion" begin
    dir = mktempdir()

    # A released date is invented when the file cannot supply one, and every
    # conversion here passes a real one so the warning does not drown the rest.
    to_msv4(gt, name; kw...) = fitsidi2msv4(
        gt.path, joinpath(dir, name);
        release_date = "2030-01-01T00:00:00.000", kw...
    )

    conversion_issues(ps) = sort!([i.message for i in check(ps)])

    @testset "the dense grid holds the values that went in" begin
        gt = synthetic_fitsidi(joinpath(dir, "grid.idifits"))
        ps = open(ProcessingSet, to_msv4(gt, "grid.ps.zarr"))

        # Spectral window and scan: two bands and two scans give four.
        @test collect(keys(ps)) == [:grid_00, :grid_01, :grid_02, :grid_03]
        @test isempty(conversion_issues(ps))

        # Band 1 of scan 1. The permutation putting the parallel hands first
        # runs the other way here: MSv4 polarization `p` came from the file's
        # `perm[p]`.
        _, _, perm = pol_axis(gt.no_stkd, gt.stk_1, gt.cdelt2)
        ms, band, times = ps[:grid_00], 1, 1:3
        vis, weight, uvw = ms[:visibility], ms[:weight], ms[:uvw]

        @test lookup(vis, Polarization) == ["RR", "RL", "LR", "LL"]
        @test lookup(vis, Frequency) == gt.frequencies[:, band]
        @test lookup(vis, Ti) == datetime2unix.(gt.datetimes[times])
        @test lookup(vis, BaselineId) == 1:length(gt.baselines)
        @test collect(ms[:baseline_antenna1_name]) ==
            [gt.antenna_names[a] for (a, _) in gt.baselines]
        @test collect(ms[:baseline_antenna2_name]) ==
            [gt.antenna_names[b] for (_, b) in gt.baselines]
        @test collect(ms[:field_name]) == fill("SRC0", 3)
        @test collect(ms[:scan_name]) == fill("1", 3)

        @test parent(vis) == gt.visibility[perm, :, band, :, times]
        @test !any(ms[:flag])

        # UU/VV/WW are seconds of light travel on disk and metres in MSv4.
        @test parent(uvw) ≈ gt.uvw[:, :, times] .* 299792458.0

        # `WEIGHTYP = "CORRELAT"` weights are validity fractions, so they are
        # thermal weights only after `2 Δν τ η²`, and they are per band on disk
        # but per channel in MSv4.
        scale = 2 * gt.ch_width * gt.inttim
        for c in 1:gt.nchan
            @test parent(weight)[:, c, :, :] ≈
                gt.weight[perm, band, :, times] .* Float32(scale)
        end

        # The frequency axis describes the one spectral window it came from.
        fmeta = metadata(lookup(vis, Frequency))
        @test fmeta[:spectral_window_name] == "band_0"
        @test fmeta[:spectral_window_intents] == ["UNSPECIFIED"]
        @test XRadio.value(fmeta[:channel_width]) == gt.ch_width
        @test XRadio.value(fmeta[:reference_frequency]) == gt.ref_freq

        tmeta = metadata(lookup(vis, Ti))
        @test (tmeta[:scale], tmeta[:format], tmeta[:units]) == ("utc", "unix", "s")
        @test XRadio.value(tmeta[:integration_time]) == gt.inttim

        # `EQUINOX = "J2000"` is `fk5`, and the UVW frame is the phase centre's.
        @test metadata(uvw)[:frame] == "fk5"

        obs = metadata(ms)[:observation_info]
        @test obs[:observer] == [gt.observer]
        @test obs[:project_UID] == gt.obscode
        @test obs[:release_date] == "2030-01-01T00:00:00.000"
        @test metadata(ms)[:processor_info][:type] == "CORRELATOR"
        @test occursin("grid.idifits", metadata(ms)[:data_groups][:base][:description])

        # The second band is the second spectral window of the same scan.
        @test lookup(ps[:grid_01][:visibility], Frequency) == gt.frequencies[:, 2]
        @test collect(ps[:grid_02][:scan_name]) == fill("2", 2)
    end

    @testset "the antennas come from ARRAY_GEOMETRY and ANTENNA" begin
        gt = synthetic_fitsidi(
            joinpath(dir, "ant.idifits");
            array_center = (10.0, -20.0, 30.0),
            mntsta = Int32[0, 4, 1], staxof = Float32[0 1.5 0; 0 0 0.5; 0 -1 0],
            diameters = Float32[25, 30, 6],
        )
        ant = branches(open(ProcessingSet, to_msv4(gt, "ant.ps.zarr"))[:ant_00])[:antenna]

        @test lookup(ant[:antenna_position], AntennaName) == gt.antenna_names
        @test collect(ant[:station_name]) == gt.antenna_names
        @test collect(ant[:telescope_name]) == fill(gt.array_name, 3)

        # `STABXYZ` is measured from the array center, which MSv4 is not.
        @test parent(ant[:antenna_position]) ≈ gt.positions .+ [10.0, -20.0, 30.0]
        @test metadata(ant[:antenna_position])[:coordinate_system] == "geocentric"

        @test collect(ant[:mount]) == ["ALT-AZ", "ALT-AZ+NASMYTH-R", "EQUATORIAL"]
        @test collect(ant[:antenna_dish_diameter]) == [25.0, 30.0, 6.0]
        @test lookup(ant[:polarization_type], ReceptorLabel) == ["pol_0", "pol_1"]
        @test parent(ant[:polarization_type]) == ["R" "R" "R"; "L" "L" "L"]
        @test all(iszero, ant[:antenna_receptor_angle])
    end

    @testset "a partition describes the antennas that observed it" begin
        # The third station drops out of the second scan entirely, so that
        # scan's partitions describe two antennas rather than three.
        gt = synthetic_fitsidi(
            joinpath(dir, "down.idifits");
            baselines = [(1, 2), (1, 3), (2, 3)],
            present = (b, k) -> !(k > 3 && b in (2, 3)),
        )
        ps = open(ProcessingSet, to_msv4(gt, "down.ps.zarr"))

        @test lookup(branches(ps[:down_00])[:antenna][:mount], AntennaName) ==
            ["AA", "BB", "CC"]
        @test lookup(branches(ps[:down_02])[:antenna][:mount], AntennaName) == ["AA", "BB"]
    end

    @testset "the fields come from SOURCE" begin
        gt = synthetic_fitsidi(joinpath(dir, "field.idifits"))
        ps = open(ProcessingSet, to_msv4(gt, "field.ps.zarr"))

        # A partition covers one scan, and a scan one source.
        fs = branches(ps[:field_00])[:field_and_source_base]
        @test lookup(fs[:source_name], FieldName) == ["SRC0"]
        @test collect(fs[:source_name]) == ["SRC0"]
        @test parent(fs[:field_phase_center_direction]) ≈
            deg2rad.([gt.ra[1]; gt.dec[1];;])
        @test parent(fs[:source_direction]) ≈ parent(fs[:field_phase_center_direction])
        @test metadata(fs[:field_phase_center_direction])[:frame] == "fk5"

        # `field_name` on the time axis names a field of this dataset.
        @test collect(ps[:field_02][:field_name]) == fill("SRC1", 2)
        @test lookup(branches(ps[:field_02])[:field_and_source_base][:source_name],
            FieldName) == ["SRC1"]

        # The data group names the sub-dataset the standard expects to find.
        @test metadata(ps[:field_00])[:data_groups][:base][:field_and_source] ==
            "field_and_source_base_xds"
    end

    @testset "a mount MSv4 cannot name is written as OTHER" begin
        # MNTSTA 2 (orbiting) and 6 (other) have no MSv4 name, and 3 (X-Y) does.
        gt = synthetic_fitsidi(joinpath(dir, "mount.idifits"); mntsta = Int32[3, 2, 6])
        ps = @test_logs (:warn, r"no `mount` name") match_mode = :any open(
            ProcessingSet, to_msv4(gt, "mount.ps.zarr")
        )
        @test collect(branches(ps[:mount_00])[:antenna][:mount]) ==
            ["X-Y", "OTHER", "OTHER"]

        gt = synthetic_fitsidi(joinpath(dir, "mnt7.idifits"); mntsta = Int32[0, 0, 7])
        @test_throws "MNTSTA 7 names no mount" to_msv4(gt, "mnt7.ps.zarr")
    end

    @testset "absent rows become padding" begin
        # One baseline drops out of one integration, which FITS-IDI records by
        # omitting the row and MSv4 by padding the cell.
        gt = synthetic_fitsidi(
            joinpath(dir, "sparse.idifits"); present = (b, k) -> !(b == 2 && k == 2)
        )
        ms = open(ProcessingSet, to_msv4(gt, "sparse.ps.zarr"))[:sparse_00]

        @test all(isnan, ms[:visibility][:, :, 2, 2])
        @test all(ms[:flag][:, :, 2, 2])
        @test all(iszero, ms[:weight][:, :, 2, 2])
        @test all(isnan, ms[:uvw][:, 2, 2])
        # Nothing else is padded.
        @test count(ms[:flag]) == size(ms[:flag], 1) * size(ms[:flag], 2)
    end

    @testset "a heterogeneous array keeps its per-baseline polarizations" begin
        gt = synthetic_fitsidi(
            joinpath(dir, "mixed.idifits");
            poltya = ["R", "R", "X"], poltyb = ["L", "L", "Y"],
        )
        ms = open(ProcessingSet, to_msv4(gt, "mixed.ps.zarr"))[:mixed_00]

        # Two of the three stations are circular, so that is what names the
        # array-wide axis; the baselines to the third measure something else,
        # which only `polarization_mixed` can say.
        @test lookup(ms[:visibility], Polarization) == ["RR", "RL", "LR", "LL"]
        mixed = ms[:polarization_mixed]
        @test lookup(mixed, BaselineId) == 1:length(gt.baselines)
        @test collect(mixed[:, 1]) == ["RR", "RL", "LR", "LL"]   # AA-AA
        @test collect(mixed[:, 3]) == ["RX", "RY", "LX", "LY"]   # AA-CC
        @test collect(mixed[:, 6]) == ["XX", "XY", "YX", "YY"]   # CC-CC
    end

    @testset "partitioning" begin
        gt = synthetic_fitsidi(
            joinpath(dir, "part.idifits");
            band_freq = [0.0],
            scans = [(; source = 1, start = 0.0, ntime = 2),
                (; source = 1, start = 5.0, ntime = 2)],
        )

        # The two scans share a source, so only the time gap separates them.
        ps = open(ProcessingSet, to_msv4(gt, "scan.ps.zarr"; gap = 2.0))
        @test collect(keys(ps)) == [:part_00, :part_01]
        @test collect(ps[:part_01][:scan_name]) == ["2", "2"]

        # Splitting by source does not see the gap, so both scans land together
        # and the `scan_name` coordinate is what distinguishes them.
        ps = open(ProcessingSet, to_msv4(gt, "source.ps.zarr"; partition = (:source,)))
        @test length(ps) == 1
        @test collect(ps[:part_00][:scan_name]) == ["1", "1", "1", "1"]

        # No partitioning beyond the mandatory spectral window.
        ps = open(ProcessingSet, to_msv4(gt, "whole.ps.zarr"; partition = ()))
        @test length(ps) == 1
        @test size(ps[:part_00][:visibility], Ti) == 4

        @test_throws "partition must be" to_msv4(gt, "bad.ps.zarr"; partition = (:band,))
    end

    @testset "weights follow WEIGHTYP" begin
        normal = synthetic_fitsidi(joinpath(dir, "normal.idifits"); weightyp = "NORMAL")
        ms = open(ProcessingSet, to_msv4(normal, "normal.ps.zarr"))[:normal_00]
        # Already `1/σ²`, so nothing is applied.
        @test ms[:weight][1, 1, 1, 1] == normal.weight[1, 1, 1, 1]

        # `:validity` writes the column through whatever the keyword says.
        correlat = synthetic_fitsidi(joinpath(dir, "corr.idifits"))
        ms = open(
            ProcessingSet, to_msv4(correlat, "validity.ps.zarr"; weight_mode = :validity)
        )[:corr_00]
        @test ms[:weight][1, 1, 1, 1] == correlat.weight[1, 1, 1, 1]

        # Efficiency scales the whole conversion, not the relative weighting.
        ms = open(
            ProcessingSet, to_msv4(correlat, "eta.ps.zarr"; weight_efficiency = 0.5)
        )[:corr_00]
        scale = 2 * correlat.ch_width * correlat.inttim * 0.25
        @test ms[:weight][1, 1, 1, 1] ≈ correlat.weight[1, 1, 1, 1] * Float32(scale)

        @test_throws "weight_mode must be" to_msv4(correlat, "wm.ps.zarr"; weight_mode = :other)
    end

    @testset "element type" begin
        gt = synthetic_fitsidi(joinpath(dir, "eltype.idifits"))
        ms = open(
            ProcessingSet, to_msv4(gt, "f64.ps.zarr"; element_type = Float64)
        )[:eltype_00]
        @test eltype(ms[:visibility]) == ComplexF64
        @test eltype(ms[:weight]) == Float64
        # UVW does not follow `element_type`: Float32 cannot resolve a VLBI
        # baseline to the fraction of a wavelength the phase needs.
        @test eltype(ms[:uvw]) == Float64

        @test_throws "element_type must be" to_msv4(gt, "int.ps.zarr"; element_type = Int)
    end

    @testset "what the file cannot say" begin
        gt = synthetic_fitsidi(joinpath(dir, "invent.idifits"))

        # A sentinel release date could be read as a real embargo policy, so it
        # is the one invented field that warns.
        dest = joinpath(dir, "invent.ps.zarr")
        @test_logs (:warn, r"carries no release date") match_mode = :any fitsidi2msv4(
            gt.path, dest
        )
        ms = open(ProcessingSet, dest)[:invent_00]
        @test metadata(ms)[:observation_info][:release_date] ==
            "1858-11-17T00:00:00.000"
        @test metadata(ms[:scan_name])[:scan_intents] == ["UNSPECIFIED"]

        ms = open(
            ProcessingSet,
            to_msv4(gt, "intents.ps.zarr"; scan_intents = ["OBSERVE_TARGET#ON_SOURCE"],
                spectral_window_intents = ["CONTINUUM"])
        )[:invent_00]
        @test metadata(ms[:scan_name])[:scan_intents] == ["OBSERVE_TARGET#ON_SOURCE"]
        @test metadata(lookup(ms[:visibility], Frequency))[:spectral_window_intents] ==
            ["CONTINUUM"]
    end

    @testset "an unconvertible file is refused" begin
        # MSv4's `UVW` frame is one of icrs, fk5 and APP, so a B1950 phase
        # center has no conforming home; calling it fk5 would displace every
        # baseline by the angle between the equinoxes.
        gt = synthetic_fitsidi(joinpath(dir, "b1950.idifits"); equinox = "B1950")
        @test_throws "Precess the observation to J2000" to_msv4(gt, "b1950.ps.zarr")

        gt = synthetic_fitsidi(joinpath(dir, "odd.idifits"); equinox = "J1975")
        @test_throws "names no MSv4 frame" to_msv4(gt, "odd.ps.zarr")

        gt = synthetic_fitsidi(joinpath(dir, "twice.idifits"))
        to_msv4(gt, "twice.ps.zarr")
        @test_throws "pass `overwrite = true`" to_msv4(gt, "twice.ps.zarr")
        @test to_msv4(gt, "twice.ps.zarr"; overwrite = true) ==
            joinpath(dir, "twice.ps.zarr")
    end

    @testset "what the standard permits that the VLBA files do not show" begin
        # AIPS Memo 114r makes SOURCE_ID conditional on there being more than one
        # source, FREQID on more than one setup, and INTTIM optional outright.
        gt = synthetic_fitsidi(
            joinpath(dir, "spare.idifits");
            band_freq = [0.0], source_names = ["ONLY"], ra = [10.0], dec = [20.0],
            scans = [(; source = 1, start = 0.0, ntime = 2)],
            source_column = false, freqid_column = false, inttim_column = false,
        )
        # MSv4's time axis must state an integration time, so a file omitting
        # INTTIM can only be converted with one supplied.
        @test_throws "pass `integration_time`" to_msv4(gt, "spare.ps.zarr")
        ms = open(
            ProcessingSet, to_msv4(gt, "spare.ps.zarr"; integration_time = 2.5)
        )[:spare_00]
        @test XRadio.value(metadata(lookup(ms[:visibility], Ti))[:integration_time]) == 2.5
        @test collect(ms[:field_name]) == ["ONLY", "ONLY"]

        # `SOURCE`, `SOURCE_ID` and `SOURCE ID` are stated to be synonyms.
        for key in (:SOURCE, :SOURCE_ID, Symbol("SOURCE ID"))
            gt = synthetic_fitsidi(joinpath(dir, "src$key.idifits"); source_key = key)
            ms = open(ProcessingSet, to_msv4(gt, "src$key.ps.zarr"))[Symbol("src$(key)_00")]
            @test collect(ms[:field_name]) == fill("SRC0", 3)
        end
    end

    @testset "the standard's spellings, not the VLBA's" begin
        # Memo 114r spells the equinox `'1950.0B'` and `'J2000'`. Reading them
        # is what is at stake here; MSv4 cannot store an fk4 observation, so a
        # B1950 file gets no further than its SOURCE table.
        gt = synthetic_fitsidi(joinpath(dir, "b50.idifits"); equinox = "1950.0B")
        @test IDI.source_frames(IDI.require_table(fits(gt.path), "SOURCE")) ==
            fill("fk4", length(gt.source_names))

        # The VLBA writes no EQUINOX column but a numeric EPOCH, which the memo
        # says must be read as an equinox rather than as an epoch.
        gt = synthetic_fitsidi(joinpath(dir, "epoch.idifits"); equinox_column = :EPOCH)
        ps = open(ProcessingSet, to_msv4(gt, "epoch.ps.zarr"))
        @test metadata(ps[:epoch_00][:uvw])[:frame] == "fk5"

        gt = synthetic_fitsidi(joinpath(dir, "b50e.idifits");
            equinox = "1950.0B", equinox_column = :EPOCH)
        @test IDI.source_frames(IDI.require_table(fits(gt.path), "SOURCE")) ==
            fill("fk4", length(gt.source_names))
    end

    @testset "---NCP is a different projection, not a spelling" begin
        # Its w axis points at the north pole rather than at the source, so
        # reading it as ---SIN would be wrong without complaining.
        gt = synthetic_fitsidi(joinpath(dir, "ncp.idifits"); uvw_suffix = "---NCP")
        @test_throws "---NCP" to_msv4(gt, "ncp.ps.zarr")

        # The VLBA archive's `-L` suffix means ---SIN in seconds despite naming
        # wavelengths, which the memo tells readers to expect.
        gt = synthetic_fitsidi(joinpath(dir, "susp.idifits"); uvw_suffix = "-L")
        ms = open(ProcessingSet, to_msv4(gt, "susp.ps.zarr"))[:susp_00]
        @test parent(ms[:uvw]) ≈ gt.uvw[:, :, 1:3] .* 299792458.0
    end

    @testset "VIS_SCAL normalizes the visibilities" begin
        # A near-real-time correlator may write accumulated sums and leave the
        # division to the reader.
        gt = synthetic_fitsidi(joinpath(dir, "scal.idifits"); vis_scal = 4.0)
        ms = open(ProcessingSet, to_msv4(gt, "scal.ps.zarr"))[:scal_00]
        _, _, perm = pol_axis(gt.no_stkd, gt.stk_1, gt.cdelt2)
        @test parent(ms[:visibility]) ≈ gt.visibility[perm, :, 1, :, 1:3] ./ 4

        # Where the weights are true weights it squares into them, since it
        # divides into an amplitude. Validity fractions are untouched.
        normal = synthetic_fitsidi(
            joinpath(dir, "scaln.idifits"); vis_scal = 4.0, weightyp = "NORMAL"
        )
        ms = open(ProcessingSet, to_msv4(normal, "scaln.ps.zarr"))[:scaln_00]
        @test ms[:weight][1, 1, 1, 1] ≈ normal.weight[1, 1, 1, 1] * 16

        gt = synthetic_fitsidi(joinpath(dir, "scal0.idifits"); vis_scal = 0.0)
        @test_throws "must be positive" to_msv4(gt, "scal0.ps.zarr")
    end

    @testset "a correlator's own scan column is believed" begin
        # Two scans of one source with no time gap between them: inference
        # cannot see the boundary, and a stated scan number can.
        knobs = (; band_freq = [0.0], source_names = ["ONLY"], ra = [10.0], dec = [20.0],
            scans = [(; source = 1, start = 0.0, ntime = 2),
                (; source = 1, start = 2.0, ntime = 2)])
        gt = synthetic_fitsidi(joinpath(dir, "noscan.idifits"); knobs...)
        ps = open(ProcessingSet, to_msv4(gt, "noscan.ps.zarr"))
        @test length(ps) == 1

        gt = synthetic_fitsidi(joinpath(dir, "scan.idifits"); knobs..., scan_column = true)
        ps = open(ProcessingSet, to_msv4(gt, "scanned.ps.zarr"))
        @test length(ps) == 2
        @test collect(ps[:scan_00][:scan_name]) == ["1", "1"]
        @test collect(ps[:scan_01][:scan_name]) == ["2", "2"]
    end

    @testset "WEIGHT is tied to the COMPLEX axis" begin
        # With MAXIS1 = 2 the weights have nowhere else to live, so the column
        # is required rather than defaulted.
        gt = synthetic_fitsidi(joinpath(dir, "noweight.idifits"); weight_column = false)
        @test_throws "requires the column" to_msv4(gt, "noweight.ps.zarr")
    end

    @testset "one array at a time" begin
        # An antenna number is unique only within its array, and the
        # 256*a1 + a2 baseline packing carries no array number.
        gt = synthetic_fitsidi(joinpath(dir, "twoarr.idifits"); narrays = 2)
        @test_throws "defines 2 arrays" to_msv4(gt, "twoarr.ps.zarr")
    end

    @testset "the data matrix layout" begin
        gt = synthetic_fitsidi(joinpath(dir, "layout.idifits"))
        uv = IDI.require_table(fits(gt.path), "UV_DATA")

        layout = IDI.uv_matrix_layout(uv)
        @test layout == (; ncomplex = 2, no_stkd = gt.no_stkd, no_chan = gt.nchan,
            no_band = gt.nband)

        # One block covers every row, and the rows are the file's own order.
        block = IDI.read_uv_block(uv, 1:6)
        @test size(block.flux) == (2, gt.no_stkd, gt.nchan, gt.nband, 6)
        @test size(block.weight) == (gt.no_stkd, gt.nband, 6)
        @test complex.(block.flux[1, :, :, :, 1], block.flux[2, :, :, :, 1]) ==
            gt.visibility[:, :, :, 1, 1]
        @test block.weight[:, :, 1] == gt.weight[:, :, 1, 1]
        @test_throws "lie outside the table's" IDI.read_uv_block(uv, 1:100)
    end
end

# Which difference of antenna positions the `UVW` of each format is.
#
# MSv4 does not say. Its `UVW` definition calls the baseline both "from
# `baseline_antenna2_id` to `baseline_antenna1_id`" and "the difference
# `POSITION2 - POSITION1`", which are opposite claims, and the geometry it adds
# — `w` towards the source, `v` in the plane through source and pole, `u`
# towards increasing longitude — fixes the axes without fixing the sign. So the
# sign is measured from data instead: rotate a baseline vector built from the
# store's own antenna positions into that frame and see which way it points.
#
# Nothing here depends on MSv2. A sign error is silent and mirrors every image
# made from a converted store, so the two formats being asserted to agree is
# what licenses `fitsidi2msv4` to pass UU/VV/WW through unscaled.

# Greenwich mean sidereal time in degrees (IAU 1982). The candidates are 180
# degrees apart, so neither UT1-UTC nor nutation matters at this precision.
function gmst_degrees(jd)
    d = jd - 2451545.0
    T = d / 36525
    return mod(280.46061837 + 360.98564736629d + 0.000387933T^2 - T^3 / 38710000, 360)
end

# The (u, v, w) an earth-fixed baseline vector has for a phase centre, by the
# geometry both standards state.
function uvw_of(baseline, ra_degrees, dec_degrees, jd)
    H = deg2rad(gmst_degrees(jd) - ra_degrees)
    d = deg2rad(dec_degrees)
    x, y, z = baseline
    return [
        x * sin(H) + y * cos(H),
        -x * sin(d) * cos(H) + y * sin(d) * sin(H) + z * cos(d),
        x * cos(d) * cos(H) - y * cos(d) * sin(H) + z * sin(d),
    ]
end

# How far along the stored vector the predicted one reaches: +1 for the same
# direction, -1 for the opposite.
alignment(predicted, stored) = dot(predicted, stored) / dot(stored, stored)

const MSV4_PATH = joinpath(
    @__DIR__, "..", "testdata", "Antennae_North.cal.lsrk.split.v8.ps.zarr"
)

if isdir(MSV4_PATH)
    @testset "MSv4 `UVW` is the first antenna's position minus the second's" begin
        ms = open(ProcessingSet, MSV4_PATH)[1]
        antenna = branches(ms)[:antenna]
        position = antenna[:antenna_position]
        of = Dict(
            String(name) => Float64.(collect(position[:, At(name)]))
                for name in lookup(position, AntennaName)
        )
        centre = branches(ms)[:field_and_source_xds_base][:field_phase_center]
        ra, dec = rad2deg(centre[1]), rad2deg(centre[2])
        ant1 = collect(ms[:baseline_antenna1_name])
        ant2 = collect(ms[:baseline_antenna2_name])
        times = collect(lookup(ms[:uvw], Ti))

        ratios = [
            alignment(
                uvw_of(
                    of[ant1[b]] .- of[ant2[b]], ra, dec, times[k] / 86400 + 2440587.5
                ),
                Float64.(collect(ms[:uvw][:, b, k]))
            ) for b in eachindex(ant1), k in eachindex(times)
        ]
        # The standard's `POSITION2 - POSITION1` would put every one of these at
        # -1. Its "from antenna 2 to antenna 1" is the clause the data follows.
        @test all(r -> isapprox(r, 1; atol = 0.01), ratios)
        @test length(ratios) == length(ant1) * length(times)
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

        # The same measurement on the other side of the conversion: FITS-IDI's
        # baseline is the first antenna's position minus the second's, so it
        # agrees with MSv4 and `fitsidi2msv4` passes UU/VV/WW through unscaled.
        stab = IDI.column(ag.data, :STABXYZ)
        nosta = Int.(IDI.column(ag.data, :NOSTA))
        of = Dict(nosta[i] => Float64.(stab[i, :]) for i in eachindex(nosta))
        src = IDI.require_table(hdus, "SOURCE").data
        ra = Dict(
            Int(s) => Float64(r)
                for (s, r) in zip(IDI.column(src, :SOURCE_ID), IDI.column(src, :RAEPO))
        )
        dec = Dict(
            Int(s) => Float64(d)
                for (s, d) in zip(IDI.column(src, :SOURCE_ID), IDI.column(src, :DECEPO))
        )
        rows = read_uv_rows(uv)
        for r in (1, 5_000, 200_000, 500_000)
            a1, a2 = baseline_antennas(rows.baseline[r])
            a1 == a2 && continue
            predicted = uvw_of(
                of[a1] .- of[a2], ra[Int(rows.source_id[r])],
                dec[Int(rows.source_id[r])], rows.jd[r] + rows.day[r]
            )
            stored = collect(rows.uvw[r, :]) .* 299792458.0
            @test isapprox(alignment(predicted, stored), 1; atol = 0.05)
        end

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
