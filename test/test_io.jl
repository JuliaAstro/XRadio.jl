@testset "synthetic data conforms" begin
    ps = synthetic_set(2)
    @test isvalid(ps)
    @test isempty(check(ps))
    @test length(ps) == 2
    @test collect(keys(ps)) == [:test_00, :test_01]
end

@testset "round trip through zarr v$fmt" for fmt in (2, 3)
    ps = synthetic_set(2)
    path = joinpath(mktempdir(), "test_v$fmt.ps.zarr")
    @test write(path, ps; zarr_format = fmt) == path
    @test isdir(path)

    back = read(ProcessingSet, path)

    @test length(back) == length(ps)
    @test collect(keys(back)) == collect(keys(ps))

    for name in keys(ps)
        a, b = ps[name], back[name]
        @testset "$name" begin
            @test sort(collect(keys(a))) == sort(collect(keys(b)))
            @test sort(collect(keys(branches(a)))) == sort(collect(keys(branches(b))))

            # The science arrays survive exactly, dimensions and all.
            for layer in (:visibility, :flag, :weight, :uvw)
                @test parent(a[layer]) == parent(b[layer])
                @test map(msv4name, dims(a[layer])) == map(msv4name, dims(b[layer]))
            end

            # Coordinates keep both their values and their attributes: the
            # integration time and channel width live there, and dropping them
            # would lose information without any array changing.
            @test frequencies(a) == frequencies(b)
            @test polarizations(a) == polarizations(b)
            @test times(a) == times(b)
            @test spectralwindow(a) == spectralwindow(b)
            fa = metadata(lookup(a[:visibility], Frequency))
            fb = metadata(lookup(b[:visibility], Frequency))
            @test XRadio.value(fa[:channel_width]) == XRadio.value(fb[:channel_width])
            ta = metadata(lookup(a[:visibility], Ti))
            tb = metadata(lookup(b[:visibility], Ti))
            @test XRadio.value(ta[:integration_time]) == XRadio.value(tb[:integration_time])

            @test antennas(a) == antennas(b)
            @test parent(a.antenna[:antenna_position]) == parent(b.antenna[:antenna_position])
            @test metadata(a)[:observation_info] == metadata(b)[:observation_info]
            @test datagroup(a) == datagroup(b)
        end
    end
end

@testset "consolidated metadata v$fmt" for fmt in (2, 3)
    ps = synthetic_set(2)
    path = joinpath(mktempdir(), "consolidated_v$fmt.ps.zarr")
    write(path, ps; zarr_format = fmt)

    # Zarr v2 keeps consolidated metadata in a `.zmetadata` file of its own;
    # v3 keeps it in a `consolidated_metadata` field of the root `zarr.json`.
    if fmt == 2
        @test isfile(joinpath(path, ".zmetadata"))
    else
        @test occursin("consolidated_metadata", read(joinpath(path, "zarr.json"), String))
    end

    # Only v2's is read back: `Zarr.jl` returns an empty group from a v3
    # consolidated read, so a v3 store is opened node by node.
    store = XRadio._zarr_group(path).storage
    @test (store isa Zarr.ConsolidatedStore) == (fmt == 2)

    # Reading through it gives the same processing set as reading the nodes.
    back = read(ProcessingSet, path)
    plain = read(ProcessingSet, Zarr.zopen(path))
    @test collect(keys(back)) == collect(keys(ps)) == collect(keys(plain))
    for name in keys(ps)
        @test metadata(back[name]) == metadata(plain[name])
        @test parent(back[name][:visibility]) == parent(ps[name][:visibility])
        @test antennas(back[name]) == antennas(ps[name])
    end
end

@testset "a store without consolidated metadata still opens" begin
    # Stores written by other tools have none, so its absence is not an error.
    ps = synthetic_set(1)
    path = joinpath(mktempdir(), "unconsolidated.ps.zarr")
    write(path, ps)
    rm(joinpath(path, ".zmetadata"))

    @test !(XRadio._zarr_group(path).storage isa Zarr.ConsolidatedStore)
    back = read(ProcessingSet, path)
    @test collect(keys(back)) == collect(keys(ps))
    @test parent(back[1][:visibility]) == parent(ps[1][:visibility])
end

@testset "a written store conforms" begin
    # Zarr v2 is the format this package can write fully to the standard; see
    # the `zarr_format` keyword of `write`.
    ps = synthetic_set(1)
    path = joinpath(mktempdir(), "conforming.ps.zarr")
    write(path, ps; zarr_format = 2)
    @test isempty(check(read(ProcessingSet, path)))
end

@testset "auxiliary coordinates are marked as coordinates" begin
    # `field_name` and `scan_name` are coordinates that are not dimensions.
    # xarray only knows that from the `coordinates` attribute on the data
    # variables, and without it reads them back as data variables — which this
    # package's own reader would not notice, because it labels from the
    # dimension names instead. So the attribute is checked directly.
    ps = synthetic_set(1)
    path = joinpath(mktempdir(), "coords.ps.zarr")
    write(path, ps)

    group = Zarr.zopen(path)["test_00"]
    listed = split(group["VISIBILITY"].attrs["coordinates"])
    for name in ("field_name", "scan_name", "baseline_antenna1_name", "baseline_antenna2_name")
        @test name in listed
    end
    # A coordinate does not list itself.
    @test !haskey(group["scan_name"].attrs, "coordinates")

    # Every non-index coordinate is named, including one whose dimensions no
    # data variable spans.
    antenna = Zarr.zopen(path)["test_00"]["antenna_xds"]
    @test "polarization_type" in split(antenna["ANTENNA_POSITION"].attrs["coordinates"])
end

@testset "lazy and eager reads agree" begin
    ps = synthetic_set(1)
    path = joinpath(mktempdir(), "lazy.ps.zarr")
    write(path, ps)

    lazy = open(ProcessingSet, path)
    eager = read(ProcessingSet, path)

    # A lazy layer is still backed by the store, but indexing it gives the
    # same numbers as the array that was read up front.
    @test parent(lazy[1][:visibility]) isa Zarr.ZArray
    @test parent(eager[1][:visibility]) isa Array
    @test lazy[1][:visibility][:, :, :, :] == parent(eager[1][:visibility])
end

@testset "element types" begin
    # The standard is itself flexible about the science arrays: any of
    # `float16`, `float32` and `float64` may carry weights and baseline
    # coordinates, and either complex width may carry visibilities.
    @test isvalid(synthetic_set(1; T = Float32, C = ComplexF64))
    @test isvalid(synthetic_set(1; T = Float16, C = ComplexF32))

    # It is not flexible about the metadata quantities, which are `float64`.
    # In memory there is no reason to refuse a narrower one, so the two checks
    # part company here: a strict check reports it, a structural one does not.
    ps = synthetic_set(1)
    ant = ps[1].antenna
    ant[:antenna_position] = map(Float32, ant[:antenna_position])

    @test !isvalid(ps)
    @test isvalid(ps; strict = false)
    issues = check(ps)
    @test length(issues) == 1
    @test occursin("ANTENNA_POSITION", only(issues).path)
    @test occursin("element type is Float32", only(issues).message)
end

@testset "selection by label" begin
    ps = synthetic_set(1; ntime = 4, nfreq = 5)
    ms = ps[1]
    vis = ms[:visibility]

    # Whichever end of the tuple a dimension sits at, selecting by label finds
    # the same data — which is what makes the reversed axis order harmless.
    @test size(vis[Polarization(At("XX"))]) == (5, 3, 4)
    @test size(vis[Frequency(1:2)]) == (2, 2, 3, 4)
    @test vis[Polarization(At("YY")), Frequency(At(2.3e11))] ==
        vis[Polarization(2), Frequency(1)]
end

@testset "processing set is a collection" begin
    ps = synthetic_set(3)
    @test length(ps) == 3
    @test !isempty(ps)
    @test eltype(collect(ps)) == MeasurementSet
    @test ps[1] === ps[:test_00] === ps.test_00
    @test haskey(ps, :test_01)
    @test length(filter(ms -> Intent("OBSERVE_TARGET#ON_SOURCE") in scan_intents(ms), ps)) == 3
    @test length(filter(ms -> isempty(scan_intents(ms)), ps)) == 0
    @test occursin("ProcessingSet", sprint(show, MIME("text/plain"), ps))
end

@testset "opening something that is not a processing set" begin
    dir = mktempdir()
    path = joinpath(dir, "plain.zarr")
    Zarr.zgroup(Zarr.DirectoryStore(path), "", 2)
    @test_throws "not an MSv4 processing set" open(ProcessingSet, path)
end

@testset "writing will not clobber" begin
    ps = synthetic_set(1)
    path = joinpath(mktempdir(), "once.ps.zarr")
    write(path, ps)
    @test_throws "pass `overwrite = true`" write(path, ps)
    @test write(path, ps; overwrite = true) == path
end

@testset "scan intents" begin
    @testset "the ASDM enumerations" begin
        # Both are closed, and both define TEST and UNSPECIFIED, which is why
        # they are separate modules.
        @test length(instances(ScanIntents.T)) == 30
        @test length(instances(SubscanIntents.T)) == 11
        @test ScanIntents.TEST isa ScanIntents.T
        @test SubscanIntents.TEST isa SubscanIntents.T
    end

    @testset "parsing and round-trip" begin
        i = Intent("OBSERVE_TARGET#ON_SOURCE")
        @test i == Intent(ScanIntents.OBSERVE_TARGET, SubscanIntents.ON_SOURCE)
        @test string(i) == "OBSERVE_TARGET#ON_SOURCE"

        # An omitted subintent is not `UNSPECIFIED`: the first says the store
        # names no subintent, the second says it names that one. Writing back
        # what was read must not turn one into the other.
        bare = Intent("UNSPECIFIED")
        @test bare.subintent === nothing
        @test string(bare) == "UNSPECIFIED"
        @test bare != Intent(ScanIntents.UNSPECIFIED, SubscanIntents.UNSPECIFIED)
        @test string(Intent(ScanIntents.UNSPECIFIED, SubscanIntents.UNSPECIFIED)) ==
            "UNSPECIFIED#UNSPECIFIED"

        @test all(s -> string(Intent(s)) == s,
            ("OBSERVE_TARGET#ON_SOURCE", "CALIBRATE_BANDPASS#UNSPECIFIED",
                "CALIBRATE_DIFFGAIN", "TEST#TEST"))
    end

    @testset "a value outside the enumerations throws" begin
        # This exact misspelling appears in the reference implementation's own
        # schema docstring, where it is silently just another string.
        @test_throws ArgumentError Intent("OBSERVE_TARGE#UNSPECIFIED")
        @test_throws "`OBSERVE_TARGE` is not an ASDM scan intent" Intent(
            "OBSERVE_TARGE#UNSPECIFIED"
        )
        @test_throws "the 30 it defines are" Intent("OBSERVE_TARGE#UNSPECIFIED")

        @test_throws "`ONSOURCE` is not an ASDM subscan intent" Intent(
            "OBSERVE_TARGET#ONSOURCE"
        )
        @test_throws "the 11 it defines are" Intent("OBSERVE_TARGET#ONSOURCE")

        @test_throws "is not a scan intent" Intent("A#B#C")
    end
end

@testset "accessors" begin
    ms = synthetic_set(1)[1]
    @test antennas(ms) == ["AA", "BB", "CC"]
    @test scans(ms) == ["1"]
    @test fields(ms) == ["field_0"]
    @test scan_intents(ms) == [Intent(ScanIntents.OBSERVE_TARGET, SubscanIntents.ON_SOURCE)]
    @test spectral_window_intents(ms) == ["CONTINUUM"]
    @test polarizations(ms) == ["XX", "YY"]
    @test spectralwindow(ms) == "spw_0"
    @test length(frequencies(ms)) == 5

    # The data group resolves which array is meant, rather than the name being
    # assumed.
    @test correlated(ms) === ms[:visibility]
    @test flags(ms) === ms[:flag]
    @test weights(ms) === ms[:weight]
    @test uvw(ms) === ms[:uvw]
    @test_throws KeyError datagroup(ms, :nonexistent)
end

@testset "baseline selection" begin
    ms = synthetic_set(1)[1]

    # A selector holding a vector compares by identity, so assertions about a
    # multi-baseline selection are made on the ids it names.
    selected(sel) = DimensionalData.val(DimensionalData.val(sel))

    @test baselines(ms) == [("AA", "BB"), ("AA", "CC"), ("BB", "CC")]

    @testset "the three key shapes agree" begin
        @test baseline(ms, ("AA", "CC")) == baseline(ms, "AA-CC") == baseline(ms, 2)
        # MSv4 records each pair once in whichever order the correlator gave,
        # so matching cannot depend on the caller guessing that order.
        @test baseline(ms, ("CC", "AA")) == baseline(ms, ("AA", "CC"))
        @test baseline(ms, "CC-AA") == baseline(ms, ("AA", "CC"))
    end

    @testset "the selector slices" begin
        vis = correlated(ms)
        @test size(vis[baseline(ms, ("AA", "BB"))]) == (2, 5, 4)
        @test vis[baseline(ms, ("AA", "BB"))] == vis[BaselineID(At(1))]
        @test size(uvw(ms)[baselines(ms, "AA")]) == (3, 2, 4)
        @test size(uvw(ms)[baselines(ms, "CC")]) == (3, 2, 4)
    end

    @testset "an absent baseline names what is present" begin
        @test_throws ArgumentError baseline(ms, ("AA", "ZZ"))
        @test_throws "no baseline AA-ZZ" baseline(ms, ("AA", "ZZ"))
        @test_throws "its baselines are over AA, BB, CC" baseline(ms, ("AA", "ZZ"))

        @test_throws "no baseline_id 9" baseline(ms, 9)
        @test_throws "3 baselines, with ids from 1 to 3" baseline(ms, 9)

        @test_throws "names antenna ZZ" baselines(ms, "ZZ")
        @test_throws "its baselines are over AA, BB, CC" baselines(ms, "ZZ")
    end

    @testset "a dash label spells exactly two names" begin
        @test_throws "is not a baseline label" baseline(ms, "AA")
        @test_throws "is not a baseline label" baseline(ms, "AA-BB-CC")
        @test_throws "Pass the pair as a tuple" baseline(ms, "AA-BB-CC")
    end

    @testset "autocorrelations" begin
        # A converted store holds them, so `auto = false` is not a rare path.
        selfed = synthetic_set(1)[1]
        base = dims(selfed, BaselineID)
        selfed[:baseline_antenna1_name] = DimArray(["AA", "AA", "BB"], (base,))
        selfed[:baseline_antenna2_name] = DimArray(["AA", "CC", "BB"], (base,))

        @test baselines(selfed) == [("AA", "AA"), ("AA", "CC"), ("BB", "BB")]
        @test baseline(selfed, ("AA", "AA")) == BaselineID(At(1))
        @test selected(baselines(selfed, "AA")) == [1, 2]
        @test selected(baselines(selfed, "AA"; auto = false)) == [2]
        # An antenna present only in its own autocorrelation selects nothing
        # once autocorrelations are dropped, which is a true answer rather than
        # an error.
        @test selected(baselines(selfed, "BB")) == [3]
        @test isempty(selected(baselines(selfed, "BB"; auto = false)))
    end

    @testset "a Measurement Set without baselines" begin
        bare = MeasurementSet(; metadata = Dict{Symbol, Any}(:type => "spectrum"))
        @test_throws "no `baseline_id` dimension" baseline(bare, 1)
        @test_throws "no `baseline_id` dimension" baselines(bare, "AA")

        noname = synthetic_set(1)[1]
        delete!(noname, :baseline_antenna1_name)
        @test_throws "names no baseline antennas" baselines(noname)
    end
end

# A store read from disk holds its antenna-name coordinates lazily, and states
# `baseline_id` from zero, so neither materialization nor 1-based ids can be
# assumed.
let path = joinpath(@__DIR__, "..", "testdata", "Antennae_North.cal.lsrk.split.v8.ps.zarr")
    if isdir(path)
        @testset "baseline selection on a store" begin
            ms = open(ProcessingSet, path)[1]
            bl = baselines(ms)
            @test length(bl) == 45
            @test allunique(bl)

            ids = collect(DimensionalData.val(dims(ms, BaselineID)))
            @test first(ids) == 0
            @test baseline(ms, first(ids)) == BaselineID(At(first(ids)))

            a, b = first(bl)
            @test baseline(ms, (b, a)) == baseline(ms, (a, b)) == BaselineID(At(first(ids)))

            @test size(uvw(ms)[baselines(ms, a)]) == (3, 9, size(uvw(ms), Ti))
            @test_throws "no baseline_id $(length(ids))" baseline(ms, length(ids))
        end
    end
end

@testset "selecting Measurement Sets by attribute" begin
    function varied(; field, spw, scan, ants)
        ms = synthetic_visibility()
        t = dims(ms, Ti)
        ms[:field_name] = DimArray(fill(field, length(t)), (t,))
        ms[:scan_name] = DimArray(fill(scan, length(t)), (t,))
        metadata(lookup(dims(ms, Frequency)))[:spectral_window_name] = spw
        branches(ms)[:antenna] = synthetic_antenna(ants)
        return ms
    end

    ps = ProcessingSet(
        OrderedCollections.OrderedDict{Symbol, MeasurementSet}(
            :a => varied(field = "M87", spw = "spw_0", scan = "1", ants = ["AA", "BB"]),
            :b => varied(field = "M87", spw = "spw_1", scan = "2", ants = ["AA", "CC"]),
            :c => varied(field = "3C273", spw = "spw_0", scan = "3", ants = ["BB", "CC"]),
        )
    )

    @testset "one attribute" begin
        @test collect(keys(select(ps; field = "M87"))) == [:a, :b]
        @test collect(keys(select(ps; spectralwindow = "spw_0"))) == [:a, :c]
        @test collect(keys(select(ps; scan = "3"))) == [:c]
        @test collect(keys(select(ps; antenna = "AA"))) == [:a, :b]
    end

    @testset "several attributes must all match" begin
        @test collect(keys(select(ps; field = "M87", spectralwindow = "spw_0"))) == [:a]
        @test isempty(select(ps; field = "3C273", scan = "1"))
    end

    @testset "a collection matches any of its values" begin
        @test collect(keys(select(ps; scan = ["1", "3"]))) == [:a, :c]
        @test collect(keys(select(ps; spectralwindow = ("spw_0", "spw_1")))) == [:a, :b, :c]
        @test collect(keys(select(ps; antenna = ["CC"]))) == [:b, :c]
    end

    @testset "the result is a processing set" begin
        got = select(ps; field = "M87")
        @test got isa ProcessingSet
        @test metadata(got)[:type] == "processing_set"
        # Selections chain, and chaining agrees with asking for both at once.
        @test keys(select(got; spectralwindow = "spw_1")) ==
            keys(select(ps; field = "M87", spectralwindow = "spw_1"))
        @test isempty(select(ps; field = "NOSUCH"))
        @test select(ps; field = "NOSUCH") isa ProcessingSet
        # No keyword asks for nothing, so nothing is excluded.
        @test keys(select(ps)) == keys(ps)
    end

    @testset "an unknown keyword is an error" begin
        # The whole risk of a keyword-driven selector is a typo silently
        # selecting everything.
        @test_throws ArgumentError select(ps; spectral_window = "spw_0")
        @test_throws "no keyword `spectral_window`" select(ps; spectral_window = "spw_0")
        @test_throws "accepts antenna, field, scan, spectralwindow" select(
            ps; spectral_window = "spw_0"
        )
        @test_throws "no keyword `source`, `freq`" select(ps; source = "M87", freq = 1)
        # `source` is not a synonym for `field`: MSv4 distinguishes them and
        # this package has no `sources` accessor to honor the distinction.
        @test_throws ArgumentError select(ps; source = "M87")
    end

    @testset "a keyword takes names" begin
        @test_throws "takes a name or a collection of names" select(ps; field = 3)
        @test_throws "got 3" select(ps; field = 3)
    end
end

@testset "transforming every Measurement Set" begin
    ps = synthetic_set(3)
    got = transform(ms -> delete!(deepcopy(ms), :weight), ps)

    @test got isa ProcessingSet
    @test keys(got) == keys(ps)
    @test metadata(got)[:type] == "processing_set"
    @test !haskey(got[1], :weight)
    # The originals are untouched; `transform` builds a new set from what `f`
    # returns rather than modifying `ps`.
    @test haskey(ps[1], :weight)

    # The result is a processing set, so it composes with the rest.
    @test length(filter(m -> true, got)) == 3
    @test keys(select(got; field = "field_0")) == keys(got)

    @testset "`f` must return a Measurement Set" begin
        @test_throws ArgumentError transform(gaincurves, ps)
        @test_throws "must return a MeasurementSet" transform(gaincurves, ps)
        @test_throws "it returned a DimStack" transform(gaincurves, ps)
        # The message points at the verb that does read values out.
        @test_throws "map(f, values(ps))" transform(gaincurves, ps)
    end

    # `map` is untouched: reading a value out of each Measurement Set still
    # gives them in order, as it does for any other iterable.
    @test map(spectralwindow, ps) == fill("spw_0", 3)
    @test map(gaincurves, values(ps)) isa Vector
end
