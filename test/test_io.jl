@testset "synthetic data conforms" begin
    ps = synthetic_set(2)
    @test isvalid(ps)
    @test isempty(check(ps))
    @test length(ps) == 2
    @test collect(keys(ps)) == [:test_00, :test_01]
end

@testset "a dimension outside the vocabulary round-trips" begin
    # Data MSv4 has no schema for brings its own dimensions, and both spellings
    # have to survive a write and a read: the one with a registered type, and
    # the one the generic `Dim` labels.
    ps = synthetic_set(1)
    gain = DimArray(
        randn(ComplexF64, 3, 2),
        (Dim{:solution_interval}([0.0, 1.0, 2.0]), Dim{:solution_basis}(["a", "b"])),
    )
    ps[1][:solution_gain] = gain

    path = joinpath(mktempdir(), "extension.ps.zarr")
    write(path, ps; validate = false)
    back = read(ProcessingSet, path)

    got = back[1][:solution_gain]
    @test parent(got) == parent(gain)
    @test map(msv4name, dims(got)) == (:solution_interval, :solution_basis)
    @test lookup(got, Dim{:solution_interval}) == [0.0, 1.0, 2.0]
    @test lookup(got, Dim{:solution_basis}) == ["a", "b"]

    # Writing back what was read reproduces the store, so a dimension's name
    # does not drift over a round trip.
    again = joinpath(mktempdir(), "extension_again.ps.zarr")
    write(again, back; validate = false)
    @test sort(collect(keys(read(ProcessingSet, again)[1]))) == sort(collect(keys(back[1])))

    # The extension is reported as a conformance issue rather than raised as an
    # error: whether the schema covers an array is a separate question from
    # whether the data can be described at all. What is reported is the
    # dimensions, not the array — an extra array is always allowed, but only
    # over dimensions the dataset already has.
    issues = check(back)
    @test any(i -> occursin("`solution_interval`", i.message), issues)
    @test any(i -> occursin("`solution_basis`", i.message), issues)
    @test !any(i -> occursin("solution_gain", i.path), issues)
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

# A copy of `ms` naming one scan throughout.
function scanned(ms::MeasurementSet, scan::AbstractString)
    out = deepcopy(ms)
    t = dims(out, Ti)
    out[:scan_name] = DimArray(
        fill(scan, length(t)), (t,); metadata = metadata(out[:scan_name])
    )
    return out
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

    @testset "an absent axis gives an empty vector of the right type" begin
        # The dimension supplies the element type when it is there; when it is
        # not, the type is what MSv4 declares the coordinate to be, so a caller
        # is not handed a `Vector{Any}` to widen everything it touches.
        bare = MeasurementSet(; metadata = Dict{Symbol, Any}(:type => "visibility"))
        @test frequencies(bare) == Float64[]
        @test times(bare) == Float64[]
        @test polarizations(bare) == String[]
        @test antennas(bare) == String[]
        @test frequencies(bare) isa Vector{Float64}
        @test times(bare) isa Vector{Float64}
        @test polarizations(bare) isa Vector{String}
        @test antennas(bare) isa Vector{String}

        # Selecting a single value drops the axis, which is the way a caller
        # meets this without building a bare Measurement Set.
        @test polarizations(ms[polarization = At("XX")]) isa Vector{String}

        # A present axis still answers with the store's own element type, which
        # a non-strict store may state more narrowly than the schema does.
        @test frequencies(ms) isa Vector{Float64}
        @test polarizations(ms) == ["XX", "YY"]
    end

    @testset "over a whole processing set" begin
        ps = synthetic_set(3)
        # Every member names the same scan and field, and each is reported once.
        @test scans(ps) == ["1"]
        @test fields(ps) == ["field_0"]

        varied = ProcessingSet(
            OrderedCollections.OrderedDict{Symbol, MeasurementSet}(
                :a => scanned(ms, "3"), :b => scanned(ms, "1"), :c => scanned(ms, "3"),
            )
        )
        # First occurrence, not sorted: the set's own order is the answer's.
        @test scans(varied) == ["3", "1"]
        @test scans(ProcessingSet(
            OrderedCollections.OrderedDict{Symbol, MeasurementSet}()
        )) == String[]
    end
end


@testset "the field-and-source dataset and the partition summary" begin
    ms = synthetic_set(1)[1]
    @test field_and_source(ms) === branches(ms)[:field_and_source_base]
    @test partition_info(ms) == (;
        spectral_window_name = "spw_0",
        spectral_window_intents = ["CONTINUUM"],
        field_name = ["field_0"],
        polarization_setup = ["XX", "YY"],
        scan_name = ["1"],
        source_name = ["source_0"],
        scan_intents = ["OBSERVE_TARGET#ON_SOURCE"],
        line_name = String[],
        data_group_name = :base,
    )

    @testset "name lists are sorted and hold each name once" begin
        out = deepcopy(ms)
        t = dims(out, Ti)
        out[:scan_name] = DimArray(
            ["3", "1", "3", "2"][mod1.(1:length(t), 4)], (t,);
            metadata = metadata(ms[:scan_name])
        )
        @test partition_info(out).scan_name == ["1", "2", "3"]
        # The summary sorts; `scans` keeps the order names first occur in.
        @test scans(out) == ["3", "1", "2"]
    end

    @testset "the data group decides which dataset is meant" begin
        out = deepcopy(ms)
        other = XRadio.Dataset(; metadata = Dict{Symbol, Any}(:type => "field_and_source"))
        other[:source_name] = DimArray(["zeta", "alpha"], (FieldName(Categorical(["f2", "f1"])),))
        branches(out)[:field_and_source_corrected] = other
        metadata(out)[:data_groups][:corrected] = merge(
            datagroup(out), Dict{Symbol, Any}(:field_and_source => "field_and_source_corrected_xds")
        )
        @test field_and_source(out, :corrected) === other
        info = partition_info(out, :corrected)
        @test info.source_name == ["alpha", "zeta"]
        @test info.field_name == ["f1", "f2"]
        @test info.data_group_name == :corrected
        @test partition_info(out).source_name == ["source_0"]
    end

    @testset "a group pointing at a dataset the set does not hold" begin
        out = deepcopy(ms)
        metadata(out)[:data_groups][:base][:field_and_source] = "field_and_source_gone_xds"
        @test_throws "points at `field_and_source_gone_xds`, which this Measurement Set does not hold" field_and_source(out)
        delete!(metadata(out)[:data_groups][:base], :field_and_source)
        @test_throws "names no `field_and_source` dataset" field_and_source(out)
        @test_throws KeyError field_and_source(ms, :nonexistent)
    end

    @testset "intents" begin
        out = deepcopy(ms)
        metadata(out)[:observation_info][:intents] = ["CALIBRATE_BANDPASS#ON_SOURCE"]
        @test partition_info(out).scan_intents == ["CALIBRATE_BANDPASS#ON_SOURCE"]

        out = deepcopy(ms)
        out[:scan_name] = DimArray(parent(ms[:scan_name]), dims(ms[:scan_name]))
        @test partition_info(out).scan_intents == ["UNSPECIFIED"]
        f = dims(out, XRadio.Frequency)
        delete!(metadata(lookup(f)), :spectral_window_intents)
        @test partition_info(out).spectral_window_intents == ["UNSPECIFIED"]
    end

    @testset "lines" begin
        out = deepcopy(ms)
        fs = field_and_source(out)
        fs[:line_name] = DimArray(
            ["SiO", "CO"], (LineLabel(Categorical(["l1", "l2"])),)
        )
        @test partition_info(out).line_name == ["CO", "SiO"]
    end
end

@testset "the antennas of one baseline" begin
    ms = synthetic_set(1)[1]

    ids = val(dims(ms, BaselineID))

    @test baseline_antennas(ms, first(ids)) == ("AA", "BB")
    @test baseline_antennas(ms, last(ids)) == ("BB", "CC")

    @testset "it inverts `baseline`" begin
        for (id, pair) in zip(ids, baselines(ms))
            @test baseline_antennas(ms, id) == pair
            @test baseline(ms, pair) == baseline(ms, id)
        end
    end

    @testset "an absent id names the range that is present" begin
        @test_throws ArgumentError baseline_antennas(ms, 7)
        @test_throws "no baseline_id 7" baseline_antennas(ms, 7)
        @test_throws "ids from 1 to 3" baseline_antennas(ms, 7)
    end

    @testset "a Measurement Set without the coordinates" begin
        bare = deepcopy(ms)
        delete!(bare, :baseline_antenna2_name)
        @test_throws "names no baseline antennas" baseline_antennas(bare, first(ids))
        # The `baseline_id` dimension going missing is the more basic fact and
        # is reported ahead of the coordinates.
        @test_throws "no `baseline_id` dimension" baseline_antennas(
            MeasurementSet(), 0
        )
    end
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

        # A `baseline_id` is the `baseline_id` dimension's own question, so it
        # is answered without the antenna-name coordinates.
        @test_throws "no baseline_id 9" baseline(ms, 9)
        @test_throws "3 baselines, with ids from 1 to 3" baseline(ms, 9)
        idonly = synthetic_set(1)[1]
        delete!(idonly, :baseline_antenna1_name)
        @test baseline(idonly, 2) == BaselineID(At(2))

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

    @testset "the answers follow the coordinates" begin
        # The accessors answer from a lookup built on the first query, so an
        # earlier answer must not survive a change to what it was built from.
        changed = synthetic_set(1)[1]
        @test baselines(changed) == [("AA", "BB"), ("AA", "CC"), ("BB", "CC")]

        base = dims(changed, BaselineID)
        changed[:baseline_antenna1_name] = DimArray(["AA", "AA", "ZZ"], (base,))
        @test baselines(changed) == [("AA", "BB"), ("AA", "CC"), ("ZZ", "CC")]
        @test baseline(changed, ("ZZ", "CC")) == BaselineID(At(3))
        @test_throws "no baseline BB-CC" baseline(changed, ("BB", "CC"))
        @test selected(baselines(changed, "AA")) == [1, 2]

        delete!(changed, :baseline_antenna1_name)
        @test_throws "names no baseline antennas" baselines(changed)

        popped = synthetic_set(1)[1]
        @test length(baselines(popped)) == 3
        pop!(popped, :baseline_antenna2_name)
        @test_throws "names no baseline antennas" baselines(popped)
    end

    @testset "a lookup answers as its Measurement Set does" begin
        bl = XRadio.baseline_lookup(ms)
        @test length(bl) == 3
        @test baselines(bl) == baselines(ms)
        @test baseline(bl, ("CC", "AA")) == baseline(ms, ("CC", "AA"))
        @test baseline(bl, "AA-CC") == baseline(ms, "AA-CC")
        @test selected(baselines(bl, "AA")) == selected(baselines(ms, "AA"))
        @test_throws "no baseline AA-ZZ" baseline(bl, ("AA", "ZZ"))

        # Holding the lookup is what makes the query concrete, which is the
        # whole reason it can be held.
        query(lookup, key) = @allocated baseline(lookup, key)
        query(bl, ("AA", "CC"))
        @test query(bl, ("AA", "CC")) == 0

        # It describes the Measurement Set it was built from, so a change to
        # that Measurement Set does not reach a lookup already taken.
        held = synthetic_set(1)[1]
        taken = XRadio.baseline_lookup(held)
        held[:baseline_antenna1_name] = DimArray(["ZZ", "AA", "BB"], (dims(held, BaselineID),))
        @test baselines(taken) == [("AA", "BB"), ("AA", "CC"), ("BB", "CC")]
        @test baselines(held) == [("ZZ", "BB"), ("AA", "CC"), ("BB", "CC")]
    end

    @testset "a repeated pair names its first baseline" begin
        doubled = synthetic_set(1)[1]
        base = dims(doubled, BaselineID)
        doubled[:baseline_antenna1_name] = DimArray(["AA", "AA", "BB"], (base,))
        doubled[:baseline_antenna2_name] = DimArray(["BB", "BB", "CC"], (base,))
        @test baseline(doubled, ("AA", "BB")) == BaselineID(At(1))
        @test selected(baselines(doubled, "AA")) == [1, 2]
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

@testset "scan and field selection" begin
    # A Measurement Set whose scans and fields both change part-way along the
    # time axis, so the two are not the same partition of it.
    function labelled(scannames, fieldnames)
        ms = synthetic_visibility(; ntime = length(scannames))
        t = dims(ms, Ti)
        ms[:scan_name] = DimArray(
            collect(scannames), (t,); metadata = metadata(ms[:scan_name])
        )
        ms[:field_name] = DimArray(collect(fieldnames), (t,))
        return ms
    end

    ms = labelled(["1", "1", "1", "2", "2", "3"], ["A", "A", "A", "A", "B", "B"])
    times = collect(val(dims(ms, Ti)))
    selected(sel) = DimensionalData.val(DimensionalData.val(sel))

    @test scans(ms) == ["1", "2", "3"]
    @test fields(ms) == ["A", "B"]

    @testset "a selector names the integrations the label covers" begin
        @test selected(scan(ms, "2")) == times[4:5]
        @test selected(field(ms, "B")) == times[5:6]
        @test size(correlated(ms)[scan(ms, "1")], Ti) == 3
        @test size(ms[field(ms, "A")][:visibility], Ti) == 4
        @test size(uvw(ms)[scan(ms, "3"), baseline(ms, ("AA", "BB"))]) == (3, 1)
    end

    @testset "a label that is not there names the ones that are" begin
        @test_throws "no scan `4`" scan(ms, "4")
        @test_throws "it names 1, 2, 3" scan(ms, "4")
        @test_throws "no field `C`" field(ms, "C")
    end

    @testset "the lookup is built once and answers as its Measurement Set does" begin
        held = labelled(["1", "1", "2", "2"], ["A", "A", "A", "A"])
        sl = XRadio.scan_lookup(held)
        @test sl === XRadio.scan_lookup(held)
        @test XRadio.labels(sl) == scans(held)
        @test length(sl) == 2

        t = collect(val(dims(held, Ti)))
        @test selected(scan(sl, t, "2")) == selected(scan(held, "2"))
        @test selected(field(XRadio.field_lookup(held), t, "A")) ==
            selected(field(held, "A"))

        # It describes the Measurement Set it was built from, so a change to
        # that Measurement Set does not reach a lookup already taken.
        held[:scan_name] = DimArray(
            ["9", "9", "9", "9"], (dims(held, Ti),); metadata = metadata(held[:scan_name])
        )
        @test XRadio.labels(sl) == ["1", "2"]
        @test scans(held) == ["9"]
    end

    @testset "the answers follow the coordinates" begin
        changed = labelled(["1", "1", "2", "2"], ["A", "A", "B", "B"])
        @test scans(changed) == ["1", "2"]

        t = dims(changed, Ti)
        changed[:scan_name] = DimArray(
            ["7", "7", "7", "8"], (t,); metadata = metadata(changed[:scan_name])
        )
        @test scans(changed) == ["7", "8"]
        @test selected(scan(changed, "7")) == collect(val(t))[1:3]
        @test_throws "no scan `1`" scan(changed, "1")

        delete!(changed, :scan_name)
        @test scans(changed) == String[]
        @test_throws "states no `scan_name`" scan(changed, "7")

        popped = labelled(["1", "2"], ["A", "A"])
        @test fields(popped) == ["A"]
        pop!(popped, :field_name)
        @test fields(popped) == String[]
    end

    @testset "a cut Measurement Set does not inherit the lookup" begin
        # Cutting rebuilds the tree, and the rebuilt one must not answer with
        # the labels of the set it came from: after a time cut those name the
        # wrong integrations, and nothing would say so.
        whole = labelled(["1", "1", "1", "2", "2", "3"], ["A", "A", "A", "A", "B", "B"])
        @test scans(whole) == ["1", "2", "3"]
        @test fields(whole) == ["A", "B"]

        t = collect(val(dims(whole, Ti)))
        cut = whole[Ti(At(t[4:6]))]
        @test getfield(cut, :scan_lookup) === nothing
        @test getfield(cut, :field_lookup) === nothing
        @test scans(cut) == ["2", "3"]
        @test selected(scan(cut, "2")) == t[4:5]
        @test_throws "no scan `1`" scan(cut, "1")

        # The set that was cut is left as it was.
        @test scans(whole) == ["1", "2", "3"]

        # Keeping whole members keeps their lookups, which nothing has
        # invalidated.
        ps = ProcessingSet(
            OrderedCollections.OrderedDict{Symbol, MeasurementSet}(:a => whole)
        )
        @test select(ps; scan = "1")[1] === whole
    end

    @testset "the answer is the caller's to keep" begin
        own = labelled(["1", "2"], ["A", "A"])
        got = scans(own)
        empty!(got)
        @test scans(own) == ["1", "2"]
    end

    @testset "a Measurement Set naming no scans" begin
        bare = MeasurementSet(; metadata = Dict{Symbol, Any}(:type => "visibility"))
        @test scans(bare) == String[]
        @test fields(bare) == String[]
        @test_throws "states no `scan_name`" XRadio.scan_lookup(bare)
        @test_throws "states no `field_name`" XRadio.field_lookup(bare)
    end

    @testset "a coordinate stated over the wrong dimension" begin
        # The positions a lookup holds index the time axis, so a coordinate
        # stated over anything else would name the wrong integrations.
        wrong = synthetic_visibility()
        wrong[:scan_name] = DimArray(["a", "b", "c"], (dims(wrong, BaselineID),))
        @test_throws "not over `time` alone" scans(wrong)
    end
end

@testset "selecting Measurement Sets by attribute" begin
    function varied(; field, spw, scan, ants, intents = ["OBSERVE_TARGET#ON_SOURCE"])
        ms = synthetic_visibility()
        t = dims(ms, Ti)
        ms[:field_name] = DimArray(fill(field, length(t)), (t,))
        ms[:scan_name] = DimArray(
            fill(scan, length(t)), (t,);
            metadata = Dict{Symbol, Any}(:scan_intents => intents)
        )
        metadata(lookup(dims(ms, Frequency)))[:spectral_window_name] = spw
        branches(ms)[:antenna] = synthetic_antenna(ants)
        return ms
    end

    ps = ProcessingSet(
        OrderedCollections.OrderedDict{Symbol, MeasurementSet}(
            :a => varied(field = "M87", spw = "spw_0", scan = "1", ants = ["AA", "BB"]),
            :b => varied(
                field = "M87", spw = "spw_1", scan = "2", ants = ["AA", "CC"],
                intents = ["CALIBRATE_BANDPASS#ON_SOURCE"],
            ),
            :c => varied(
                field = "3C273", spw = "spw_0", scan = "3", ants = ["BB", "CC"],
                intents = ["CALIBRATE_BANDPASS"],
            ),
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

    @testset "by intent" begin
        # A wanted intent naming no subintent asks only what the scan was for,
        # so it matches whichever subintent the store states.
        @test collect(keys(select(ps; intent = "CALIBRATE_BANDPASS"))) == [:b, :c]
        # One that names a subintent asks for both halves.
        @test collect(keys(select(ps; intent = "CALIBRATE_BANDPASS#ON_SOURCE"))) == [:b]
        @test collect(keys(select(ps; intent = "OBSERVE_TARGET"))) == [:a]

        # The enumeration member and the parsed intent say the same thing as
        # the string does.
        @test keys(select(ps; intent = ScanIntents.CALIBRATE_BANDPASS)) ==
            keys(select(ps; intent = "CALIBRATE_BANDPASS"))
        @test keys(select(ps; intent = Intent(ScanIntents.OBSERVE_TARGET))) ==
            keys(select(ps; intent = "OBSERVE_TARGET"))

        @test collect(keys(select(
            ps; intent = ["OBSERVE_TARGET", "CALIBRATE_BANDPASS#ON_SOURCE"]
        ))) == [:a, :b]
        @test isempty(select(ps; intent = "CALIBRATE_POINTING"))
        # `intent` narrows alongside the other keywords rather than replacing
        # them.
        @test collect(keys(select(ps; intent = "CALIBRATE_BANDPASS", field = "M87"))) ==
            [:b]

        # An intent outside the ASDM enumerations is an error, not a name that
        # quietly matches nothing.
        @test_throws ArgumentError select(ps; intent = "NOT_AN_INTENT")
        @test_throws "is not an ASDM scan intent" select(ps; intent = "NOT_AN_INTENT")
    end

    @testset "an unknown keyword is an error" begin
        # The whole risk of a keyword-driven selector is a typo silently
        # selecting everything.
        @test_throws ArgumentError select(ps; spectral_window = "spw_0")
        @test_throws "no keyword `spectral_window`" select(ps; spectral_window = "spw_0")
        @test_throws "accepts antenna, field, intent, scan, spectralwindow" select(
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

@testset "cutting inside one Measurement Set" begin
    ms = synthetic_set(1)[1]
    t = val(dims(ms, Ti))

    # What `select`'s docstring points at for cutting a Measurement Set, as
    # against choosing whole ones out of a processing set.
    @test length(dims(ms[baselines(ms, "AA")], BaselineID)) == 2
    @test length(dims(ms[Ti(t[1] .. t[2])], Ti)) == 2
    # A scalar selector drops the dimension, on the tree as on a layer.
    @test !hasdim(correlated(ms)[baseline(ms, ("AA", "CC"))], BaselineID)
    @test !hasdim(ms[baseline(ms, ("AA", "CC"))], BaselineID)

    # A one-element collection keeps it instead.
    id = first(val(dims(ms, BaselineID)))
    @test length(dims(ms[BaselineID(At([id]))], BaselineID)) == 1

    @testset "a coordinate left over the dropped dimension holds one value" begin
        # MSv4 keeps such a coordinate as a scalar, and that scalar is the
        # answer to what was selected — the dimension goes, the value does not.
        cut = ms[baseline(ms, ("AA", "CC"))]
        @test ndims(cut[:baseline_antenna1_name]) == 0
        @test cut[:baseline_antenna1_name][] == "AA"
        @test cut[:baseline_antenna2_name][] == "CC"
        # A layer spanning more than the dropped dimension keeps the rest.
        @test size(cut[:visibility]) == (2, 5, size(correlated(ms), Ti))

        # Attributes survive the reduction, which the scan intents need.
        tcut = ms[Ti(At(t[1]))]
        @test ndims(tcut[:scan_name]) == 0
        @test scans(tcut) == ["1"]
        @test scan_intents(tcut) ==
            [Intent(ScanIntents.OBSERVE_TARGET, SubscanIntents.ON_SOURCE)]
    end

    @testset "a dimension is named by type or by the store's own name" begin
        # `DimensionalData` resolves a keyword through `name2dim`, which would
        # label `baseline_id` as `Dim{:baseline_id}` and match nothing; MSv4
        # names dimensions as the store does, and both spellings must cut.
        @test parent(ms[baseline_id = At([id])][:visibility]) ==
            parent(ms[BaselineID(At([id]))][:visibility])
        @test parent(ms[time = At(t[1])][:visibility]) ==
            parent(ms[Ti(At(t[1]))][:visibility])
        @test !hasdim(ms[baseline_id = At(id)], BaselineID)
        @test length(dims(ms[frequency = At(frequencies(ms)[1:2])], Frequency)) == 2
        @test !hasdim(ms[polarization = At("XX")], Polarization)

        # A dimension outside the standard has no type of its own, so the
        # keyword is its only spelling.
        own = deepcopy(ms)
        own[:solution] = DimArray([1.0, 2.0], (Dim{:solution_interval}([0.0, 1.0]),))
        @test own[solution_interval = At(1.0)][:solution][] == 2.0
    end

    @testset "a dimension nothing carries is an error" begin
        @test_throws "no `nosuchdim` dimension" ms[nosuchdim = At(3)]
        @test_throws "no `nosuchdim` dimension" ms[Dim{:nosuchdim}(At(3))]
        # The message names what is there, as the store spells it.
        @test_throws "`baseline_id`" ms[nosuchdim = At(3)]

        # A dimension only a sub-dataset carries is not absent.
        @test hasdim(branches(ms)[:antenna], AntennaName)
        @test !hasdim(ms, AntennaName)
        cut = ms[antenna_name = At(["AA"])]
        @test length(dims(branches(cut)[:antenna], AntennaName)) == 1
        # and cutting on it leaves the correlated data alone.
        @test size(correlated(cut)) == size(correlated(ms))
    end

    @testset "selection does not warn about nodes that lack the dimension" begin
        # Every dimension used to be applied to every node, so the `antenna`
        # branch warned for having no `baseline_id` — and the warning that
        # means something was lost among those that do not.
        @test_nowarn ms[BaselineID(At([id]))]
        @test_nowarn ms[Ti(t[1] .. t[2])]
        @test_nowarn ms[baseline_id = At([id])]
    end
end

@testset "merging processing sets" begin
    ps = synthetic_set(4)
    first_two = filter(ms -> true, ProcessingSet(
        OrderedCollections.OrderedDict(k => ps[k] for k in (:test_00, :test_01))
    ))
    last_two = filter(ms -> true, ProcessingSet(
        OrderedCollections.OrderedDict(k => ps[k] for k in (:test_02, :test_03))
    ))

    got = merge(first_two, last_two)
    @test got isa ProcessingSet
    @test collect(keys(got)) == [:test_00, :test_01, :test_02, :test_03]
    @test metadata(got)[:type] == "processing_set"
    # The Measurement Sets are shared, not copied.
    @test got[:test_00] === ps[:test_00]

    @test keys(merge(first_two)) == keys(first_two)
    @test collect(keys(merge(last_two, first_two))) ==
        [:test_02, :test_03, :test_00, :test_01]

    @testset "a shared name is an error" begin
        # Two Measurement Sets under one name are two different observations,
        # so neither may win silently.
        @test_throws ArgumentError merge(ps, first_two)
        @test_throws "names a Measurement Set `test_00`" merge(ps, first_two)
        @test_throws ArgumentError merge(first_two, last_two, first_two)
    end

    @testset "attributes are the union, and must agree" begin
        # A processing set's attributes are its store's root attributes, which
        # `write` puts back verbatim, so one argument's cannot stand for all.
        metadata(first_two)[:observer] = "someone"
        metadata(last_two)[:telescope] = "TEST"
        got = merge(first_two, last_two)
        @test got isa ProcessingSet
        @test metadata(got)[:observer] == "someone"
        @test metadata(got)[:telescope] == "TEST"
        @test metadata(got)[:type] == "processing_set"
        # The arguments are left as they were.
        @test !haskey(metadata(first_two), :telescope)

        metadata(last_two)[:observer] = "someone else"
        @test_throws ArgumentError merge(first_two, last_two)
        @test_throws "disagree about `observer`" merge(first_two, last_two)
        @test_throws "\"someone\" against \"someone else\"" merge(first_two, last_two)
    end

    # Merging the parts of a selection puts the set back together.
    @test keys(merge(
        select(ps; scan = "1"), filter(ms -> false, ps)
    )) == keys(ps)
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
