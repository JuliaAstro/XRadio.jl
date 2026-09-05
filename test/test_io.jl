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

@testset "a written store conforms" begin
    # Zarr v2 is the format this package can write fully to the standard; see
    # the `zarr_format` keyword of `write`.
    ps = synthetic_set(1)
    path = joinpath(mktempdir(), "conforming.ps.zarr")
    write(path, ps; zarr_format = 2)
    @test isempty(check(read(ProcessingSet, path)))
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
    @test length(filter(ms -> "OBSERVE_TARGET#ON_SOURCE" in intents(ms), ps)) == 3
    @test length(filter(ms -> isempty(intents(ms)), ps)) == 0
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

@testset "accessors" begin
    ms = synthetic_set(1)[1]
    @test antennas(ms) == ["AA", "BB", "CC"]
    @test scans(ms) == ["1"]
    @test fields(ms) == ["field_0"]
    @test intents(ms) == ["OBSERVE_TARGET#ON_SOURCE"]
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
