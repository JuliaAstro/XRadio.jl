@testset "the checker catches" begin
    @testset "a missing required variable" begin
        ps = synthetic_set(1)
        delete!(DimensionalData.data(ps[1]), :weight)
        delete!(DimensionalData.layerdims(ps[1]), :weight)
        delete!(DimensionalData.layermetadata(ps[1]), :weight)
        issues = check(ps)
        @test any(i -> occursin("WEIGHT", i.path) && occursin("required", i.message), issues)
    end

    @testset "a missing optional variable" begin
        # Optional variables are absent from the synthetic data, and that is
        # not a fault.
        ps = synthetic_set(1)
        @test !haskey(ps[1], :time_centroid)
        @test isempty(check(ps))
    end

    @testset "an unknown dataset type" begin
        ps = synthetic_set(1)
        metadata(ps[1])[:type] = "not_a_dataset"
        @test any(i -> occursin("names no MSv4 dataset", i.message), check(ps))
    end

    @testset "a missing attribute" begin
        ps = synthetic_set(1)
        delete!(metadata(ps[1]), :creator)
        @test any(i -> occursin("`creator` is required", i.message), check(ps))
    end

    @testset "an attribute outside its permitted values" begin
        ps = synthetic_set(1)
        metadata(lookup(ps[1][:visibility], Ti))[:scale] = "galactic"
        issues = check(ps)
        @test any(i -> occursin("`scale` is \"galactic\"", i.message), issues)
    end

    @testset "a variable the schema does not define" begin
        ps = synthetic_set(1)
        ps[1][:nonsense] = ps[1][:weight]
        @test any(i -> occursin("not defined by the `visibility` schema", i.message), check(ps))
    end

    @testset "an extra version of a science array" begin
        # Data groups exist so a dataset can hold several versions of its
        # science arrays, so these must not be reported as unknown.
        ps = synthetic_set(1)
        ps[1][:visibility_corrected] = ps[1][:visibility]
        @test isempty(check(ps))
    end

    @testset "dimensions in the wrong order" begin
        ps = synthetic_set(1)
        vis = ps[1][:visibility]
        ps[1][:visibility] = permutedims(vis, (2, 1, 3, 4))
        @test any(i -> occursin("spans", i.message), check(ps))
    end

    @testset "a data group naming a variable that is not there" begin
        ps = synthetic_set(1)
        datagroup(ps[1])[:weight] = "WEIGHT_MISSING"
        issues = check(ps)
        @test any(i -> occursin("which the dataset does not have", i.message), issues)
    end

    @testset "a data group with no base" begin
        ps = synthetic_set(1)
        groups = datagroups(ps[1])
        groups[:other] = pop!(groups, :base)
        @test any(i -> occursin("no `base` group", i.message), check(ps))
    end

    @testset "no antenna sub-dataset" begin
        ps = synthetic_set(1)
        delete!(branches(ps[1]), :antenna)
        @test any(i -> occursin("no `antenna_xds`", i.message), check(ps))
    end

    @testset "a processing set with the wrong type" begin
        ps = synthetic_set(1)
        metadata(ps)[:type] = "something_else"
        @test any(i -> occursin("expected \"processing_set\"", i.message), check(ps))
    end
end

@testset "validate" begin
    ps = synthetic_set(1)
    @test validate(ps) === ps

    delete!(metadata(ps[1]), :creator)
    @test_throws SchemaError validate(ps)

    err = try
        validate(ps)
    catch e
        e
    end
    # The error names every departure, not just the first.
    @test occursin("MSv4 standard", sprint(showerror, err))
    @test occursin("creator", sprint(showerror, err))

    # Writing validates by default, so a faulty set cannot be written by
    # accident — but can be, deliberately.
    path = joinpath(mktempdir(), "bad.ps.zarr")
    @test_throws SchemaError write(path, ps)
    @test write(path, ps; validate = false, overwrite = true) == path
end
