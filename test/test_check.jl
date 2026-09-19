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
        @test any(i -> occursin("names no known dataset", i.message), check(ps))
    end

    @testset "a missing attribute" begin
        ps = synthetic_set(1)
        delete!(metadata(ps[1]), :creator)
        @test any(i -> occursin("`creator` is required", i.message), check(ps))
    end

    @testset "an attribute present but empty" begin
        # An attribute whose value is nothing is missing, not wrong.
        ps = synthetic_set(1)
        metadata(ps[1])[:creator] = nothing
        @test any(i -> occursin("`creator` is required", i.message), check(ps))
    end

    @testset "an attribute of the wrong type" begin
        ps = synthetic_set(1)
        metadata(ps[1])[:creation_date] = 2026
        @test any(i -> occursin("`creation_date` is a Int", i.message), check(ps))

        ps = synthetic_set(1)
        metadata(ps[1])[:observation_info][:observer] = "Nobody"
        @test any(i -> occursin("expected a list of strings", i.message), check(ps))

        ps = synthetic_set(1)
        metadata(ps[1])[:processor_info] = "CORRELATOR"
        @test any(i -> occursin("expected a dictionary", i.message), check(ps))
    end

    @testset "an attribute outside its permitted values" begin
        ps = synthetic_set(1)
        metadata(lookup(ps[1][:visibility], Ti))[:scale] = "galactic"
        issues = check(ps)
        @test any(i -> occursin("`scale` is \"galactic\"", i.message), issues)
    end

    @testset "a restricted attribute holding something that is not a string" begin
        ps = synthetic_set(1)
        metadata(lookup(ps[1][:visibility], Ti))[:scale] = 3
        @test any(i -> occursin("`scale`", i.message), check(ps))
    end

    @testset "a fault inside a nested metadata dictionary" begin
        ps = synthetic_set(1)
        delete!(metadata(ps[1])[:observation_info], :project_UID)
        issues = check(ps)
        @test any(
            i -> occursin("observation_info", i.path) && occursin("`project_UID`", i.message),
            issues
        )
    end

    @testset "a fault inside an array-valued attribute" begin
        ps = synthetic_set(1)
        meta = metadata(lookup(ps[1][:visibility], Ti))
        meta[:integration_time] = XRadio.Measure(10, Dict{Symbol, Any}(:units => "s"))
        issues = check(ps)
        # Its element type is wrong, and it is missing its `type` attribute.
        @test any(i -> occursin("integration_time", i.path) && occursin("element type", i.message), issues)
        @test any(i -> occursin("integration_time", i.path) && occursin("`type`", i.message), issues)
    end

    @testset "a data group with no base" begin
        ps = synthetic_set(1)
        groups = datagroups(ps[1])
        groups[:other] = pop!(groups, :base)
        @test any(i -> occursin("`base` is required", i.message), check(ps))
    end

    @testset "a data group missing a field of its own" begin
        ps = synthetic_set(1)
        delete!(datagroup(ps[1]), :description)
        @test any(i -> occursin("`description` is required", i.message), check(ps))
    end

    @testset "dimensions in the wrong order" begin
        ps = synthetic_set(1)
        vis = ps[1][:visibility]
        ps[1][:visibility] = permutedims(vis, (2, 1, 3, 4))
        @test any(i -> occursin("spans", i.message), check(ps))
    end

    @testset "a dimension the schema does not define" begin
        ps = synthetic_set(1)
        ps[1][:nonsense] = DimArray(zeros(2), (Dim{:nonsense_label}(1:2),))
        @test any(i -> occursin("`nonsense_label`", i.message), check(ps))
    end

    @testset "an array carrying no attributes at all" begin
        # Metadata that is not a dictionary has every attribute missing, rather
        # than passing the whole schema unexamined.
        ps = synthetic_set(1)
        ant = branches(ps[1])[:antenna]
        pos = ant[:antenna_position]
        ant[:antenna_position] = DimArray(parent(pos), dims(pos))
        issues = check(ps)
        @test any(
            i -> occursin("ANTENNA_POSITION", i.path) && occursin("`units` is required", i.message),
            issues
        )
    end

    @testset "a dimension with no index coordinate" begin
        # A dimension a store gave no values is numbered instead, which MSv4
        # counts as a coordinate that is not there rather than as one present.
        ds = Dataset(; metadata = Dict{Symbol, Any}(:type => "field_and_source"))
        ds[:source_name] = DimArray(["source_0"], (FieldName(Lookups.NoLookup(Base.OneTo(1))),))
        issues = check(ds)
        @test any(i -> occursin("field_name", i.path) && occursin("required", i.message), issues)
    end

    @testset "a fault in a sub-dataset nested more than one level deep" begin
        ps = synthetic_set(1)
        inner = branches(ps[1])[:field_and_source_base]
        branches(inner)[:extra] = synthetic_field_and_source()
        delete!(metadata(branches(inner)[:extra]), :type)
        @test any(i -> occursin("extra", i.path), check(ps))
    end

    @testset "a processing set with the wrong type" begin
        ps = synthetic_set(1)
        metadata(ps)[:type] = "something_else"
        @test any(i -> occursin("expected \"processing_set\"", i.message), check(ps))
    end
end

@testset "the checker accepts what the reference accepts" begin
    @testset "a variable the schema does not define" begin
        # Extra data variables and coordinates are always allowed, as long as
        # they span dimensions the dataset already has.
        ps = synthetic_set(1)
        ps[1][:nonsense] = ps[1][:weight]
        @test isempty(check(ps))
    end

    @testset "an attribute the schema does not define" begin
        ps = synthetic_set(1)
        metadata(ps[1])[:nonsense] = "whatever"
        @test isempty(check(ps))
    end

    @testset "an extra version of a science array" begin
        # Data groups exist so a dataset can hold several versions of its
        # science arrays, so these must not be reported as unknown — and each
        # version is checked against the specification it is a version of.
        ps = synthetic_set(1)
        ps[1][:visibility_corrected] = ps[1][:visibility]
        @test isempty(check(ps))

        ps[1][:visibility_corrected] = permutedims(ps[1][:visibility], (2, 1, 3, 4))
        @test any(i -> occursin("VISIBILITY_CORRECTED", i.path), check(ps))
    end

    @testset "no antenna sub-dataset" begin
        # MSv4 declares no sub-dataset as a required field of any dataset, so
        # neither does the checker.
        ps = synthetic_set(1)
        delete!(branches(ps[1]), :antenna)
        @test isempty(check(ps))
    end

    @testset "a data group naming a variable that is not there" begin
        # A data group's fields are checked for presence and type; the names
        # they hold are not resolved against the dataset.
        ps = synthetic_set(1)
        datagroup(ps[1])[:weight] = "WEIGHT_MISSING"
        @test isempty(check(ps))
    end

    @testset "a branch holding no data" begin
        # A node with neither variables nor dimensions is a branch of the tree
        # rather than a dataset, and has nothing to conform to.
        ps = synthetic_set(1)
        branches(ps[1])[:empty] = Dataset()
        @test isempty(check(ps))
    end

    @testset "a processing set holding nothing" begin
        @test isempty(check(synthetic_set(0)))
    end

    @testset "a data group named anything else" begin
        ps = synthetic_set(1)
        datagroups(ps[1])[:imaging] = Dict{Symbol, Any}(:whatever => 1)
        @test isempty(check(ps))
    end
end

@testset "a dataset type of one's own" begin
    # MSv4 extends at the level of whole sub-datasets: a node whose type the
    # standard does not define conforms once the caller says what it must look
    # like. Nothing is registered globally, so the same store checks the same
    # way whatever else has run.
    spec = XRadio.DatasetSpec(
        :calibration, "calibration";
        coords = [XRadio.ArraySpec("antenna_name", (:antenna_name,), String)],
        vars = [XRadio.ArraySpec("GAIN", (:antenna_name,), ComplexF64)],
        attrs = [XRadio.AttrSpec(:type, :str; values = ["calibration"])],
    )

    ps = synthetic_set(1)
    cal = Dataset(; metadata = Dict{Symbol, Any}(:type => "calibration"))
    cal[:gain] = DimArray(
        randn(ComplexF64, 3), (AntennaName(Lookups.Categorical(["AA", "BB", "CC"])),)
    )
    branches(ps[1])[:calibration] = cal

    @test any(i -> occursin("names no known dataset", i.message), check(ps))
    @test isempty(check(ps; schemas = [spec]))

    # The schema also decides what each layer is called on disk, so a store
    # written with it reads back under the schema's names.
    path = joinpath(mktempdir(), "extended.ps.zarr")
    write(path, ps; schemas = [spec])
    back = read(ProcessingSet, path)
    @test haskey(branches(back[1]), :calibration)
    @test haskey(branches(back[1])[:calibration], :gain)
    @test isempty(check(back; schemas = [spec]))
end

@testset "extending a standard schema" begin
    # `extend` is how a caller has its own entries checked rather than merely
    # accepted. The worked example is what Gustavo needs of an MSv4 store:
    # values radio data carries that the standard has no field for.
    eop = Dict{Symbol, Any}(
        :gst_iat0 => 100.0, :earth_rot_rate => 360.985, :ut1utc => -0.1,
        :polarx => 0.2, :polary => 0.3, :datutc => 37.0,
        :xyzhand => "RIGHT", :poltype => "APPROX",
    )
    orientation = XRadio.DictSpec(
        :earth_orientation,
        [
            XRadio.AttrSpec(f, :float) for f in
                (:gst_iat0, :earth_rot_rate, :ut1utc, :polarx, :polary, :datutc)
        ]
    )
    spec = XRadio.extend(
        XRadio.VISIBILITY_SCHEMA;
        attrs = [
            XRadio.AttrSpec(
                :earth_orientation, :dict; optional = true, nested = orientation
            ),
        ],
        coords = [
            XRadio.ArraySpec(
                "frequency", (:frequency,), Float64;
                attrs = [XRadio.AttrSpec(:sideband, :int; optional = true)]
            ),
            XRadio.ArraySpec("sub_scan_name", (:time,), String; optional = true),
        ],
    )

    @testset "the standard schema is untouched" begin
        @test !haskey(XRadio.VISIBILITY_SCHEMA, :sub_scan_name)
        @test !any(a -> a.name === :earth_orientation, XRadio.VISIBILITY_SCHEMA.attrs)
        @test !any(
            a -> a.name === :sideband, XRadio.VISIBILITY_SCHEMA[:frequency].attrs
        )
    end

    @testset "an entry already defined keeps everything but its attributes" begin
        freq = spec[:frequency]
        std = XRadio.VISIBILITY_SCHEMA[:frequency]
        @test freq.dims == std.dims
        @test freq.eltypes == std.eltypes
        @test any(a -> a.name === :sideband, freq.attrs)
        # The standard's own attributes survive, so an extension cannot
        # silently stop `spectral_window_name` from being required.
        @test any(a -> a.name === :spectral_window_name, freq.attrs)
    end

    @testset "restating a standard entry throws" begin
        @test_throws "may add attributes to it, but not restate it" XRadio.extend(
            XRadio.VISIBILITY_SCHEMA;
            coords = [XRadio.ArraySpec("frequency", (:time,), Float64)]
        )
    end

    @testset "the store is accepted either way, and checked only with the schema" begin
        ps = synthetic_set(1)
        ntime = length(dims(ps[1][:visibility], Ti))
        metadata(ps[1])[:earth_orientation] = copy(eop)
        metadata(lookup(dims(ps[1][:visibility], Frequency)))[:sideband] = 1
        ps[1][:sub_scan_name] = DimArray(
            fill("a", ntime), (dims(ps[1][:visibility], Ti),)
        )
        @test isempty(check(ps))
        @test isempty(check(ps; schemas = [spec]))

        # Wrong, and reported only by the schema that knows what it should be.
        # The standard check stays exactly as permissive as the reference.
        delete!(metadata(ps[1])[:earth_orientation], :ut1utc)
        metadata(lookup(dims(ps[1][:visibility], Frequency)))[:sideband] = "USB"
        @test isempty(check(ps))
        issues = check(ps; schemas = [spec])
        @test any(
            i -> occursin("earth_orientation", i.path) && occursin("`ut1utc`", i.message),
            issues
        )
        @test any(i -> occursin("`sideband` is a String", i.message), issues)
    end

    @testset "an added coordinate is written under its own name" begin
        ps = synthetic_set(1)
        ntime = length(dims(ps[1][:visibility], Ti))
        ps[1][:sub_scan_name] = DimArray(
            fill("a", ntime), (dims(ps[1][:visibility], Ti),)
        )

        path = joinpath(mktempdir(), "extended.ps.zarr")
        write(path, ps; schemas = [spec])
        @test isfile(joinpath(path, "test_00", "sub_scan_name", ".zarray")) ||
            isdir(joinpath(path, "test_00", "sub_scan_name"))
        back = read(ProcessingSet, path)
        @test back[1][:sub_scan_name] == ps[1][:sub_scan_name]

        # Without the schema, the same layer is written as a data variable.
        plain = joinpath(mktempdir(), "plain.ps.zarr")
        write(plain, ps)
        @test isdir(joinpath(plain, "test_00", "SUB_SCAN_NAME"))
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
