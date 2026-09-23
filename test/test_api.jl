@testset "public API surface" begin
    exported = Set(filter(n -> Base.isexported(XRadio, n), names(XRadio)))

    @testset "exports resolve" begin
        for name in exported
            @test isdefined(XRadio, name)
            @test !startswith(string(name), "_")
        end
    end

    @testset "every MSv4 dimension is exported" begin
        # A dimension type a caller cannot name is a dimension they cannot
        # select along, so the vocabulary and the export list move together.
        # `time` is DimensionalData's own `Ti`, which arrives by re-export.
        for D in XRadio.MSV4_DIMS
            @test nameof(D) in exported
        end
    end

    @testset "acronyms are capitalized" begin
        # `BaselineID`, `ElementID`, `UVWLabel` — one spelling rule for all of
        # them, so a caller need not remember which acronym got which case.
        @test BaselineID === XRadio.dimtype(:baseline_id)
        @test ElementID === XRadio.dimtype(:element_id)
        @test UVWLabel === XRadio.dimtype(:uvw_label)
        for D in XRadio.MSV4_DIMS
            @test !occursin(r"Id\b|Uvw", string(nameof(D)))
        end
    end

    @testset "public but not exported" begin
        # Conformance checking and the gain-curve interface are part of what
        # the package promises without being what a caller types.
        public_only = (
            :check, :validate, :SchemaIssue, :SchemaError,
            :gain_variable, :gain_quantity, :msv4_gain_curve,
            :msv4_mount, :mount_from_msv4, :axis_offset, :set_mounts!,
        )
        for name in public_only
            @test isdefined(XRadio, name)
            @test !(name in exported)
        end
        @static if VERSION >= v"1.11"
            for name in public_only
                @test Base.ispublic(XRadio, name)
            end
        end
    end

    @testset "converter entry points" begin
        # File in, store out: the conversion streams to disk, so it is its own
        # verb rather than a `read` that would have to hold the result.
        @test :fitsidi2msv4 in exported
        @test :fitsidi_gain_curves in exported
        @test hasmethod(fitsidi2msv4, Tuple{AbstractString, AbstractString})
        @test hasmethod(fitsidi_gain_curves, Tuple{AbstractString})
    end
end
