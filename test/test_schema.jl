@testset "schema matches the reference" begin
    refs = reference_datasets()

    # Every dataset type the reference defines must be one we implement, and
    # vice versa. `DopplerXds` is exported empty by the reference and carries
    # no type attribute, so it contributes nothing here.
    @testset "dataset coverage" begin
        for (type, _) in refs
            @test XRadio.schema(type) isa XRadio.DatasetSpec
        end
        for spec in values(XRadio.SCHEMAS), type in spec.type
            @test haskey(refs, type)
        end
    end

    for (type, json) in sort(collect(refs); by = first)
        spec = XRadio.schema(type)
        ref_arrays = reference_arrays(json)

        @testset "$type" begin
            @testset "arrays present" begin
                ours = Set(a.diskname for a in XRadio.arrays(spec))
                theirs = Set(keys(ref_arrays))
                @test setdiff(theirs, ours) == Set{String}()
                @test setdiff(ours, theirs) == Set{String}()
            end

            for a in XRadio.arrays(spec)
                haskey(ref_arrays, a.diskname) || continue
                ref = ref_arrays[a.diskname]
                @testset "$(a.diskname)" begin
                    @test a.optional == Bool(get(ref, :optional, false))
                    @test a.eltypes == reference_eltype(ref)
                    @test spec_dims(a) == reference_dims(ref)

                    ref_at = reference_attrs(ref)
                    ours_at = Dict(x.name => spec_attr(x) for x in a.attrs)
                    @test setdiff(keys(ref_at), keys(ours_at)) == Set{Symbol}()
                    @test setdiff(keys(ours_at), keys(ref_at)) == Set{Symbol}()
                    for (name, r) in ref_at
                        haskey(ours_at, name) || continue
                        @test ours_at[name] == r
                    end
                end
            end

            @testset "dataset attributes" begin
                ref_at = reference_attrs(json)
                ours_at = Dict(x.name => spec_attr(x) for x in spec.attrs)
                @test setdiff(keys(ref_at), keys(ours_at)) == Set{Symbol}()
                @test setdiff(keys(ours_at), keys(ref_at)) == Set{Symbol}()
                for (name, r) in ref_at
                    haskey(ours_at, name) || continue
                    @test ours_at[name] == r
                end
            end
        end
    end
end

@testset "dimension vocabulary" begin
    # Every dimension any schema names must be in the vocabulary, and every
    # dimension in the vocabulary must be used by some schema.
    used = Set{Symbol}()
    for spec in values(XRadio.SCHEMAS)
        union!(used, XRadio.dimensions(spec))
    end
    @test setdiff(used, Set(keys(XRadio.MSV4_DIMS))) == Set{Symbol}()
    @test setdiff(Set(keys(XRadio.MSV4_DIMS)), used) == Set{Symbol}()

    # `dimtype` and `msv4name` are inverse.
    for name in keys(XRadio.MSV4_DIMS)
        @test XRadio.msv4name(XRadio.dimtype(name)) === name
    end
    @test_throws ArgumentError XRadio.dimtype(:not_a_dimension)
end

@testset "relaxed element types" begin
    # The standard's dtypes widen to their kind, so an in-memory dataset is
    # free to use any float or complex representation.
    # `==` rather than `===`: two separately built `UnionAll`s are equal types
    # but not the same object.
    @test XRadio.relaxed_eltype(Union{Float16, Float32, Float64}) == AbstractFloat
    @test XRadio.relaxed_eltype(Union{ComplexF32, ComplexF64}) == Complex{<:AbstractFloat}
    @test XRadio.relaxed_eltype(Union{Int32, Int64}) == Integer
    @test XRadio.relaxed_eltype(String) == AbstractString
    @test Float16 <: XRadio.relaxed_eltype(Float64)
    # Flags mix booleans with unsigned bit fields, and must keep both.
    flags = XRadio.relaxed_eltype(Union{Bool, UInt8, UInt16, UInt32, UInt64})
    @test Bool <: flags
    @test UInt8 <: flags
    @test !(Int8 <: flags)
end

@testset "node names" begin
    @test XRadio.nodekey("antenna_xds") === :antenna
    @test XRadio.nodekey("field_and_source_base_xds") === :field_and_source_base
    @test XRadio.nodename(:antenna) == "antenna_xds"
    # A node without the suffix keeps its name rather than losing characters.
    @test XRadio.nodekey("antenna") === :antenna
end

@testset "julian names" begin
    @test XRadio.julianame("VISIBILITY") === :visibility
    @test XRadio.julianame("baseline_id") === :baseline_id
    @test XRadio.VISIBILITY_SCHEMA[:visibility].diskname == "VISIBILITY"
    @test haskey(XRadio.VISIBILITY_SCHEMA, :uvw)
    @test !haskey(XRadio.SPECTRUM_SCHEMA, :uvw)
end
