using XRadio
using Test

using DimensionalData
using DimensionalData: branches
using FITSFiles
using OrderedCollections
using Zarr

using XRadio: check, validate, SchemaError, SchemaIssue
using XRadio: msv4name, isvalid

include("reference_schema.jl")
include("synthetic.jl")

@testset "XRadio" begin
    include("test_schema.jl")
    include("test_io.jl")
    include("test_check.jl")
    include("test_fitsidi.jl")
end
