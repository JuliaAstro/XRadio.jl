using XRadio
using Test

using DimensionalData
using LinearAlgebra: dot
using DimensionalData: branches
using FITSFiles
using OrderedCollections
using Zarr

using XRadio: check, validate, SchemaError, SchemaIssue
using XRadio: msv4name, isvalid

include("reference_schema.jl")
include("synthetic.jl")
include("synthetic_fitsidi.jl")

@testset "XRadio" begin
    include("test_api.jl")
    include("test_schema.jl")
    include("test_io.jl")
    include("test_check.jl")
    include("test_build.jl")
    include("test_average.jl")
    include("test_mount.jl")
    include("test_gaincurve.jl")
    include("test_fitsidi.jl")
end
