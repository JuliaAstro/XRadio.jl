module XRadioFITSFilesExt

using Dates
using FITSFiles

using XRadio
using DimensionalData: DimStack, Dim, metadata
using DimensionalData.Lookups: Categorical, Sampled

include("fitsidi.jl")
include("convert.jl")

end # module
