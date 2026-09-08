module XRadioFITSFilesExt

using Dates
using FITSFiles

using XRadio
using DimensionalData: DimStack, Dim, branches
using DimensionalData.Lookups: Categorical, Sampled

include("fitsidi.jl")
include("convert.jl")

end # module
