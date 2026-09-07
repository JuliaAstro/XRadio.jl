"""
    fitsidi2msv4(path, dest; kwargs...) -> dest

Convert the FITS-IDI file at `path` into an MSv4 processing set written to
`dest`, and return `dest`.

FITS-IDI (AIPS Memo 114) is the interchange format VLBI correlators write, and
this is the whole of XRadio's interface to it: the FITS-IDI tables, and the
judgement needed to map them onto MSv4, stay inside the extension rather than
becoming part of the package's vocabulary.

The conversion writes each partition as it is built rather than assembling a
[`ProcessingSet`](@ref) and writing it afterwards, so peak memory is one
partition regardless of how large the input is. Correlator output routinely
runs to tens of gigabytes, which is more than a machine can be assumed to
hold.

Requires `FITSFiles`:

```julia
using XRadio, FITSFiles
fitsidi2msv4("observation.idifits", "observation.ps.zarr")
ps = open(ProcessingSet, "observation.ps.zarr")
```
"""
function fitsidi2msv4 end
