"""
    XRadio

A Julia implementation of the Measurement Set v4 (MSv4) standard for radio
interferometry and single-dish data.

An MSv4 stores one observation's correlated data as labelled n-dimensional
arrays rather than as tables, together with everything needed to interpret them
— the antennas, the field and source, the pointing, the weather. Several of
them are collected into a *processing set*, which on disk is a single Zarr
store. `XRadio` reads and writes those stores, and checks data against the
standard.

The reference implementation is the Python package
[xradio](https://github.com/casangi/xradio); this package targets the same
on-disk format with a Julia interface.

# Reading

```julia
ps = open(ProcessingSet, "track.ps.zarr")   # lazy: metadata only
ms = ps[1]                                  # one Measurement Set v4
ms[:visibility]                             # the science array
```
"""
module XRadio

using Dates
using Reexport

# A MeasurementSet is a DimTree and a Dataset a DimStack, so selecting, slicing
# and inspecting one is done with DimensionalData's vocabulary — `At`, `Near`,
# `dims`, `lookup`, `metadata`.
@reexport using DimensionalData

using DimensionalData: @dim, Dimension, DimArray, DimStack, DimTree, branches, val
using DimensionalData.Lookups: Lookups
using OrderedCollections: OrderedDict

import JSON3
import Zarr

export ProcessingSet, MeasurementSet, Dataset, Measure
export fitsidi2msv4
export datagroups, datagroup, correlated, flags, weights, uvw
export antennas, scans, fields, intents, frequencies, polarizations, times, spectralwindow
export BaselineId, Frequency, Polarization, UVWLabel, AntennaName, StationName,
    FieldName, ReceptorLabel, ElementId, CartesianPosLabel, CartesianPosLabelLocal,
    EllipsoidDirLabel, EllipsoidDisLabel, SkyDirLabel, SkyDisLabel, LocalSkyDirLabel,
    LineLabel, ToneLabel, PolyTerm, NPolynomial, TimeEphemeris, TimePointing,
    TimeSystemCal, TimeWeather, TimePhaseCal, FrequencySystemCal

include("dimensions.jl")
include("schema/spec.jl")
include("schema/measures.jl")
include("schema/infodicts.jl")
include("schema/correlated.jl")
include("schema/subdatasets.jl")
include("schema/registry.jl")
include("types.jl")
include("check.jl")
include("zarr/read.jl")
include("zarr/write.jl")
include("accessors.jl")
include("fitsidi.jl")

end # module
