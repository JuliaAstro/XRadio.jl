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

using Compat: @compat
using Dates
using EnumX: @enumx
using Reexport

# A MeasurementSet is a DimTree and a Dataset a DimStack, so selecting, slicing
# and inspecting one is done with DimensionalData's vocabulary — `At`, `Near`,
# `dims`, `lookup`, `metadata`.
@reexport using DimensionalData

using DimensionalData: @dim, Dimension, DimArray, DimStack, DimTree, branches, val
using DimensionalData.Lookups: Lookups, Categorical, Sampled
using OrderedCollections: OrderedDict

import JSON
import Zarr

export ProcessingSet, MeasurementSet, Dataset, Measure
export AbstractMount, Mount, MountAltAz, MountEquatorial, MountNasmythR, MountNasmythL,
    MountXY, MountOrbiting, MountOther
export AbstractGainCurve, GainPolynomial, TabulatedGain, SphericalHarmonicGain,
    AtDeclination, gain
export fitsidi2msv4, fitsidi_gain_curves
export datagroups, datagroup, correlated, flags, weights, uvw
export antennas, scans, fields, scan_intents, spectral_window_intents,
    Intent, ScanIntents, SubscanIntents,
    frequencies, polarizations, times, spectralwindow,
    gaincurves, baseline, baselines, select, transform
export BaselineID, Frequency, Polarization, UVWLabel, AntennaName, StationName,
    FieldName, ReceptorLabel, ElementID, CartesianPosLabel, CartesianPosLabelLocal,
    EllipsoidDirLabel, EllipsoidDisLabel, SkyDirLabel, SkyDisLabel, LocalSkyDirLabel,
    LineLabel, ToneLabel, PolyTerm, NPolynomial, TimeEphemeris, TimePointing,
    TimeSystemCal, TimeWeather, TimePhaseCal, FrequencySystemCal

# Public, but not exported: reachable as `XRadio.name`, and part of what the
# package promises. Checking a store against the standard is not what most work
# with one involves, and the gain-curve interface is what an outside curve type
# implements rather than what a caller types, so neither earns a name in the
# caller's namespace.
@compat public check, validate, SchemaIssue, SchemaError
@compat public gain_variable, gain_quantity, msv4_gain_curve

include("dimensions.jl")
include("schema/spec.jl")
include("schema/measures.jl")
include("schema/infodicts.jl")
include("schema/correlated.jl")
include("schema/subdatasets.jl")
include("schema/registry.jl")
include("intents.jl")
include("types.jl")
include("build.jl")
include("mount.jl")
include("gaincurve.jl")
include("check.jl")
include("zarr/read.jl")
include("zarr/write.jl")
include("accessors.jl")
include("fitsidi.jl")

end # module
