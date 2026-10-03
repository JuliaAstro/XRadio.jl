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

# A MeasurementSet and a Dataset are both DimTrees, so selecting, slicing and
# inspecting one is done with DimensionalData's vocabulary — `At`, `Near`,
# `dims`, `lookup`, `metadata`.
@reexport using DimensionalData

using DimensionalData: @dim, Dimension, DimArray, DimStack, DimTree, branches, val
using DimensionalData.Lookups: Lookups, Categorical, Sampled
using IntervalSets: AbstractInterval, leftendpoint
using OrderedCollections: OrderedDict

import JSON
import Zarr

export ProcessingSet, MeasurementSet, Dataset, Measure
export AbstractMount, Mount, MountAltAz, MountEquatorial, MountNasmythR, MountNasmythL,
    MountBWGR, MountBWGL, MountRotator, MountXY, MountOrbiting, MountHALCA, MountOther
export AbstractGainCurve, GainPolynomial, TabulatedGain, SphericalHarmonicGain,
    AtDeclination, gain
export fitsidi2msv4, fitsidi_gain_curves
export datagroups, datagroup, add_datagroup!, correlated, flags, weights, uvw
export antennas, mounts, antenna_positions, receptor_angles, polarization_types,
    scans, fields, sources, scan_intents, spectral_window_intents,
    Intent, ScanIntents, SubscanIntents,
    frequencies, polarizations, times, spectralwindow,
    gaincurves, system_temperatures, baseline, baselines, baseline_antennas, scan, field,
    field_and_source, partition_info,
    query, transform, ByScan, ChannelBins, ByChannel, average
export BaselineID, Frequency, Polarization, UVWLabel, AntennaName, StationName,
    FieldName, ReceptorLabel, ElementID, CartesianPosLabel, CartesianPosLabelLocal,
    EllipsoidDirLabel, EllipsoidDisLabel, SkyDirLabel, SkyDisLabel, LocalSkyDirLabel,
    LineLabel, ToneLabel, PolyTerm, NPolynomial, TimeEphemeris, TimePointing,
    TimeSystemCal, TimeWeather, TimePhaseCal, FrequencySystemCal

# Public, but not exported: reachable as `XRadio.name`, and part of what the
# package promises. Checking a store against the standard is not what most work
# with one involves, the gain-curve interface is what an outside curve type
# implements rather than what a caller types, and the lookups are what a solver
# holds to keep a hot loop off the Measurement Set, so none earns a name in the
# caller's namespace.
@compat public check, validate, SchemaIssue, SchemaError
@compat public DatasetSpec, ArraySpec, AttrSpec, DictSpec, schema_types, extend
# `extend` builds on a schema, so the schemas themselves are part of the
# extension surface.
@compat public SCHEMAS, VISIBILITY_SCHEMA, SPECTRUM_SCHEMA, ANTENNA_SCHEMA,
    POINTING_SCHEMA, WEATHER_SCHEMA, SYSTEM_CALIBRATION_SCHEMA, GAIN_CURVE_SCHEMA,
    PHASE_CALIBRATION_SCHEMA, PHASED_ARRAY_SCHEMA, FIELD_SOURCE_SCHEMA,
    FIELD_SOURCE_EPHEMERIS_SCHEMA
@compat public gain_variable, gain_quantity, msv4_gain_curve
@compat public set_mounts!, mount_from_msv4, msv4_mount, axis_offset
@compat public baseline_lookup, BaselineLookup
@compat public subdataset, measurement_set
@compat public scan_lookup, field_lookup, LabelLookup
@compat public averaging_axis

include("dimensions.jl")
include("schema/spec.jl")
include("schema/measures.jl")
include("schema/infodicts.jl")
include("schema/correlated.jl")
include("schema/subdatasets.jl")
include("schema/registry.jl")
include("intents.jl")
include("baselines.jl")
include("labels.jl")
include("types.jl")
include("indexing.jl")
include("build.jl")
include("mount.jl")
include("gaincurve.jl")
include("check.jl")
include("zarr/read.jl")
include("zarr/chunks.jl")
include("zarr/write.jl")
include("accessors.jl")
include("grouping.jl")
include("averaging.jl")
include("fitsidi.jl")
include("testing.jl")

end # module
