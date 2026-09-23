# XRadio.jl

A Julia implementation of the **Measurement Set v4 (MSv4)** standard for radio
interferometry and single-dish data.

An MSv4 stores one observation's correlated data as labeled n-dimensional
arrays rather than as tables, together with everything needed to interpret
them: the antennas, the field and source, the pointing, the weather. Several
are collected into a *processing set*, which on disk is a single Zarr store.
XRadio reads and writes those stores and checks data against the standard.

The reference implementation is the Python package
[xradio](https://github.com/casangi/xradio); this package targets the same
on-disk format with a Julia interface.

```julia
using XRadio

ps = open(ProcessingSet, "track.ps.zarr")   # lazy: metadata only
ms = ps[1]                                   # one Measurement Set v4
correlated(ms)                               # the visibilities
```

- [Selecting data](@ref) covers choosing Measurement Sets from a processing set
  and cutting inside one.
- [Coming from xradio](@ref) maps the Python package's calls to this one's.
- [Reference](@ref) lists the public API.
