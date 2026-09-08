# XRadio.jl

A Julia implementation of the **Measurement Set v4 (MSv4)** standard for radio
interferometry and single-dish data.

An MSv4 stores one observation's correlated data as labelled n-dimensional
arrays rather than as tables, together with everything needed to interpret them:
the antennas, the field and source, the pointing, the weather. Several are
collected into a *processing set*, which on disk is a single Zarr store.
`XRadio` reads and writes those stores and checks data against the standard.

The reference implementation is the Python package
[xradio](https://github.com/casangi/xradio); this package targets the same
on-disk format with a Julia interface. Stores written here are read by `xradio`,
and stores written by `xradio` are read here.

`XRadio` is early and unregistered. The API will change.

## Reading

```julia
using XRadio

ps = open(ProcessingSet, "track.ps.zarr")   # lazy: metadata only
ms = ps[1]                                   # one Measurement Set v4

correlated(ms)          # the visibilities, via the base data group
antennas(ms)            # ["DV02_A015", "DV06_T704", …]
frequencies(ms)         # channel centres, in Hz
scans(ms), fields(ms)
```

`read(ProcessingSet, path)` does the same but brings every array into memory.

## Selecting

A `MeasurementSet` is a `DimensionalData.DimTree`, so everything that library
offers applies:

```julia
vis = ms[:visibility]

vis[Polarization(At("XX"))]          # one polarization product
vis[Frequency(1:64)]                 # the first 64 channels
ms[Ti(At(t))]                        # a time step, across every layer

ms.antenna[:antenna_position]        # a sub-dataset
```

## Writing

```julia
write("out.ps.zarr", ps)
write("out.ps.zarr", ps; chunks = Dict(:time => 64, :frequency => 128))
```

Data is checked against the standard before anything is written, so a store this
produces conforms or the call fails.

## Checking

Conformance checking is public but not exported, since reading and writing is
what most work needs:

```julia
using XRadio: check, validate

issues = check(ps)          # every departure from the standard
isvalid(ps)                 # or just whether there are any
validate(ps)                # or raise them as an error
```

`check(ps; strict = false)` relaxes only the element types, accepting any type of
the same kind as the standard names — see *Element types* below.

`check` returns `XRadio.SchemaIssue`s; `validate` throws an
`XRadio.SchemaError` carrying them.

## Things worth knowing

### Names

MSv4 spells data variables in upper snake case (`VISIBILITY`, `ANTENNA_POSITION`)
and coordinates in lower snake case (`time`, `baseline_id`). Those are the
on-disk identifiers. In memory `XRadio` uses the lower-case spelling throughout —
`ms[:visibility]`, `ms.antenna[:antenna_position]` — and translates on the way in
and out. Sub-dataset nodes lose the `_xds` suffix, which names the xarray type
they are stored as and says nothing about the data: `antenna_xds` is
`ms.antenna`.

### Dimension order is reversed, and that is fine

Zarr stores arrays in C order, so the first dimension MSv4 names is the
*slowest* varying. Julia arrays are column-major, and `Zarr.jl` maps the store
onto them without copying, which puts the fastest-varying dimension first.
`XRadio` keeps that order:

```julia
size(ms[:visibility])   # (polarization, frequency, baseline_id, time)
                        # xradio reports (time, baseline_id, frequency, polarization)
```

Permuting back would either copy the whole array or destroy its contiguity, and
neither is worth it when every dimension is labelled — `vis[Frequency(At(f))]`
selects the same data whichever end of the tuple frequency sits at. The store
itself is written in the correct MSv4 order, so this is a difference in how the
data is presented, not in what is on disk.

### Element types

The standard fixes the concrete types a conforming *store* may hold. It is
already flexible for the science arrays — weights and `UVW` may be any of
`float16`, `float32` or `float64`, visibilities either complex width — but it
pins the metadata quantities to `float64`.

In memory `XRadio` imposes nothing. Containers are parametric, so a `BFloat16`
visibility or a dual number carried through a derivative is held and operated on
without complaint. The two checks part company there: `check(ps)` reports an
element type the standard does not permit on disk, and `check(ps; strict = false)`
accepts any type of the same kind.

### Zarr versions

Both occur in the wild: MSv4 stores written before 2026 are Zarr v2, and current
tooling writes v3. `XRadio` reads both, including the `dimension_names` field
that v3 keeps in the array metadata and `Zarr.jl` does not surface.

Writing defaults to **v2**, because `Zarr.jl` cannot yet write the fixed-width
strings the standard specifies for v3. Pass `zarr_format = 3` for a v3 store;
it reads back correctly everywhere, but its string arrays are variable-length
and a strict check against the standard calls their dtype wrong.

Every store this package writes carries **consolidated metadata**, gathered as a
finalize step, so opening one costs a single read rather than one per node — a
processing set partitioned per scan runs to hundreds of nodes. Consolidated
metadata is not part of the Zarr v3 specification; it is written to a v3 store
anyway, because the reference implementation reads it, but `XRadio` does not
read it back there. A store that has none, as one written by another tool may,
opens node by node instead.

Four `Zarr.jl` gaps sit behind this, all in Zarr v3 and all worth fixing
upstream: `dimension_names` is parsed but discarded on read and unavailable on
write; complex arrays are written with a `data_type` and `fill_value` that
neither the specification nor `Zarr.jl`'s own reader accepts; fixed-width
strings can be read but not written; and consolidated metadata can be written
but reads back as an empty group. Reproductions for all four, depending only on
`Zarr`, are in `zarr-mwes/`.

## Conformance

The schema tables in `src/schema/` are written by hand, so they read as Julia
rather than as transcribed Python. To keep that from drifting from the standard,
`test/reference/*.json` holds the machine-readable schema exported from the
reference implementation, and the test suite compares every dataset, variable,
dimension tuple, dtype set, optionality flag and permitted attribute value
against it.

The other direction is checked by hand rather than in CI, since it needs a
Python environment: a store written here passes `xradio`'s own
`check_datatree` with no issues, and `xradio` reads it with the dimensions in
MSv4 order and the data unchanged.

That check has been run at full scale on real correlator output. A 24.73 GB
VLBA FITS-IDI file (project BT164A, Q band, ten antennas, four spectral
windows, 752,160 `UV_DATA` rows) converts to a processing set of 144
Measurement Sets — four spectral windows by thirty-six inferred scans — and
`check_datatree` reports no issues across its 1004 datasets. Three properties
the schema checker cannot see hold as well: the `time`×`baseline_id` cells sum
to exactly the file's row count in each spectral window, so the conversion
neither loses nor duplicates a row; none of the 3,080,847,360 visibility cells
is padding, which is what partitioning by scan is for; and the phase center
reads back as 3C273's catalogue position, which exercises the
degrees-to-radians conversion against a number that can be looked up.

Regenerate the reference against a new `xradio` release with:

```python
from xradio.schema.export import export_schema_json_file
import xradio.measurement_set.schema as s
export_schema_json_file(s.VisibilityXds, "test/reference/VisibilityXds.json")
```

## Not included

- **MSv2 conversion.** MSv4 is reached by converting once, with `xradio`'s
  `convert_msv2_to_processing_set`, which needs casacore. There is no Julia
  path to that and no plan for one.
- **Images.** `xradio` also defines an image data model; this package covers
  the Measurement Set only.
- **Attributes nested inside attributes.** `check` validates an array's own
  attributes but does not descend into the labelled arrays some of them hold,
  such as the frequency axis's `reference_frequency`.

## Status against the standard

Implements MSv4 schema version 4.0.0: the `visibility`, `radiometer` and
`spectrum` correlated datasets, and the `antenna`, `pointing`, `weather`,
`system_calibration`, `gain_curve`, `phase_calibration`, `phased_array`,
`field_and_source` and `field_and_source_ephemeris` sub-datasets.
