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
scan_intents(ms)        # [Intent("OBSERVE_TARGET#ON_SOURCE")]
```

Scan intents come back as `Intent`s rather than strings. Both halves are closed
ASDM enumerations, so they are compared against names the compiler checks
instead of strings a typo makes silently unequal:

```julia
Intent(ScanIntents.OBSERVE_TARGET, SubscanIntents.ON_SOURCE) in scan_intents(ms)
```

An intent outside those enumerations throws when read, naming the value.
`spectral_window_intents(ms)` stays strings, since that vocabulary is not
enumerated.

`read(ProcessingSet, path)` does the same but brings every array into memory.

## Selecting

`select` narrows a processing set to the Measurement Sets stating the
attributes asked for, by `antenna`, `field`, `scan` or `spectralwindow`:

```julia
select(ps; field = "M87")                       # one source
select(ps; field = "M87", spectralwindow = "spw_2")
select(ps; scan = ["3", "4", "5"])              # any of three scans
```

A keyword takes one name or a collection of them, several keywords must all
match, and an unrecognized keyword is an error rather than a selection that
quietly keeps everything. The result is a `ProcessingSet`, so selections chain;
`filter(f, ps)` remains for predicates the keywords do not reach.

`transform(f, ps)` turns each Measurement Set into another one, keeping the
names, so the result is a processing set again. Reading a value out of each
instead gives them in order, as for any collection:

```julia
transform(ms -> delete!(deepcopy(ms), :weight), ps)   # ProcessingSet
map(gaincurves, values(ps))                          # Vector, in ps's order
```

Within one Measurement Set, a `MeasurementSet` is a `DimensionalData.DimTree`,
so everything that library offers applies:

```julia
vis = ms[:visibility]

vis[Polarization(At("XX"))]          # one polarization product
vis[Frequency(1:64)]                 # the first 64 channels
ms[Ti(At(t))]                        # a time step, across every layer

ms.antenna[:antenna_position]        # a sub-dataset
```

Baselines are named by their two antennas, rather than by the integer
`baseline_id` the store indexes them with:

```julia
baselines(ms)                                   # [("DV02_A015", "DV06_T704"), …]
vis[baseline(ms, ("DV02_A015", "DV06_T704"))]   # one baseline, either order
vis[baseline(ms, "DV02_A015-DV06_T704")]        # the same pair as one label
vis[baselines(ms, "DV02_A015"; auto = false)]   # one antenna's cross-correlations
```

Both return a selector rather than data, so they compose with everything above
and leave the choice of layer to the caller.

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

### Checking entries of your own

MSv4 is extensible, and a conforming checker accepts what it does not recognize:
extra attributes unconditionally, and extra coordinates and data variables as
long as they span dimensions the dataset already has. `check` does the same, so
data of your own passes without XRadio knowing anything about it — and without
being examined.

To have it examined, say what it must look like. `XRadio.extend` copies a
standard specification with entries added; `schemas` puts the result in place of
the standard one for that `type`, in `check`, `validate`, `isvalid` and `write`
alike:

```julia
const MY_VISIBILITY = XRadio.extend(
    XRadio.VISIBILITY_SCHEMA;
    attrs = [XRadio.AttrSpec(:earth_orientation, :dict; optional = true)],
    coords = [
        # An entry the standard already defines gains these attributes and
        # keeps everything else.
        XRadio.ArraySpec(
            "frequency", (:frequency,), Float64;
            attrs = [XRadio.AttrSpec(:sideband, :int; optional = true)]
        ),
        XRadio.ArraySpec("sub_scan_name", (:time,), String; optional = true),
    ],
)

check(ps; schemas = [MY_VISIBILITY])
write(path, ps; schemas = [MY_VISIBILITY])
```

Passing it to `write` is what gives an added coordinate its lower-case on-disk
name and lists it among the `coordinates` of the variables it labels; a layer no
specification names is written upper case, as MSv4 spells its data variables.

Nothing is registered globally, so a store checks the same way whatever else has
been loaded. A whole dataset of a kind the standard does not define is the same
mechanism: pass a `DatasetSpec` whose `type` is your own.

## Converting FITS-IDI

FITS-IDI (AIPS Memo 114) is the interchange format VLBI correlators write.
Loading `FITSFiles` alongside `XRadio` enables the converter:

```julia
using XRadio, FITSFiles

fitsidi2msv4("observation.idifits", "observation.ps.zarr")
ps = open(ProcessingSet, "observation.ps.zarr")
```

Each partition is written as it is built, so peak memory is one partition
whatever the input weighs. A 24.73 GB file converts in about 4 GB.

### What maps to what

| FITS-IDI | MSv4 |
| --- | --- |
| `UV_DATA` | the `visibility` dataset: row-major `(time, baseline)` becomes a dense `(time, baseline_id, frequency, polarization)` array |
| `ARRAY_GEOMETRY` + `ANTENNA` | `antenna_xds` — MSv4 merges antenna and feed, as FITS-IDI does |
| `SOURCE` | `field_and_source_base_xds` |
| `FREQUENCY` | the `frequency` coordinate's attributes |
| `FLAG` | the `FLAG` array, expanded from ranges |
| `SYSTEM_TEMPERATURE` | `system_calibration_xds` |
| `PHASE-CAL` | `phase_calibration_xds` |
| `WEATHER` | `weather_xds` |
| `GAIN_CURVE` | `gain_curve_xds` |

One Measurement Set covers a single spectral window, so the `BAND` axis is
always split. A partition carries only the antennas its own baselines name and
the fields its own rows observed, and carries a calibration or environment
sub-dataset only where measurements cover its time span — a scan earlier than
the first weather reading gets no `weather_xds`.

Units become MSv4's throughout: seconds to meters for `UVW`, degrees to
radians, centigrade to kelvin, millibars to pascals.

### Scans, and how they are inferred

AIPS Memo 114r defines no scan column, so by default a scan boundary is
inferred where the source changes **or** where more than `gap = 60.0` seconds
separate consecutive timestamps. A file that carries a `SCAN` column anyway is
believed instead, since it states what an inference can only guess.

Scans matter because `partition = (:scan,)` is the default. Splitting by scan
makes each partition's baseline axis describe antennas that actually observed
together, so no padding is needed; `(:source,)` splits only where the source
changes, and `()` puts a whole observation in one partition per spectral
window, which only a small file can afford.

### What the converter invents

MSv4 asks for several things FITS-IDI does not record. Each is listed here with
the value used when the file is silent, and whether it says so:

| MSv4 field | Value when the file does not say | Warns |
| --- | --- | --- |
| `observation_info.release_date` | `1858-11-17T00:00:00.000`, the MJD zero point | yes |
| `VISIBILITY.units` | `uncalib`, when the `FLUX` column has no `TUNIT` | yes |
| `gain_curve_xds` `measured_date` | the MJD zero point | yes |
| `antenna_xds` `mount` | `ORBITING` or `OTHER`, for `MNTSTA` 2 and 6, which MSv4 reserves no name for | no |
| `polarization` | the majority receptor basis, for an array that does not share one | yes |
| `scan_intents` | `["UNSPECIFIED"]` | no |
| `spectral_window_intents` | `["UNSPECIFIED"]` | no |
| `gain_curve_xds` `GAIN_CURVE_INTERVAL` | the observation's whole duration, since FITS-IDI states one curve per file | no |
| `processor_info.sub_type` | `UNKNOWN`, when there is no `CORRELAT` keyword | no |
| `scan_name` | the inferred segmentation above | no |
| `spectral_window_name` | `spw_0`, `spw_1`, … from the `BAND` index | no |

`release_date`, `scan_intents` and `spectral_window_intents` are keywords of
the same name, so a caller who knows the real values can state them.
`scan_intents` takes `Intent`s, or strings that parse as one, so an intent
outside the ASDM enumerations is rejected rather than written. Only
`release_date` warns, because a sentinel date could be mistaken for a real
embargo.

`INTTIM` is a fourth case, but the converter refuses to guess it: MSv4's time
axis requires an integration time and FITS-IDI makes the column optional, so a
file without one needs the `integration_time` keyword. Inferring it from the
spacing of the timestamps would present a guess as a measurement.

### What is dropped

Five FITS-IDI tables have no MSv4 home. Where the file holds one, it is named
in a warning and its content does not reach the store: `INTERFEROMETER_MODEL`,
`MODEL_COMPS`, `BASELINE`, `CALIBRATION` and `BANDPASS`. A table the standard
does not define at all — a correlator's own extension, such as DiFX's `CALC` —
is ignored without comment.

Three columns of tables that are otherwise carried have no MSv4 field, and are
dropped silently: `PHASE-CAL`'s `PC_RATE` and `STATE`, and `WEATHER`'s
`PRECIPITATION`.

`ARRAY_GEOMETRY`'s `STAXOF` has no MSv4 field either, but is kept: the antenna
dataset carries `ANTENNA_AXIS_OFFSET`, `MOUNT_PARALLACTIC` and
`MOUNT_ELEVATION` as extension variables over `antenna_name`, so `mounts(ms)`
gives each mount back with its axis offset. A reader that does not know them
ignores them, as MSv4 permits.

A gain curve that is not a polynomial in elevation or in zenith angle — a
tabulated curve, or a spherical-harmonic expansion — cannot be written at all,
because MSv4 models only the polynomial. Such a file gets no `gain_curve_xds`
and a warning saying so. `fitsidi_gain_curves` reads every curve
whatever its shape, so one the store cannot hold is still reachable from the
file.

### What changes value

Two conversions rescale the data rather than relabeling it. `WEIGHTYP =
"CORRELAT"` means the `WEIGHT` column holds validity fractions rather than the
`1/σ²` MSv4 wants, so each is converted by `f · 2Δν · τ · η²`, with `η` the
correlator efficiency from the `weight_efficiency` keyword; `weight_mode`
overrides the choice of whether to convert at all. And a `VIS_SCAL` keyword,
which lets a near-real-time correlator write accumulated sums, divides into
every visibility amplitude.

Flagging leaves the weights alone: `FLAG` says a datum must not be used, and
`WEIGHT` says how good it would have been.

### Where the converter stops

It fails rather than guessing, in every case with a message naming what it
found:

- **A B1950 phase center.** MSv4 restricts `UVW.frame` to `icrs`, `fk5` and
  `APP`, so an fk4 direction has no conforming label. The error says to precess
  the observation to J2000 first.
- **`UU/VV/WW---NCP`.** Its `w` axis points at the north pole rather than at
  the source, so reading it as `---SIN` would be wrong without complaining.
- **More than one `ARRAY_GEOMETRY` table.** An antenna number is unique only
  within its array, and the `256·a1 + a2` `BASELINE` packing carries no array,
  so two numberings would go through one map.
- **More than one row in `FREQUENCY`.** Several frequency setups would be
  silently mixed onto one axis.
- **No `INTTIM` column and no `integration_time` keyword**, as above.
- **A non-positive `VIS_SCAL`**, which divides into every amplitude.

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
