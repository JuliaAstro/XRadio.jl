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

One Measurement Set covers a single spectral window, so a file's `BAND` axis is
always split. It is split by scan as well by default, because a scan's baseline
axis then describes antennas that actually observed together and needs no
padding. AIPS Memo 114r defines no scan column, so scans are normally
*inferred*: a scan ends where the source changes or where more than `gap`
seconds pass between consecutive timestamps. A correlator that writes a `SCAN`
column anyway states what an inference can only guess, and it is used instead.

Three MSv4 fields have no FITS-IDI source and are given defaults, each
overridable by the keyword of the same name: `scan_intents`,
`spectral_window_intents` and `release_date`. Only `release_date` warns, since
a sentinel date could be mistaken for a real embargo policy.

Each Measurement Set carries the antennas its own baselines name and the fields
its own rows observed, so a scan that ran on a subset of the array describes
that subset rather than the whole. `STAXOF` has no field in MSv4's `antenna`
dataset and is dropped; every other column of `ARRAY_GEOMETRY`, `ANTENNA` and
`SOURCE` that MSv4 defines a home for is carried.

A `FLAG` table is expanded onto the dense `FLAG` array. Its rows are *ranges* —
a source, an antenna pair with `0` as a wildcard, a frequency setup, a time
range, a band selector, a channel range and a polarization selector — and
nothing in `UV_DATA` repeats what they say, so a store written without them
would present flagged data as good. Flagging leaves the weights alone: `FLAG`
says the datum must not be used and `WEIGHT` says how good it would have been.

# Keywords
- `partition = (:scan,)`: how to split beyond the mandatory spectral window.
  `(:source,)` splits where the source changes but not across time gaps, and
  `()` puts a whole observation in one partition per spectral window — which
  only a small file can afford.
- `gap = 60.0`: seconds between timestamps above which a scan is taken to have
  ended.
- `scan_intents = ["UNSPECIFIED"]`, `spectral_window_intents = ["UNSPECIFIED"]`.
- `release_date`: ISO 8601 date on which the data may become public. Defaults
  to the MJD zero point, with a warning.
- `weight_mode = :auto`: how the on-disk `WEIGHT` column becomes an MSv4
  weight. `:auto` follows the `WEIGHTYP` keyword, which usually says the column
  holds validity fractions rather than the `1/σ²` MSv4 wants and so needs the
  radiometer conversion `f · 2Δν · τ · η²`. `:radiometer` converts regardless,
  `:validity` writes the column through unchanged.
- `weight_efficiency = 1.0`: the correlator efficiency `η` of that conversion.
  It sets the absolute weight scale, not the relative weighting.
- `integration_time`: seconds, required only for a file with no `INTTIM`
  column. FITS-IDI makes that column optional and MSv4's time axis does not, so
  the value has to come from somewhere; inferring it from the spacing of the
  timestamps would present a guess as a measurement.
- `element_type = Float32`: element type of the weights, and the real type
  behind the complex visibilities. `UVW` is always `Float64`: a VLBI baseline
  is around 1e7 metres and fixes a phase that turns over in a wavelength, which
  `Float32` cannot resolve.
- `name`: base of the Measurement Set names, which are `name_00`, `name_01` and
  so on. Defaults to the input file's basename.
- `zarr_format = 2`, `compressor`, `chunks`, `overwrite`, `validate`: as for
  [`write`](@ref). `chunks` defaults to 128 times by 64 channels, since a
  converted store is large enough that one chunk per array is unusable. Each
  partition is checked as it is built, so a departure from the standard stops
  the conversion where it arose rather than at the end.

Requires `FITSFiles`:

```julia
using XRadio, FITSFiles
fitsidi2msv4("observation.idifits", "observation.ps.zarr")
ps = open(ProcessingSet, "observation.ps.zarr")
```
"""
function fitsidi2msv4 end
