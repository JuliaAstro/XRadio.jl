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

The four optional calibration and environment tables become sub-datasets of
every partition they have something to say about: `SYSTEM_TEMPERATURE` becomes
`system_calibration_xds`, `PHASE-CAL` becomes `phase_calibration_xds`,
`WEATHER` becomes `weather_xds` and `GAIN_CURVE` becomes `gain_curve_xds`.
FITS-IDI samples the first three on the instrument's own clock, and each keeps
that clock rather than being interpolated onto the visibilities' time axis. A
partition takes the rows whose interval overlaps its own time span, so a scan
is described by the measurements covering it, and carries no such sub-dataset
where nothing does. Weather measured for the site as a whole — `ANTENNA_NO` of
zero — is recorded against every station.

Units become MSv4's: kelvin, pascals and radians rather than FITS-IDI's
centigrade, millibars and degrees. A measurement the file writes as `NaN` stays
`NaN`, which is what the standard means by one it does not have.

`VISIBILITY` is the exception, because MSv4 fixes the units of every other array
to one permitted value but leaves this one a free string with `Jy` as its
default. Correlator output is usually not in janskys, so the store says what the
`FLUX` column's `TUNIT` says — `UNCALIB` for uncalibrated correlation
coefficients, which is what a VLBI correlator writes. `TUNIT` is optional in
FITS, and a file that omits it gets `uncalib` and a warning, since the
alternative is to claim a calibration the file never asserted.

Four things these tables carry have no MSv4 field and are dropped:
`PHASE-CAL`'s `PC_RATE` and `STATE`, `WEATHER`'s `PRECIPITATION`, and any
`GAIN_CURVE` curve that is not a polynomial in elevation or in zenith angle —
MSv4 models only the polynomial, so a file whose curves are tabulated or are
spherical-harmonic expansions gets no `gain_curve_xds` and is told so.
[`fitsidi_gain_curves`](@ref) reads all of them whatever their shape, so a
curve the store cannot hold is still reachable from the file. MSv4 also asks
when a gain curve was measured and FITS-IDI does not say, so `measured_date` is
the MJD zero point, with a warning. Its `GAIN_CURVE_INTERVAL` is the
observation's own duration, since FITS-IDI states one curve for a whole file,
and its coefficients are against the angle in **degrees**, which is the
convention casacore and AIPS write and MSv4 inherits without stating.

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

"""
    fitsidi_gain_curves(path) -> DimStack

Every gain curve the FITS-IDI file at `path` records, over the antennas,
receptors and spectral windows it describes them for.

The result has two layers — `curve`, holding an
[`XRadio.AbstractGainCurve`](@ref), and `sensitivity`, its `SENS` in kelvin per
jansky, which turns a dimensionless curve into a flux scale — selected the same
way as anything else in XRadio:

```julia
curves = fitsidi_gain_curves("observation.idifits")
curves[:curve][AntennaName(At("BR")), ReceptorLabel(At("pol_0")), Frequency(Near(43e9))]
```

The `frequency` axis holds one entry per spectral window rather than per
channel, since that is the resolution a gain curve is stated at, and each entry
is that window's reference frequency — the same value a converted store carries
on its own frequency axis. One Measurement Set covers one spectral window, so
[`fitsidi2msv4`](@ref) splits across what this keeps together. A combination
the file states no curve for holds `nothing` rather than an empty curve.

This reads every shape AIPS Memo 114r defines — a polynomial in an angle, a
table of measured values, and the spherical-harmonic expansion in hour angle
and co-declination the Green Bank 140-foot used — and so it reads curves that
[`fitsidi2msv4`](@ref) cannot write. MSv4's `gain_curve_xds` holds a polynomial
in elevation or in zenith angle and nothing else, so a converted store keeps
only those; this is how the rest stay reachable. [`gaincurves`](@ref) reads what
a converted store does keep, in the same shape minus the frequency axis.

Angles are radians, where FITS-IDI writes degrees. A polynomial's coefficients
are scaled to match, so `gain(curve, elevation_in_radians)` is the gain
directly.

Requires `FITSFiles`:

```julia
using XRadio, FITSFiles
curves = fitsidi_gain_curves("observation.idifits")
g = gain(first(curves[:curve]), deg2rad(45))   # gain at 45 degrees elevation
```
"""
function fitsidi_gain_curves end
