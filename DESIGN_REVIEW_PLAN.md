# FITS-IDI → MSv4 converter plan

<!-- Preserve chunk IDs and status values. -->

## Metadata
- **Kind**: `design`
- **Package**: XRadio
- **Written**: 2026-09-05
- **Current version**: 0.1.0
- **Issue**: n/a
- **Goal**: read a FITS-IDI file and write it as a conforming MSv4 processing set

## Stated values

> XRadio must actually meet the MSv4 standard — conformance is checkable, and
> the hand-written schema tables are cross-checked against the reference
> implementation's own export so they cannot drift silently.

Author's working preferences, gathered while building the package:

- **Fail fast over silently continuing.** Unexpected input should surface, not
  be papered over. But do not litter the code with checks that earn nothing.
- **Julian, not transcribed Python.** Follow Julia conventions — `open`/`read`/
  `write` on types, dispatch rather than accessor objects, lower-case keys in
  memory — even where xradio does otherwise. The standard governs the bytes on
  disk, not the shape of the Julia API.
- **Generic in element type.** Containers are parametric; any float or complex
  representation can be held and operated on. Concrete dtypes are a property of
  a conforming *store*, enforced only when checking one.
- **Independent of Gustavo.** Port from it freely — same author, MIT — but do
  not depend on it.
- **Unreleased alpha.** The API breaks freely; no deprecation shims are owed.
- **Verify rather than assume.** Where the specification and real data can
  disagree, measure the real data. Several conversions here are silent when
  wrong.

## Release strategy
- **Pre-breaking-release**: `no` (unregistered; nothing to protect)
- **Inter-cluster releases**: `no` (one terminal release when the converter works)

## Decisions
<!-- Answers to `decide` chunks land here, tagged with the chunk ID. -->

#### CHUNK-002 — partitioning

Default to **spectral window + scan**. On the Q-band file that is 144
partitions of 270 MB at 100% fill, against 4 partitions of 6.81 GB at 90.4%
fill for spectral window alone; peak working memory drops 25-fold and the
9.6% of `NaN` padding disappears entirely. Each partition's baseline axis then
describes antennas that actually observed together, which the spw-only scheme
cannot express — baselines per scan run from 10 to 55.

- Keyword: `partition = (:scan,)`, mirroring xradio's `partition_scheme`.
  `:source` remains available but is not a default; on this file it is wholly
  subsumed by `:scan` (both give 36 scans).
- Node names: `<file basename>_<zero-padded index>`, matching xradio, so a
  store holding both converted and xradio-written data stays coherent.
- 144 nodes means several thousand metadata files, so CHUNK-015 (consolidated
  metadata) must land before CHUNK-012 validates a real conversion — the
  dependency is recorded there.

#### CHUNK-003 — invented metadata

Only three MSv4 fields lack a FITS-IDI source, and each takes the value the
reference implementation already writes:

| field | value | overridable |
| --- | --- | --- |
| `spectral_window_intents` | `["UNSPECIFIED"]` | keyword |
| `scan_intents` | `["UNSPECIFIED"]` | keyword |
| `observation_info.release_date` | `"1858-11-17T00:00:00.000"` | keyword |

Warn for `release_date` only. A sentinel date could be mistaken for a real
embargo policy, whereas `["UNSPECIFIED"]` is what xradio itself writes from
MSv2 and so is normal content rather than a shortfall. The source file name
goes into `data_groups.base.description` as provenance.

`observation_info.observer` comes from the **primary header's** `OBSERVER`,
falling back to `UV_DATA`'s. The two disagree in the test file — primary says
`"BT164A"`, `UV_DATA` says `"PLUTO"` — and the primary header is the
file-level authority, which is where AIPS Memo 114 places the keyword. This
duplicates `project_UID` (also `"BT164A"`, from `OBSCODE`) on VLBA data, which
is an honest reflection of the file rather than a defect to paper over.

#### CHUNK-004 — time and frames

**Time**: `format = "unix"`, `scale = "utc"`, `units = "s"`.

The plan recommended `mjd` on the grounds that it avoids an epoch shift and
keeps precision. Both grounds were wrong, and testing rather than assuming is
what showed it. MSv4 pins the time coordinate's `units` to `"s"` and defines
decoding as `Time(data * u.Unit(units), format=format, scale=scale)`. Under
astropy 8.0.1, `Time(60617.5 * u.s, format="mjd")` yields **1858-11-17**
rather than 2024-11-03 — writing MJD in days under `units = "s"` is wrong by
166 years and raises nothing. `format = "mjd"` therefore obliges seconds since
MJD 0, so:

- the epoch shift is unavoidable either way;
- the precision claim inverts. MJD-in-seconds is 5.27e9 for this observation
  against unix's 1.77e9, so its `Float64` ulp is 9.5e-7 s against 2.4e-7 s —
  four times coarser.

`unix` is additionally the MSv4 default, XRadio's own schema default, and what
existing stores carry. `Dates.datetime2unix` supplies POSIX seconds from
`RDATE` with no new dependency.

**Frames**: `EQUINOX` maps literally — `"J2000"` → `"fk5"`, which is exactly
xradio's `Ref_map` in `_utils/schema.py`. Anything else is an error rather
than a silent default to `icrs`; the difference is ~25 mas and invisible in
the output. `UVW.frame` takes the phase center's frame, following xradio's
`fix_uvw_frame`, which distrusts any separately declared UVW frame.

Source directions come from `RAEPO`/`DECEPO` (degrees → radians).
`RAAPP`/`DECAPP` are **unusable**: DiFX wrote them as verbatim copies of the
mean places, so they are not apparent places despite the column names.

**No time-scale dependency.** `TIMSYS = "UTC"` throughout, and the file
supplies `IATUTC = 37.0` and `UT1UTC = 0.077593` explicitly if ever needed.
Converting `RDATE` + `DATE` + `TIME` to UTC seconds is arithmetic on values
FITS-IDI already provides.

#### CHUNK-005 — polarization is labelled by receptor, not by basis

Follow **Gustavo's** representation rather than MSv4's. The `UV_DATA` `STOKES`
axis is decoded to generic correlation products — `PP`, `QQ`, `PQ`, `QP`, by
receptor index — and the polarization basis stays **per antenna**, from the
`ANTENNA` table's `POLTYA`/`POLTYB`. The circular block (-1..-4) and the linear
block (-5..-8) therefore decode identically; they differ only in what the
antennas' receptors carry.

MSv4's array-wide `polarization` coordinate of `"RR"`, `"RL"`, … presumes every
antenna shares a basis, and a VLBI array need not: mixing circular- and
linear-feed stations is ordinary, and the baseline between them measures `RX`,
`LY`, `RY`, `LX` — products no single array-wide label can name. Labelling the
axis by basis is only correct for a homogeneous array and is silently wrong
otherwise, which is exactly the failure mode the stated values call out.

`baseline_pol_labels(header, nosta1, nosta2)` resolves the products to actual
polarizations for one baseline, which is the level at which the question has an
answer. On the homogeneous VLBA file it returns the familiar
`["RR", "LL", "RL", "LR"]`.

This is a representation choice inside the reader, not a licence to write a
non-conforming store. What CHUNK-007 puts in the MSv4 `polarization` coordinate
is still MSv4's problem, and a mixed-basis array will need an explicit decision
there rather than a silent mislabelling.

#### CHUNK-005 — no intermediate model

There is **no FITS-IDI data model**. `fitsidi2msv4` builds MSv4 datasets
straight from the file's tables, and the extension holds only what that needs:
access primitives for cards and columns, the arithmetic the format states
rather than stores (`pol_axis`, `pol_labels`, `baseline_antennas`,
`channel_frequencies`, `scan_ranges`), and one bulk reader for `UV_DATA`'s
per-row columns.

An earlier pass built eight structs mirroring the tables. It was a
representation for a FITS-IDI *library*, not for a one-way converter: nothing
but the conversion ever reads it, and roughly sixteen of its fields had no MSv4
target at all — proper motion and parallax, which `FIELD_SOURCE_SCHEMA` does
not define; the earth-orientation keywords CHUNK-004 had already ruled out;
`RAAPP`/`DECAPP`, which the same chunk found unusable.

The standing preference is **simple, readable code**. Structs are welcome where
they make the conversion clearer and should be introduced then, driven by the
call sites that want them, rather than anticipated. `read_uv_rows` returns a
`NamedTuple` of columns for the same reason.

Consequence for sequencing: with nothing between the file and the MSv4 arrays,
the conversion cannot be tested without a FITS file. **CHUNK-006 is now a hard
prerequisite for CHUNK-007**, not an optional convenience.

#### CHUNK-005 — where the FITS-IDI code lives

**Entirely in the extension.** XRadio implements MSv4; FITS-IDI is not MSv4,
and a package that carries FITS-IDI vocabulary in its own namespace has widened
its subject. `src/` declares one function, `fitsidi2msv4`, and everything the
conversion needs on the way stays behind `using FITSFiles`.

The cost is that unit tests reach the reader through
`Base.get_extension(XRadio, :XRadioFITSFilesExt)` rather than importing it
directly. That is a documented mechanism and the extension is loaded in the
test environment regardless, so it buys the scope boundary cheaply.

`fitsidi2msv4(path, dest)` **writes to disk and returns `dest`** rather than
returning a `ProcessingSet`. Each partition is written as it is built, so peak
memory is one partition — 270 MB on the validation file — no matter how large
the input. Returning an in-memory set would put the whole 24.73 GB observation
through memory to reach a store the caller then has to write anyway.

`fitsidi2msv4` is the author's name for it. It is not the `read`/`write`-on-a-type
form the plan's stated values otherwise ask for; CHUNK-011 owns the public API
and can settle that.

#### CHUNK-005 — which FITSFiles

Depend on the **development branch**, `ptiede/FITSFiles.jl` at
`ptiede-uvfixes` (four commits ahead of `JuliaAstro/FITSFiles.jl` `main`,
including "Added new lazy reading features using disk arrays and fixed perf").
Author's instruction, and there is no alternative that is not worse.

Registered `0.3.2` and the dev checkout report the same version but are not the
same code. Registered `0.3.2` has no `LazyStructuredData` at all: `hdu.data`
materializes the whole HDU on **every** property access, so merely touching
`f[9].data` on the Q-band file reads all 24.73 GB. Everything CHUNK-001 and
CHUNK-005 measured about the FITSFiles API was measured against the dev
checkout.

The `[compat]` bound is `FITSFiles = "0.3.2"`, which is the closest
expressible thing and is *not* sufficient — the registered release of that
version cannot run the extension. It becomes true once the four commits are
merged and released; until then a manifest pinning the branch is required, and
XRadio cannot be registered.

#### CHUNK-014 — streaming architecture

**Deferred, and no longer a blocker.** With spectral window + scan as the
CHUNK-002 default, partitions are 270 MB, so CHUNK-007 builds one partition at
a time in memory and writes it through the existing `write` before starting the
next. Peak memory is one partition, and no new machinery is needed to get
there. The chunk stays open as an optimization: a correlator emitting one long
scan brings the problem straight back, because a single partition is then too
large to hold.

The eventual shape, when something demands it, is a lazy `ProcessingSet` with
a streaming `write` — it makes the existing lazy read path useful for the
first time and keeps one `write` implementation rather than two. A
time-chunked Zarr default (128 times x all baselines x 64 channels x all
polarizations) should be set with CHUNK-007 regardless, since the current
default of one chunk per array is wrong at any of these sizes.

## Why this and not MSv2

MSv4 is reached from MSv2 by `xradio`'s `convert_msv2_to_processing_set`, which
needs casacore and has no Julia path. FITS-IDI is the other route into VLBI data,
it is a pure FITS format that Julia can read today, and MSv4 was explicitly
designed with VLBI in mind — its `gain_curve`, `phase_calibration` and
`system_calibration` sub-datasets exist for exactly the tables FITS-IDI carries.
So this is the conversion that Julia can do end to end and Python currently
cannot.

## What already exists

`Gustavo.jl`'s FITS-IDI reader (`ext/GustavoFITSFilesExt/`, ~4200 lines) is the
same author's MIT-licensed code and is the reference to port from, **not** to
depend on — XRadio stays independent of Gustavo (decision from the session that
created this package).

| File | Lines | Relevance |
| --- | --- | --- |
| `fitsidi_read.jl` | 1834 | The reader. Table parsing, dense regridding, lazy `DiskArrays` chunks, scan segmentation, flag application, weight semantics |
| `fitsidi_write.jl` | 478 | A FITS-IDI writer — the model for building test fixtures |
| `fitsidi_apriori.jl` | 114 | `SYSTEM_TEMPERATURE`, `GAIN_CURVE` |
| `fitsidi_phasecal.jl` | 74 | `PHASE-CAL` |

Directly reusable in shape (all need rewriting against MSv4 rather than `UVSet`):

- `_build_idi_antenna_table` (ARRAY_GEOMETRY + ANTENNA)
- `_build_idi_freq_setups` (FREQUENCY)
- `_build_idi_sources` (SOURCE)
- `_build_idi_flags` / `_filter_flags_for_leaf` (FLAG)
- `_idi_stokes_perm` (AIPS STOKES codes → labels)
- `_idi_segment_scans` (scan boundaries from source changes and time gaps)
- `_fill_vis_dense!` / `_fill_weights_dense!` (row → dense grid with NaN padding)

`FITSFiles.jl` 0.3.2 is registered and is what Gustavo uses.

## Table mapping

This is the substance of the work. FITS-IDI table (AIPS Memo 114) → MSv4:

| FITS-IDI | MSv4 | Notes |
| --- | --- | --- |
| `UV_DATA` | main `visibility` dataset | Row-major `(time, baseline)` → dense `(time, baseline_id, frequency, polarization)` |
| `ARRAY_GEOMETRY` + `ANTENNA` | `antenna_xds` | MSv4 already merges antenna and feed, as FITS-IDI does |
| `SOURCE` | `field_and_source_*_xds` | |
| `FREQUENCY` | `frequency` coordinate attributes | One MSv4 per spectral window; `spectral_window_name` from it |
| `FLAG` | `FLAG` variable | FITS-IDI flags are ranges to apply, not a per-datum array |
| `GAIN_CURVE` | `gain_curve_xds` | |
| `SYSTEM_TEMPERATURE` | `system_calibration_xds` | |
| `PHASE-CAL` | `phase_calibration_xds` | |
| `WEATHER` | `weather_xds` | Not handled by Gustavo; new work |
| `INTERFEROMETER_MODEL` | — | MSv4 lists this as future schema work; drop with a warning |
| `BASELINE`, `CALIBRATION`, `BANDPASS`, `MODEL_COMPS` | — | No MSv4 home; drop with a warning |

## The seven real problems

Ordered by how much design they need, not by code volume.

1. **Partitioning.** One MSv4 covers a single spectral window, polarization
   setup, observation mode, processor and beam. FITS-IDI's `BAND` axis — counted
   by `NO_BAND` — is what MSv4 calls a **spectral window**; it is not a receiver
   band, despite the `Qband`/`Wband` in the filenames, which name the receiver.
   One `UV_DATA` holds all `NO_BAND` spectral windows in a single array, plus
   many sources and scans, so splitting by spectral window is mandatory.
   Splitting further is a choice measured in the CHUNK-001 addendum.

2. **Scans.** MSv4 requires a `scan_name` coordinate. FITS-IDI has no scan
   column — unlike MSv2's `SCAN_NUMBER` — so scans must be *inferred* from
   source changes and time gaps. That is a physical judgement baked into the
   output, and the gap threshold has to be a documented, overridable parameter.
   Gustavo's `_idi_segment_scans` is the starting point.

3. **Metadata MSv4 requires that FITS-IDI does not carry.** `scan_intents`,
   `spectral_window_intents`, `observation_info.project_UID` and
   `release_date` are all non-optional in MSv4 and absent from FITS-IDI. Each
   needs a defensible default, a keyword to override it, and a docstring saying
   what was invented. This is the single biggest correctness-of-provenance risk:
   a converter that silently fabricates `intents` produces data that *looks*
   authoritative and is not.

4. **Weights.** Confirmed by CHUNK-001 to be worse than assumed. The test file
   has `WEIGHTYP = "CORRELAT"`, so weights are correlation-coefficient
   normalized rather than `1/σ²` — the conversion is required, not a corner
   case. And the `WEIGHT` column is `NO_BAND × NO_STKD`, *not* per channel,
   while MSv4's `WEIGHT` has the full science-array shape, so it must also be
   broadcast across `frequency`. Gustavo's `weight_mode` handling
   (`:auto`/`:validity`/`:radiometer`) encodes hard-won knowledge here and
   should be ported rather than re-derived.

5. **Time.** FITS-IDI works in days/hours relative to `RDATE`; MSv4 wants
   seconds since an epoch with an explicit `scale` and `format`. Writing
   `format: "unix"`, `scale: "utc"` means a real UTC→Unix conversion, and for
   VLBI the precision matters. Consider `format: "mjd"` instead, which avoids
   the epoch shift entirely — see CHUNK-004.

6. **UVW and the other unit conversions.** `UU/VV/WW---SIN` carry
   `TUNIT = SECONDS` (verified), so multiply by `c` for MSv4's metres. The sign
   convention (antenna1→antenna2) is the one thing headers cannot settle and
   must be checked against real `BASELINE` values. Several other conversions
   ride along: `RAEPO`/`DECEPO` and `WIND_DIRECTION` in degrees → radians,
   `TEMPERATURE`/`DEWPOINT` in centigrade → kelvin, `PRESSURE` in millibars →
   pascals. Also settle the frame attribute — MSv4 permits `icrs`/`fk5`/`APP`.

7. **Dense regridding, at a size that will not fit in memory.** FITS-IDI rows
   are `(time, baseline)` pairs with dropouts; MSv4 is a dense rectangular grid
   padded with `NaN`. The plan originally assumed the converter could be
   simpler than Gustavo's lazy reader because it writes once — CHUNK-001 shows
   otherwise, at ~10 GB per band partition. See CHUNK-014.

## Test strategy

Two real VLBA files now sit under `testdata/` (CHUNK-001), but at 24.73 GB the
usable one is an integration fixture, not a unit-test one. The synthetic route
is still needed: a minimal, test-only FITS-IDI writer paralleling
`test/synthetic.jl`, which already builds a conforming MSv4 from nothing.
Gustavo takes the same approach, writing fixtures with its own `write_fitsidi`.

Three layers, mirroring what is already in place:

1. **Synthetic round-trip.** Build a small FITS-IDI in the test suite, convert,
   and assert `isempty(check(ps))` — the converter's output must pass the
   existing schema checker with no issues.
2. **Value round-trip.** Assert visibilities, UVW, weights, flags, antenna
   positions and frequencies survive with the values put in, including the
   axis-order reversal.
3. **Real file, by hand.** Convert the Q-band file and run `xradio`'s
   `check_datatree` on the result, as was done for the writer. Because it is
   24.73 GB, prefer converting a row range or a single spectral window; a full conversion
   is a deliberate long-running job, not part of any suite.

## Chunks

### CHUNK-001: preflight
- **Kind**: `preflight`
- **Description**: Establish the baseline and, critically, **obtain a real FITS-IDI file**. Record its HDU list and headers here — every later chunk is guessing without it. Verify the `UU---SIN` units, the UVW sign convention, and whether `NO_BAND > 1`.
- **Status**: `complete`
- **Notes**: Baseline is green. Tree clean at `bfc45d5`, version 0.1.0 (untracked: `zarr-mwes/`, `testdata/` is git-ignored). Tests: `julia --project=test test/runtests.jl` — 1532 passing, 0 failures. `Test.detect_ambiguities(XRadio; recursive = true)`: **0**. Measured on Julia 1.11.9 from the shell, **not** via the Julia MCP session — `kaimon` failed to connect for the whole of the session that wrote this plan, and MCP tools register only at session start. Check it is up before starting, since the skill assumes it.

  Two VLBA BT164A files under `testdata/` (git-ignored; `testdata/fetch-all.sh` re-fetches, `testdata/SHA1SUMS` verifies). Surveyed with the **dev checkout** of FITSFiles at `/home/ptiede/.julia/dev/FITSFiles`, which reads them without complaint. Note its API is `cards["KEY"]`/`haskey`, not the `card_value` helper Gustavo defines locally.

  **Q-band — `VLBA_BT164A_bt164aQband_...idifits`, 24.73 GB, complete.** The full-coverage case: every table the mapping needs is present.

  | HDU | rows | | HDU | rows |
  | --- | --- | --- | --- | --- |
  | ARRAY_GEOMETRY | 10 | | UV_DATA | 752,160 |
  | SOURCE | 3 | | FLAG | 7,245 |
  | ANTENNA | 10 | | SYSTEM_TEMPERATURE | 2,043 |
  | FREQUENCY | 1 | | PHASE-CAL | 2,368 |
  | INTERFEROMETER_MODEL | 1,265 | | WEATHER | 899 |
  | MODEL_COMPS | 1,265 | | GAIN_CURVE | 10 |
  | CALC | 5 | | | |

  `UV_DATA`: `NO_STKD=4`, `STK_1=-1` (RR, LL, RL, LR — circular, as expected for VLBA), `NO_BAND=4`, `NO_CHAN=256`, `REF_FREQ=42.912 GHz`, `CHAN_BW=500 kHz`, `WEIGHTYP="CORRELAT"`, `DATE-OBS=2025-12-08`, `TELESCOP=VLBA`, `OBSCODE=BT164A`. Columns: `UU/VV/WW---SIN [SECONDS]`, `DATE [DAYS]`, `TIME [DAYS]`, `BASELINE`, `FILTER`, `SOURCE`, `FREQID`, `INTTIM [SECONDS]`, `WEIGHT 16E`, `GATEID 0J` (empty), `FLUX 8192E [UNCALIB]`. Data-matrix axes `CTYPE1..4 = COMPLEX, STOKES (CRVAL=-1, CDELT=-1), FREQ, BAND`.

  `ARRAY_GEOMETRY`: `FRAME=GEOCENTRIC`, `TIMSYS=UTC`, `RDATE=2025-12-08`; columns `ANNAME`, `STABXYZ [METERS]`, `NOSTA`, `MNTSTA`, `STAXOF [METERS]`, `DIAMETER [METERS]` — a direct fit to `antenna_xds`.
  `FREQUENCY`: one row, `BANDFREQ`, `CH_WIDTH`, `TOTAL_BANDWIDTH`, `SIDEBAND`, `BB_CHAN`, each of length `NO_BAND`.
  `SOURCE`: `RAEPO`/`DECEPO [DEGREES]`, `EQUINOX`, `RAAPP`/`DECAPP`, fluxes, velocities.
  `WEATHER`: `TEMPERATURE [CENTIGRADE]`, `PRESSURE [MILLIBARS]`, `DEWPOINT [CENTIGRADE]`, `WIND_VELOCITY [M/SEC]`, `WIND_DIRECTION [DEGREES]`, `WVR_H2O`, `IONOS_ELECTRON`.
  `SYSTEM_TEMPERATURE`: `TSYS_1/2 [K]`, `TANT_1/2 [K]` only. `GAIN_CURVE`: `TYPE`, `NTERM`, `X_VAL`, `Y_VAL`, `GAIN`, `SENS [K/JY]`. `PHASE-CAL`: `CABLE_CAL [SECONDS]`, `PC_FREQ [Hz]`, `PC_REAL`, `PC_IMAG`, `PC_RATE`. `FLAG`: `SOURCE_ID`, `ARRAY`, `ANTS`, `FREQID`, `TIMERANG [DAYS]`, `BANDS`, `CHANS`, `PFLAGS`, `REASON`, `SEVERITY`.
  The `_1`/`_2` suffixes throughout are **per receptor** (R and L) and map onto MSv4's `receptor_label`.

  **W-band — `VLBA_BT164A_bt164aWbandrecor_...idifits`: INCOMPLETE.** Its `UV_DATA` header declares 47.43 GB and only 3.09 GB is on disk; a `wget -c` is still running. `NO_BAND=16`, `NO_CHAN=64`, 1,434,170 rows, and it carries **only** ARRAY_GEOMETRY, SOURCE, ANTENNA, FREQUENCY, INTERFEROMETER_MODEL, CALC, MODEL_COMPS — no FLAG, SYSTEM_TEMPERATURE, PHASE-CAL, WEATHER or GAIN_CURVE. Once complete it is the useful *minimal* counterpart, and its 16 spectral windows stress partitioning harder than Q-band's 4. Re-verify against `SHA1SUMS` before relying on it.

  **Answers to the three questions this chunk existed to settle:**
  1. `UU/VV/WW---SIN` carry `TUNIT = SECONDS`. Multiply by `c` for MSv4 metres — confirmed, not assumed.
  2. `NO_BAND` is 4 and 16 — these count **spectral windows**, not receiver bands. Partitioning by spectral window is **mandatory**, not hypothetical.
  3. Sign convention still unverified — it needs actual `BASELINE` values and antenna positions, not headers. Carried into CHUNK-007.


  **Addendum — scan and baseline structure, measured.** Read the scalar columns
  of `UV_DATA` (all within the first 48 bytes of each 32,880-byte row, so ~23 s
  for all 752,160 rows) and segmented scans on a source change or a time gap
  greater than ten median steps:

  | | |
  | --- | --- |
  | unique times | 15,120 over 7.64 h, `INTTIM` 0.5 s |
  | baselines | 55 (10 antennas: 45 cross + 10 auto) |
  | sources | 3 |
  | scans | **36**, by time gap alone *and* by gap-plus-source |
  | times per scan | min 240, median 600, max 600 (5-minute scans) |

  Three results that bear on partitioning:

  - **Source changes always coincide with scan boundaries here**, so `:scan`
    subsumes `:source` for this file — both give 36. Splitting by source as
    well as by scan would add nothing.
  - **Every source was seen by all 55 baselines.** The hypothesis that
    telescopes cannot all see the same sources does not hold for this VLBA
    track, though it remains true of wider arrays.
  - **Baselines per scan range from 10 to 55.** The heterogeneity is
    per *scan*, not per source — antennas drop out within the track, not by
    target.

  Partition sizes for `VISIBILITY` in one spectral window:

  | scheme | per partition | fill | nodes |
  | --- | --- | --- | --- |
  | spectral window only | 6.81 GB | 90.4% | 4 |
  | spw + source | 4.87 GB | — | 12 |
  | **spw + scan** | **270 MB** | **100.0%** | 144 |

  Splitting by scan makes every partition **exactly dense** — the 9.6% of NaN
  padding disappears, saving ~2.6 GB across the track — and drops peak working
  memory 25-fold, which very largely dissolves CHUNK-014.

  **Three findings that were not anticipated:**
  - **`WEIGHT` is `16E` = `NO_BAND × NO_STKD`, not per channel.** MSv4 `WEIGHT` has the full science-array shape, so the weight must be broadcast across `frequency`.
  - **`WEIGHTYP = "CORRELAT"`.** Weights are correlation-coefficient normalized, not `1/σ²`. Converting them is required, not optional — Gustavo's `weight_mode = :radiometer` path is what this needs.
  - **`FLUX` is in `UNCALIB` units,** so `VISIBILITY.units` must not claim `Jy`.

### CHUNK-002: decide-partitioning
- **Kind**: `decide`
- **Description**: How does the user control the split from one FITS-IDI into N Measurement Sets? Questions: (a) what is the default — spectral window alone, or spectral window plus scan? (b) what is the keyword — `partition = (:scan,)` mirroring xradio's `partition_scheme`, or a more Julian predicate? (c) how are the resulting MSv4 nodes named, given MSv4 names must be unique within a store and xradio uses `<basename>_<zero-padded index>`?
  **Recommendation, revised after the CHUNK-001 addendum measured it**: default to **spectral window + scan**. On the test file that gives 144 partitions of 270 MB at 100% fill, against 4 partitions of 6.81 GB at 90.4% fill for spectral window alone. It is denser, 25x smaller in peak memory, and each partition's baseline axis then describes antennas that actually observed together — baselines per scan run from 10 to 55. `:source` stays available but is subsumed by `:scan` on this data and should not be a default. Names `<file basename>_<NN>`, matching xradio so a mixed store stays coherent.
  **Cost to weigh**: 144 nodes means ~4,300 metadata files, and this package does not yet write consolidated metadata — Python warns about its absence when opening our stores today. At 4 nodes that is irrelevant; at 144 it governs how long `open` takes, especially over a network. `Zarr.consolidate_metadata` exists (unexported) for both Zarr versions, so this is a small addition to the writer, but it should land with this decision rather than after it.
- **Status**: `complete`
- **Depends on**: CHUNK-001
- **Notes**: spectral window + scan, `partition = (:scan,)`, names
  `<basename>_<NN>`. Recorded under Decisions. Consolidated metadata (CHUNK-015) is now a
  prerequisite of CHUNK-012 rather than an independent improvement.

### CHUNK-003: decide-invented-metadata
- **Kind**: `decide`
- **Description**: Settle what the converter puts in the MSv4 fields FITS-IDI cannot supply, and how loudly.
  **Largely answered by reading what the reference implementation does** (`_utils/_msv2/msv4_info_dicts.py`, `conversion.py`) and what the FITS-IDI primary header actually carries. The set of genuinely invented fields is far smaller than this chunk assumed:

  | MSv4 field | FITS-IDI source | xradio's own MSv2 behaviour |
  | --- | --- | --- |
  | `observation_info.observer` | primary `OBSERVER` (`"BT164A"` — the obs code, not a person; VLBA convention) | `OBSERVATION.OBSERVER` |
  | `observation_info.project_UID` | `UV_DATA.OBSCODE` (`"BT164A"`) | `OBSERVATION.PROJECT`, renamed — xradio does **not** synthesize a UID either |
  | `observation_info.observing_log` | `HISTORY` cards (correlator, version, per-job spans) | `OBSERVATION.LOG` |
  | `processor_info.type` | `"CORRELATOR"` | PROCESSOR table |
  | `processor_info.sub_type` | primary `CORRELAT` (`"DIFX"`) | PROCESSOR table |
  | `spectral_window_name` | `"spw_<n>"` | same `"spw_<id>"` fallback when `SPECTRAL_WINDOW.NAME` is empty |
  | `observation_info.release_date` | **none** | `OBSERVATION.RELEASE_DATE` |
  | `scan_intents` | **none** — FITS-IDI has no `STATE`/`OBS_MODE` analogue | `STATE.OBS_MODE`, comma-split |
  | `spectral_window_intents` | **none** | **hardcoded `["UNSPECIFIED"]`**, even from MSv2 (`conversion.py:543`) |

  So only three fields lack a source, and each has a precedent rather than needing invention:
  - `spectral_window_intents` → `["UNSPECIFIED"]`, which is exactly what xradio writes from MSv2. Not a compromise; it is the reference behaviour.
  - `scan_intents` → `["UNSPECIFIED"]`, matching xradio's own accessor default (`measurement_set_xdt.py:233`).
  - `release_date` → `"1858-11-17T00:00:00.000"`, the MJD-zero "unset" epoch that real converted MSv4 data already carries (the ALMA `Antennae_North` fixture has exactly this).

  Remaining questions: (a) should each of the three be overridable by keyword? (b) is one summarizing warning right, or should a sentinel be silent given it matches the reference implementation? (c) does the source file name go into `data_groups.base.description` as provenance?
  **Recommendation**: keywords for all three; warn only for `release_date` (the one where a sentinel could be mistaken for real policy) and stay silent on the two intents, since `["UNSPECIFIED"]` is the standard's normal content rather than a shortfall; record the source file in `data_groups.base.description`.
  **Aside worth reporting upstream**: `conversion.py:1084` reads `partition_info["OBS_MODE"][0]` unconditionally, but `OBS_MODE` is only added when a non-empty `STATE` table exists (`partition_queries.py:124-142`). An MSv2 without a `STATE` table — common for VLBI — appears to raise `KeyError` there. Worth confirming before filing.
- **Status**: `complete`
- **Depends on**: CHUNK-001
- **Notes**: keywords for all three sentinels; warn on `release_date`
  only. Recorded under Decisions. One thing this chunk did not anticipate: the primary
  header and `UV_DATA` carry **different** `OBSERVER` values (`"BT164A"` vs `"PLUTO"`).
  Resolved in favour of the primary header.

### CHUNK-004: decide-time-and-frames
- **Kind**: `decide`
- **Description**: Settle the time and coordinate-frame conversions. Questions: (a) `format = "unix"` (needs a real UTC→Unix conversion, and leap seconds) or `format = "mjd"` (FITS-IDI is already MJD-adjacent via `RDATE`, avoiding an epoch shift and a precision loss on ~1e9-second values held in `Float64`)? MSv4 permits both. (b) which `frame` for `UVW` and which for the source directions, given FITS-IDI's `EQUINOX`/`RADECSYS`? (c) does any of this need a time-scale dependency (AstroTime.jl, ERFA.jl), or is it arithmetic on values FITS-IDI already provides? **Recommendation**: `mjd`, which keeps VLBI-relevant precision and needs no epoch arithmetic; revisit only if a consumer demands `unix`.
- **Status**: `complete`
- **Depends on**: CHUNK-001
- **Notes**: `format = "unix"`, **reversing** this chunk's own `mjd`
  recommendation — `units = "s"` is mandatory, which forces MJD to be written in seconds and
  destroys both arguments for it. Frames map `EQUINOX` literally, `"J2000"` → `"fk5"`.
  `RAAPP`/`DECAPP` are copies of the mean places and must not be used. No new dependency.
  Full reasoning and the astropy measurement are under Decisions.

### CHUNK-005: implement-fitsidi-reader
- **Kind**: `implement`
- **Description**: A `XRadioFITSFilesExt` package extension on `FITSFiles` (weak dependency, as Gustavo does) holding what it takes to read a FITS-IDI file: access primitives for cards and columns, the arithmetic the format states rather than stores (the STOKES axis, `BASELINE` packing, channel frequencies, scan inference), and a bulk reader for `UV_DATA`'s per-row columns. Port from Gustavo's `_build_idi_*` helpers. **No MSv4 types and no intermediate data model** — the conversion in CHUNK-007 builds MSv4 datasets straight from the tables.
- **Status**: `complete`
- **Depends on**: CHUNK-001
- **Notes**: Landed as `ext/XRadioFITSFilesExt/fitsidi.jl` (346 lines) — access
  primitives, the format's own arithmetic, and one bulk reader for `UV_DATA`'s
  per-row columns. **No intermediate data model**; see `Decisions` for that, for
  where the code lives, for the polarization representation, and for which
  FITSFiles to depend on.

  Verified on the Q-band file: 10 stations, 4 bands x 256 channels x 4 pols,
  `WEIGHTYP = 'CORRELAT'`, `SORT = 'T*'`, 3 sources (3C273, M87, 3C279), 752160
  `UV_DATA` rows read in ~23 s. `scan_ranges` recovers **36** scans, which is
  exactly the number of `HISTORY DIFXJOB` entries in the primary header — an
  independent confirmation of the gap heuristic, and the test asserts the two
  agree. 36 scans x 4 bands = the 144 partitions CHUNK-002 predicted.

  `UV_DATA` row layout on this file: 752160 rows of 32880 bytes; UU/VV/WW, DATE,
  TIME, BASELINE, FILTER, SOURCE, FREQID, INTTIM all inside the first 48 bytes,
  then `WEIGHT` (16 = NO_STKD x NO_BAND), a zero-length `GATEID`, and `FLUX`
  (8192 = 2 x 4 x 256 x 4). So `read_uv_rows` makes one strided pass collecting
  48-byte prefixes into a single 36 MB buffer and decodes columns out of memory;
  column-at-a-time access would traverse 24.73 GB once per column.

  Not carried over from the discarded model, and now CHUNK-007's to write when it
  needs them: the `FLAG` table parse (CHUNK-009) and the `MNTSTA` mount-code
  lookup, which survives as `IDI_MOUNTS`.

### CHUNK-006: implement-test-fitsidi-writer
- **Kind**: `implement`
- **Description**: A minimal FITS-IDI writer under `test/`, enough to build fixtures with known values: `ARRAY_GEOMETRY`, `ANTENNA`, `SOURCE`, `FREQUENCY`, `UV_DATA`, and optionally `FLAG`. Not public API and not in `src/`. Model on Gustavo's `fitsidi_write.jl`, keeping only what a fixture needs. Parallels the existing `test/synthetic.jl`.
- **Status**: `not-started`
- **Depends on**: CHUNK-005
- **Notes**: **Load-bearing, not a convenience.** With no intermediate model between
  the file and the MSv4 arrays there is nothing else to assert on, so a fixture
  `.idifits` is the only way to test the conversion without the 24 GB validation
  file. Test-only throughout: under `test/`, never in `src/`, never public API.

### CHUNK-007: implement-main-dataset-conversion
- **Kind**: `implement`
- **Description**: Rows → the dense `(time, baseline_id, frequency, polarization)` grid, in XRadio's reversed Julia axis order, with `NaN` padding for absent baselines. Produce `VISIBILITY`, `FLAG`, `WEIGHT`, `UVW` and the coordinates `time`, `baseline_id`, `frequency`, `polarization`, `field_name`, `scan_name`, `baseline_antenna1_name`, `baseline_antenna2_name`. Apply the CHUNK-002 partitioning, the CHUNK-003 defaults, the CHUNK-004 time convention, the weight semantics from `WEIGHTYP`, and the UVW seconds→metres conversion. Gate: the result passes `check` with no issues.
- **Status**: `not-started`
- **Depends on**: CHUNK-002, CHUNK-003, CHUNK-004, CHUNK-005, CHUNK-006
- **Notes**: `fitsidi2msv4` writes each partition to the destination store as it is
  built rather than assembling a `ProcessingSet` in memory, so the `check` gate runs
  against the store read back from disk. See `Decisions`.

### CHUNK-008: implement-antenna-and-field-source
- **Kind**: `implement`
- **Description**: `antenna_xds` from ARRAY_GEOMETRY + ANTENNA (positions, mounts, dish diameters, receptor labels, polarization types), and `field_and_source_*_xds` from SOURCE. These two are the sub-datasets MSv4 effectively requires — `check` reports a missing `antenna_xds`, and every data group names a `field_and_source`.
- **Status**: `not-started`
- **Depends on**: CHUNK-007
- **Notes**:

### CHUNK-009: implement-flag-table
- **Kind**: `implement`
- **Description**: FITS-IDI `FLAG` rows are *ranges* (`SOURCE_ID`, `ANTS` with 0 as wildcard, `FREQID`, `TIMERANG` in days relative to `RDATE`, `BANDS`, `CHANS`, `PFLAGS`) to be applied, not a per-datum array. Expand them onto the dense `FLAG` grid. Port `_build_idi_flags`/`_filter_flags_for_leaf`. Note that Gustavo's `write_fitsidi` emits no FLAG table, so its tests append one by hand with `FITSFiles: HDU, Bintable, Card` (`test/test_fitsidi.jl:612-625`) — that recipe is directly reusable for CHUNK-006 and is the only worked example of building this table in Julia.
- **Status**: `not-started`
- **Depends on**: CHUNK-007
- **Notes**:

### CHUNK-010: implement-calibration-sub-datasets
- **Kind**: `implement`
- **Description**: `system_calibration_xds` (SYSTEM_TEMPERATURE), `gain_curve_xds` (GAIN_CURVE), `phase_calibration_xds` (PHASE-CAL), `weather_xds` (WEATHER). Each optional and skipped when its table is absent. This is where FITS-IDI→MSv4 earns its keep over MSv2→MSv4 for VLBI, and where XRadio's schema tables get their first real exercise beyond the correlated datasets.
- **Status**: `not-started`
- **Depends on**: CHUNK-008
- **Notes**:

### CHUNK-011: implement-public-api
- **Kind**: `implement`
- **Description**: The user-facing entry point. Shape: `read(ProcessingSet, path; format = :fitsidi, ...)` reusing the existing verb, or a distinct `convert_fitsidi(in, out)`. Settle against the `open`/`read`/`write` convention already established. Include the keywords from CHUNK-002/003, and make `write` the natural next call so `read(...) |> ps -> write(out, ps)` is the whole workflow.
- **Status**: `not-started`
- **Depends on**: CHUNK-010
- **Notes**: Three export-surface changes already landed, ahead of this chunk and at
  the author's direction, so the shape this chunk settles is narrower than it was:
  - `check`, `validate`, `SchemaIssue` and `SchemaError` are **public but not
    exported**. Conformance checking is not what most work with a store involves,
    so it does not earn a name in the caller's namespace. The README shows
    `using XRadio: check, validate`.
  - DimensionalData is **re-exported** (`@reexport`, hence the `Reexport`
    dependency), because a `MeasurementSet` is a `DimTree` and nothing can be
    selected or sliced without `At`, `Near`, `dims` and `lookup`.
  - The dimension type is `UVWLabel`, not `UvwLabel`. The MSv4 dimension name
    `uvw_label` is unchanged; only the Julia spelling moved.

  `BaselineId` and `ElementId` still spell their acronym `Id` rather than `ID`,
  which is the same inconsistency; settle both here.

### CHUNK-012: implement-real-file-validation
- **Kind**: `implement`
- **Description**: Convert the real file from CHUNK-001, run `xradio`'s `check_datatree` on the output, and record the result in the README's conformance section alongside the existing writer evidence. Use the venv recipe in the `msv4-test-fixtures` memory. Any issue found here is a converter bug, not a fixture artefact.
- **Status**: `not-started`
- **Depends on**: CHUNK-011, CHUNK-015
- **Notes**: CHUNK-015 is a dependency because scan-level partitioning produces ~144 nodes;
  without consolidated metadata, opening the result in Python is slow enough to obscure
  whether `check_datatree` is reporting a conformance problem or a latency one.

### CHUNK-013: implement-docs
- **Kind**: `implement`
- **Description**: README section on the converter: the table mapping, every invented field and its default, the scan-inference threshold, and the tables dropped with a warning. A user must be able to tell from the docs what in the output came from their file and what the converter made up.
- **Status**: `not-started`
- **Depends on**: CHUNK-012
- **Notes**:

### CHUNK-014: decide-streaming-architecture
- **Kind**: `decide`
- **Description**: **Raised by CHUNK-001 and not anticipated when this plan was written.** The Q-band file is 24.73 GB of `UV_DATA`, and splitting it by spectral window gives four partitions of roughly 6.2 GB of `VISIBILITY` alone (13,676 times x 55 baselines x 256 channels x 4 polarizations x 8 bytes), plus ~3 GB of `WEIGHT` and ~0.8 GB of `FLAG` each. `write(path, ps)` as it stands takes a fully materialized `ProcessingSet`, so the obvious implementation would need ~10 GB resident per partition and would not run on an ordinary machine.
  Questions: (a) does the converter bypass the in-memory `ProcessingSet` and write Zarr chunk by chunk as it reads FITS-IDI rows? (b) or does `ProcessingSet` learn to hold lazy arrays, with `write` streaming them — which would also let `open(ProcessingSet, ...; lazy = true)` round-trip without materializing? (c) either way, what is the default `chunks` for a converted store, given the current default of one chunk per array is exactly wrong at this size?
  **Recommendation**: (b). The reader already produces lazy `Zarr.ZArray`-backed layers, so a `write` that streams whatever array it is handed makes the lazy read path useful for the first time and keeps one write implementation rather than two. A time-chunked default (say 128 times x all baselines x 64 channels x all polarizations) follows from how the data is read and how it is later used.
- **Status**: `complete`
- **Depends on**: CHUNK-001, CHUNK-002
- **Notes**: the CHUNK-001 addendum found that partitioning by scan gives 270 MB partitions rather than 6.81 GB ones. If CHUNK-002 adopts that default, this chunk stops being a blocker and becomes an optimization — a 270 MB partition materializes comfortably. Keep it open, because a correlator producing one long scan would bring the problem straight back, but it need not gate CHUNK-007 in that case.

  **Decided: deferred.** The scan-level partitioning adopted in
  CHUNK-002 confirms the reading above — 270 MB partitions materialize comfortably, so this
  does not gate CHUNK-007. Left open as an optimization against single-long-scan data.
  Recorded under Decisions.

### CHUNK-015: implement-consolidated-metadata
- **Kind**: `implement`
- **Description**: Write Zarr consolidated metadata when closing a store, and use it on `open`. Independent of FITS-IDI — it improves every store this package writes — but the converter is what makes it matter, since per-scan partitioning produces on the order of 144 nodes and several thousand metadata files. xarray already warns when opening our stores that it fell back to non-consolidated metadata.

  **How the reference implementation does it** (`convert_msv2_to_processing_set.py:256-264`): a single finalize step after every partition is written, not per partition —

  ```python
  root_group = zarr.open(out_file, mode="r+", zarr_format=ZARR_FORMAT)
  root_group.attrs["type"] = "processing_set"
  zarr.consolidate_metadata(root_group.store)
  ```

  Consolidated metadata is **not** part of the Zarr v3 specification. xradio writes it there regardless and silences the resulting warning in `src/xradio/__init__.py`, commented "XRADIO intentionally writes/reads consolidated metadata" — a deliberate trade of spec purity for open latency. Their read side is inconsistent: `open_processing_set` passes no `consolidated=` argument and takes xarray's try-then-warn default, while `load_processing_set` passes `consolidated=False` outright (their own notes call that a performance bug over S3). We should not copy that inconsistency.

  **What `Zarr.jl` can actually do**, measured:

  | | write | read back |
  | --- | --- | --- |
  | v2 | `consolidate_metadata(g)` | root attrs, children and node attrs all intact |
  | v3 | `consolidate_metadata(g.storage, g.path, ZarrFormat(3))` writes a correct `consolidated_metadata` into the root `zarr.json` | **empty group** — no attributes, no children |

  So v3 consolidated metadata can be written for zarr-python's benefit but must not be read back by us. Also note `consolidate_metadata(::ZGroup)` does not dispatch on the group's format: on a v3 store it writes a v2-style `.zmetadata` and leaves `zarr.json` alone. Reproduced in `zarr-mwes/4_consolidated_metadata_v3.jl`.

  Consequences for the implementation: consolidate as a finalize step in `write`; on `open`, use consolidated metadata for v2 and fall back for v3 until Zarr.jl is fixed; do not fail when it is absent, since existing stores lack it.
- **Status**: `not-started`
- **Depends on**: CHUNK-002
- **Notes**:

## Open risks

- **Scale.** 24.73 GB in. Scan-level partitioning (CHUNK-002) cuts this to 270 MB
  per partition, so nothing needs to stream and CHUNK-014 no longer gates the main
  conversion. The constraint returns for any correlator that emits one long scan. A unit
  test still cannot use these files; the synthetic writer (CHUNK-006) is not optional.
- **The UVW sign convention is still unverified.** It is the one CHUNK-001
  question headers could not answer, and getting it wrong is silent — the
  output would pass the schema checker and be conjugated. Check it against
  actual `BASELINE` values and `STABXYZ` positions in CHUNK-007.
- **Unit conversions are pervasive and individually trivial**, which is what
  makes them easy to miss: seconds→metres (UVW), degrees→radians (`RAEPO`,
  `DECEPO`, `WIND_DIRECTION`), centigrade→kelvin (`TEMPERATURE`, `DEWPOINT`),
  millibars→pascals (`PRESSURE`). MSv4 states units in every array's
  attributes, so `check` catches a wrong *declaration* but never a wrong
  *value*. Each conversion needs a test with a known number.
- **A mixed-basis array has no correct MSv4 `polarization` coordinate.** The
  reader now represents this honestly — generic products, per-antenna basis —
  but CHUNK-007 still has to write a single array-wide coordinate to conform,
  and for a heterogeneous array any value it picks is wrong for some baseline.
  Decide explicitly there; the validation files are both homogeneously
  circular, so nothing in the test suite will catch a silent choice.
- **Scan inference and invented metadata are judgement, not translation.** They
  belong in the docs and in the store's provenance, not only in the code.
- **Trailing whitespace and NUL bytes in FITS strings.** Untrimmed, they reach MSv4
  coordinate labels and comparison sites alike, and every symptom is silent. Recorded in
  CHUNK-005's Notes.
- **The FITSFiles dependency is not satisfiable from the registry.** The
  extension needs the lazy-reading API that exists only on `ptiede-uvfixes`;
  the registered `0.3.2` shares the version number but not the code, and its
  `hdu.data` materializes an entire HDU per access. Until those commits are
  merged and released, `test/Manifest.toml` must pin the branch and XRadio
  cannot be registered. Recorded under Decisions.
- **FLAG expansion is the fiddliest table.** Six independent range axes with
  wildcards, and no fixture generator writes one — Gustavo's tests hand-build
  the HDU. Budget for it accordingly.

## Open Questions
<!-- Problems noticed in passing. Recorded rather than fixed in-flight. -->

- `intents(ms)` in `src/accessors.jl` reads `observation_info["intents"]`, which
  the MSv4 schema does not define — the standard puts intents on the `scan_name`
  coordinate (`scan_intents`) and on `frequency` (`spectral_window_intents`).
  It works on older fixtures because xradio 0.0.58 wrote that key, but current
  xradio has the line commented out. The accessor should read `scan_intents`
  and treat `observation_info["intents"]` as a legacy fallback.
- `check` does not descend into attributes that are themselves labelled arrays,
  so a malformed `reference_frequency` or `integration_time` passes. Noted in
  the README as a known limitation.

## Session ledger

- 2026-09-05 CHUNK-001 (preflight) → next: CHUNK-002, CHUNK-003, CHUNK-004, CHUNK-014 (all `decide`, batchable), then CHUNK-005
- 2026-09-06 CHUNK-002, CHUNK-003, CHUNK-004, CHUNK-014 (all `decide`) → next: CHUNK-005
- 2026-09-06 CHUNK-005 (implement-fitsidi-reader) → next: CHUNK-006
