# FITS-IDI → MSv4 converter plan

<!-- Preserve chunk IDs and status values. -->

## Metadata
- **Kind**: `feature`
- **Package**: XRadio
- **Written**: 2026-09-05
- **Current version**: 0.1.0 (unregistered, API breaks freely)
- **Goal**: read a FITS-IDI file and write it as a conforming MSv4 processing set

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
- **Notes**: Two VLBA BT164A files under `testdata/` (git-ignored; `testdata/fetch-all.sh` re-fetches, `testdata/SHA1SUMS` verifies). Surveyed with the **dev checkout** of FITSFiles at `/home/ptiede/.julia/dev/FITSFiles`, which reads them without complaint. Note its API is `cards["KEY"]`/`haskey`, not the `card_value` helper Gustavo defines locally.

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
- **Status**: `todo`
- **Depends on**: CHUNK-001

### CHUNK-003: decide-invented-metadata
- **Kind**: `decide`
- **Description**: Settle what the converter puts in the MSv4 fields FITS-IDI cannot supply, and how loudly. Fields: `scan_intents`, `spectral_window_intents`, `observation_info.project_UID`, `release_date`, `processor_info.{type,sub_type}`, `data_groups.base.description`. Questions: (a) invent a documented placeholder, or refuse to convert without the user supplying it? (b) if invented, is a warning emitted per conversion, or is a record written into the store (e.g. `creator` or a `description` naming the source file and this package)? **Recommendation**: placeholders with a single summarizing warning listing every invented field, plus provenance in `data_groups.base.description` naming the FITS-IDI file — fail-fast is right for corrupt input, but refusing to convert a valid FITS-IDI because it lacks an ALMA project code would make the converter useless for its actual audience.
- **Status**: `todo`
- **Depends on**: CHUNK-001

### CHUNK-004: decide-time-and-frames
- **Kind**: `decide`
- **Description**: Settle the time and coordinate-frame conversions. Questions: (a) `format = "unix"` (needs a real UTC→Unix conversion, and leap seconds) or `format = "mjd"` (FITS-IDI is already MJD-adjacent via `RDATE`, avoiding an epoch shift and a precision loss on ~1e9-second values held in `Float64`)? MSv4 permits both. (b) which `frame` for `UVW` and which for the source directions, given FITS-IDI's `EQUINOX`/`RADECSYS`? (c) does any of this need a time-scale dependency (AstroTime.jl, ERFA.jl), or is it arithmetic on values FITS-IDI already provides? **Recommendation**: `mjd`, which keeps VLBI-relevant precision and needs no epoch arithmetic; revisit only if a consumer demands `unix`.
- **Status**: `todo`
- **Depends on**: CHUNK-001

### CHUNK-005: implement-fitsidi-reader
- **Kind**: `implement`
- **Description**: A `XRadioFITSFilesExt` package extension on `FITSFiles` (weak dependency, as Gustavo does) that parses the FITS-IDI tables into plain Julia structures — antennas, frequency setups, sources, flags, and a row iterator over `UV_DATA`. **No MSv4 types yet**: this chunk's output is intermediate, so it can be tested against a real file independently of the conversion. Port from Gustavo's `_build_idi_*` helpers. Handle the STOKES code permutation here.
- **Status**: `todo`
- **Depends on**: CHUNK-001

### CHUNK-006: implement-test-fitsidi-writer
- **Kind**: `implement`
- **Description**: A minimal FITS-IDI writer under `test/`, enough to build fixtures with known values: `ARRAY_GEOMETRY`, `ANTENNA`, `SOURCE`, `FREQUENCY`, `UV_DATA`, and optionally `FLAG`. Not public API and not in `src/`. Model on Gustavo's `fitsidi_write.jl`, keeping only what a fixture needs. Parallels the existing `test/synthetic.jl`.
- **Status**: `todo`
- **Depends on**: CHUNK-005

### CHUNK-007: implement-main-dataset-conversion
- **Kind**: `implement`
- **Description**: Rows → the dense `(time, baseline_id, frequency, polarization)` grid, in XRadio's reversed Julia axis order, with `NaN` padding for absent baselines. Produce `VISIBILITY`, `FLAG`, `WEIGHT`, `UVW` and the coordinates `time`, `baseline_id`, `frequency`, `polarization`, `field_name`, `scan_name`, `baseline_antenna1_name`, `baseline_antenna2_name`. Apply the CHUNK-002 partitioning, the CHUNK-003 defaults, the CHUNK-004 time convention, the weight semantics from `WEIGHTYP`, and the UVW seconds→metres conversion. Gate: the result passes `check` with no issues.
- **Status**: `todo`
- **Depends on**: CHUNK-002, CHUNK-003, CHUNK-004, CHUNK-005, CHUNK-006, CHUNK-014

### CHUNK-008: implement-antenna-and-field-source
- **Kind**: `implement`
- **Description**: `antenna_xds` from ARRAY_GEOMETRY + ANTENNA (positions, mounts, dish diameters, receptor labels, polarization types), and `field_and_source_*_xds` from SOURCE. These two are the sub-datasets MSv4 effectively requires — `check` reports a missing `antenna_xds`, and every data group names a `field_and_source`.
- **Status**: `todo`
- **Depends on**: CHUNK-007

### CHUNK-009: implement-flag-table
- **Kind**: `implement`
- **Description**: FITS-IDI `FLAG` rows are *ranges* (`SOURCE_ID`, `ANTS` with 0 as wildcard, `FREQID`, `TIMERANG` in days relative to `RDATE`, `BANDS`, `CHANS`, `PFLAGS`) to be applied, not a per-datum array. Expand them onto the dense `FLAG` grid. Port `_build_idi_flags`/`_filter_flags_for_leaf`. Note that Gustavo's `write_fitsidi` emits no FLAG table, so its tests append one by hand with `FITSFiles: HDU, Bintable, Card` (`test/test_fitsidi.jl:612-625`) — that recipe is directly reusable for CHUNK-006 and is the only worked example of building this table in Julia.
- **Status**: `todo`
- **Depends on**: CHUNK-007

### CHUNK-010: implement-calibration-sub-datasets
- **Kind**: `implement`
- **Description**: `system_calibration_xds` (SYSTEM_TEMPERATURE), `gain_curve_xds` (GAIN_CURVE), `phase_calibration_xds` (PHASE-CAL), `weather_xds` (WEATHER). Each optional and skipped when its table is absent. This is where FITS-IDI→MSv4 earns its keep over MSv2→MSv4 for VLBI, and where XRadio's schema tables get their first real exercise beyond the correlated datasets.
- **Status**: `todo`
- **Depends on**: CHUNK-008

### CHUNK-011: implement-public-api
- **Kind**: `implement`
- **Description**: The user-facing entry point. Shape: `read(ProcessingSet, path; format = :fitsidi, ...)` reusing the existing verb, or a distinct `convert_fitsidi(in, out)`. Settle against the `open`/`read`/`write` convention already established. Include the keywords from CHUNK-002/003, and make `write` the natural next call so `read(...) |> ps -> write(out, ps)` is the whole workflow.
- **Status**: `todo`
- **Depends on**: CHUNK-010

### CHUNK-012: implement-real-file-validation
- **Kind**: `implement`
- **Description**: Convert the real file from CHUNK-001, run `xradio`'s `check_datatree` on the output, and record the result in the README's conformance section alongside the existing writer evidence. Use the venv recipe in the `msv4-test-fixtures` memory. Any issue found here is a converter bug, not a fixture artefact.
- **Status**: `todo`
- **Depends on**: CHUNK-011

### CHUNK-013: implement-docs
- **Kind**: `implement`
- **Description**: README section on the converter: the table mapping, every invented field and its default, the scan-inference threshold, and the tables dropped with a warning. A user must be able to tell from the docs what in the output came from their file and what the converter made up.
- **Status**: `todo`
- **Depends on**: CHUNK-012

### CHUNK-015: implement-consolidated-metadata
- **Kind**: `implement`
- **Description**: Write Zarr consolidated metadata when closing a store, via `Zarr.consolidate_metadata` / `consolidate_metadata_v3`, and read it back on `open`. Independent of FITS-IDI — it improves every store this package writes — but the converter is what makes it matter, since a per-scan partitioning produces on the order of 144 nodes and several thousand metadata files. xarray currently warns when opening our stores that it fell back to non-consolidated metadata.
- **Status**: `todo`
- **Depends on**: CHUNK-002

### CHUNK-014: decide-streaming-architecture
- **Kind**: `decide`
- **Description**: **Raised by CHUNK-001 and not anticipated when this plan was written.** The Q-band file is 24.73 GB of `UV_DATA`, and splitting it by spectral window gives four partitions of roughly 6.2 GB of `VISIBILITY` alone (13,676 times x 55 baselines x 256 channels x 4 polarizations x 8 bytes), plus ~3 GB of `WEIGHT` and ~0.8 GB of `FLAG` each. `write(path, ps)` as it stands takes a fully materialized `ProcessingSet`, so the obvious implementation would need ~10 GB resident per partition and would not run on an ordinary machine.
  Questions: (a) does the converter bypass the in-memory `ProcessingSet` and write Zarr chunk by chunk as it reads FITS-IDI rows? (b) or does `ProcessingSet` learn to hold lazy arrays, with `write` streaming them — which would also let `open(ProcessingSet, ...; lazy = true)` round-trip without materializing? (c) either way, what is the default `chunks` for a converted store, given the current default of one chunk per array is exactly wrong at this size?
  **Recommendation**: (b). The reader already produces lazy `Zarr.ZArray`-backed layers, so a `write` that streams whatever array it is handed makes the lazy read path useful for the first time and keeps one write implementation rather than two. A time-chunked default (say 128 times x all baselines x 64 channels x all polarizations) follows from how the data is read and how it is later used.
- **Status**: `todo`
- **Depends on**: CHUNK-001, CHUNK-002
- **Note**: the CHUNK-001 addendum found that partitioning by scan gives 270 MB partitions rather than 6.81 GB ones. If CHUNK-002 adopts that default, this chunk stops being a blocker and becomes an optimization — a 270 MB partition materializes comfortably. Keep it open, because a correlator producing one long scan would bring the problem straight back, but it need not gate CHUNK-007 in that case.

## Open risks

- **Scale.** Settled by CHUNK-001 and now the dominant constraint: 24.73 GB in,
  ~10 GB per spectral-window partition out. Nothing in the current package streams, so
  CHUNK-014 gates the main conversion. A unit test cannot use these files;
  the synthetic writer (CHUNK-006) is not optional.
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
- **Scan inference and invented metadata are judgement, not translation.** They
  belong in the docs and in the store's provenance, not only in the code.
- **FLAG expansion is the fiddliest table.** Six independent range axes with
  wildcards, and no fixture generator writes one — Gustavo's tests hand-build
  the HDU. Budget for it accordingly.
