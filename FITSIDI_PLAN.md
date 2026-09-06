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
| `FREQUENCY` | `frequency` coordinate attributes | One MSv4 per band; `spectral_window_name` from the band |
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
   setup, observation mode, processor and beam. A FITS-IDI `UV_DATA` holds
   `NO_BAND` bands in one array and many sources and scans. Splitting by band is
   mandatory; splitting further (by source, by scan) is a user choice. Needs a
   keyword mirroring `xradio`'s `partition_scheme`.

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

4. **Weights.** FITS-IDI `WEIGHT` semantics depend on the `WEIGHTYP` header;
   MSv4 `WEIGHT` is strictly `1/σ²`. Gustavo's `weight_mode` handling
   (`:auto`/`:validity`/`:radiometer`) encodes hard-won knowledge here and
   should be ported rather than re-derived.

5. **Time.** FITS-IDI works in days/hours relative to `RDATE`; MSv4 wants
   seconds since an epoch with an explicit `scale` and `format`. Writing
   `format: "unix"`, `scale: "utc"` means a real UTC→Unix conversion, and for
   VLBI the precision matters. Consider `format: "mjd"` instead, which avoids
   the epoch shift entirely — see CHUNK-004.

6. **UVW.** FITS-IDI `UU---SIN`/`VV---SIN`/`WW---SIN` are in **seconds** (light
   travel time); MSv4 `UVW` is in **metres**. Multiply by `c`. Also settle the
   frame attribute — FITS-IDI is typically apparent-of-date, MSv4 permits
   `icrs`/`fk5`/`APP`. Verify the sign convention (antenna1→antenna2) against a
   real file rather than assuming.

7. **Dense regridding.** FITS-IDI rows are `(time, baseline)` pairs with
   dropouts; MSv4 is a dense rectangular grid padded with `NaN`. Gustavo does
   this well and lazily; the converter can be simpler because it writes once.

## Test strategy

There is **no FITS-IDI file on this machine** and none in Gustavo's fixtures —
Gustavo's tests build one by writing a synthetic `UVSet` with its own
`write_fitsidi`. XRadio needs the same trick: a minimal, test-only FITS-IDI
writer, exactly paralleling `test/synthetic.jl`, which already builds a
conforming MSv4 from nothing.

Three layers, mirroring what is already in place:

1. **Synthetic round-trip.** Build a small FITS-IDI in the test suite, convert,
   and assert `isempty(check(ps))` — the converter's output must pass the
   existing schema checker with no issues.
2. **Value round-trip.** Assert visibilities, UVW, weights, flags, antenna
   positions and frequencies survive with the values put in, including the
   axis-order reversal.
3. **Real file, by hand.** Convert one real VLBI FITS-IDI and run `xradio`'s
   `check_datatree` on the result, as was done for the writer. Needs a file —
   see CHUNK-001.

## Chunks

### CHUNK-001: preflight
- **Kind**: `preflight`
- **Description**: Establish the baseline and, critically, **obtain a real FITS-IDI file**. Confirm XRadio's suite is green at HEAD. Find a public VLBI FITS-IDI (EHT data releases, VLBA archive, or one of the author's own tracks) and record its path, size, and which optional tables it carries. Dump its HDU list and headers into this chunk's Notes — every later chunk is guessing without it. Also verify against that file: the `UU---SIN` units (seconds vs metres), the UVW sign convention, and whether `NO_BAND > 1`.
- **Status**: `todo`

### CHUNK-002: decide-partitioning
- **Kind**: `decide`
- **Description**: How does the user control the split from one FITS-IDI into N Measurement Sets? Questions: (a) is band-splitting the only mandatory axis, or does anything else in FITS-IDI force a split? (b) what is the keyword — `partition = [:source, :scan]` mirroring xradio's `partition_scheme`, or a more Julian predicate? (c) how are the resulting MSv4 nodes named, given MSv4 names must be unique within a store and xradio uses `<basename>_<zero-padded index>`? **Recommendation**: band is mandatory; `partition = ()` by default with `:source` and `:scan` as opt-in; names `<file basename>_<NN>`, matching xradio so a mixed store stays coherent.
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
- **Depends on**: CHUNK-002, CHUNK-003, CHUNK-004, CHUNK-005, CHUNK-006

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

## Open risks

- **No real FITS-IDI file yet.** CHUNK-001 is a hard gate: several conversions
  (UVW units and sign, `WEIGHTYP` in practice, which optional tables real data
  carries) cannot be settled from the spec alone, and getting them wrong is
  silent — the output would pass the schema checker and be physically wrong.
- **Scan inference and invented metadata are judgement, not translation.** They
  belong in the docs and in the store's provenance, not only in the code.
- **FLAG expansion is the fiddliest table.** Six independent range axes with
  wildcards, and no fixture generator writes one — Gustavo's tests hand-build
  the HDU. Budget for it accordingly.
