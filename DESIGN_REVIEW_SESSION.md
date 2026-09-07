# Session Handoff — 2026-09-06

## Plan
DESIGN_REVIEW_PLAN.md — XRadio 0.1.0 (FITS-IDI → MSv4 converter)

## What was just completed
CHUNK-005: implement-fitsidi-reader. One file, `ext/XRadioFITSFilesExt/fitsidi.jl`
(346 lines): access primitives for cards and columns, the arithmetic FITS-IDI
states rather than stores, and one bulk reader for `UV_DATA`'s per-row columns.
`src/` gains a single declaration, `fitsidi2msv4`, and no FITS-IDI vocabulary.

There is **no intermediate data model**. An earlier pass in this session built
eight structs mirroring the tables and it was rejected as a representation for a
FITS-IDI library rather than a one-way converter — sixteen of its fields had no
MSv4 target. The standing preference is simple, readable code; structs are
welcome where they make a call site clearer and should be introduced then,
driven by need rather than anticipated.

## Key decisions / shim choices
- **Depend on the dev FITSFiles**, per the author's instruction mid-session.
  Registered `0.3.2` and `~/.julia/dev/FITSFiles` (branch `ptiede-uvfixes`,
  four commits ahead of upstream `main`) report the same version and are not
  the same code. The registered one has no `LazyStructuredData`: `hdu.data`
  materializes the whole HDU on **every** access, so touching `f[9].data` on
  the Q-band file reads all 24.73 GB. `test/Manifest.toml` now `dev`s the
  branch. Full note under the plan's `Decisions`.
- **Nothing FITS-IDI goes in `src/`.** XRadio implements MSv4; FITS-IDI
  vocabulary in its namespace would widen the subject. The single `src/`
  declaration is `fitsidi2msv4`, which CHUNK-007 will give a method. Tests
  reach the reader through `Base.get_extension(XRadio, :XRadioFITSFilesExt)`.
  This also means CHUNK-007 through CHUNK-009 are extension code.
- **`fitsidi2msv4(path, dest)` writes to disk and returns `dest`**, rather than
  returning a `ProcessingSet`. Partitions are written as they are built, so
  peak memory is one partition (270 MB here) whatever the input size. The
  declaration says so; CHUNK-007 has to honour it.
- **Polarization follows Gustavo, not MSv4.** `pol_axis` decodes the `STOKES`
  axis to generic products `PP`/`QQ`/`PQ`/`QP` named by receptor index; the
  basis stays per antenna, in `POLTYA`/`POLTYB`. MSv4's array-wide
  `"RR"`/`"RL"` coordinate presumes a homogeneous array, and a mixed
  circular/linear baseline measures `RX`, `LY`, `RY`, `LX`, which no array-wide
  label names. `pol_labels(products, basis1, basis2)` resolves it per baseline.
  Nothing is pre-permuted: `pol_axis` returns the file's order plus the
  permutation, and each consumer applies it once.
- **`scan_ranges` is a free function over `IDIRows`, not `scans`.** `scans` and
  `times` are already exported accessors returning MSv4 scan names and the time
  coordinate; adding methods with unrelated return types would have made those
  generics mean two things. Hence `scan_ranges` and `julian_dates`.
- **The row reader bypasses FITSFiles' column access deliberately** and so
  applies neither `TSCAL`/`TZERO` nor multi-element layout. Both would be
  silent, so `_scalar` rejects any column carrying either.

## API changes made alongside, at the author's direction
Not CHUNK-005 work, but they touch the package's export surface and CHUNK-011
inherits them:
- `check`, `validate`, `SchemaIssue`, `SchemaError` are no longer exported —
  public, but conformance checking is not what most work with a store involves.
  `test/runtests.jl` imports them explicitly; the README shows the import.
- DimensionalData is re-exported via `Reexport` (new dependency, compat `"1"`).
  A `MeasurementSet` is a `DimTree`, so `At`, `Near`, `dims` and `lookup` are
  needed for any use at all. XRadio's 49 exports and DimensionalData's 70 have
  zero name overlap, so nothing shadows.
- `UvwLabel` is now `UVWLabel`. The MSv4 dimension name `uvw_label` is
  untouched; the reference JSON under `test/reference/` mentions
  `UvwLabelArray`, which is xradio's own Python schema name and must stay.

## State of the codebase
- Files modified: `Project.toml` (weakdep, extension, `Reexport`, compat),
  `test/Project.toml`, `test/runtests.jl`, `test/synthetic.jl`, `README.md`,
  `src/XRadio.jl`, `src/dimensions.jl`, `DESIGN_REVIEW_PLAN.md`,
  `DESIGN_REVIEW_SESSION.md`
- Files added: `src/fitsidi.jl` (the `fitsidi2msv4` declaration and nothing
  else), `ext/XRadioFITSFilesExt/XRadioFITSFilesExt.jl`,
  `ext/XRadioFITSFilesExt/fitsidi.jl`, `test/test_fitsidi.jl`
- Test suite: pass — 1588, up from the 1532 baseline
- Ambiguity count: 0 (baseline 0)
- Staged but uncommitted: no — nothing is staged; the working tree carries
  everything for review

## Cluster status
No clusters in this plan. 6 of 15 chunks complete.

## Next chunk
CHUNK-006: implement-test-fitsidi-writer — a minimal FITS-IDI writer under
`test/`, enough to build fixtures with known values. Model on Gustavo's
`ext/GustavoFITSFilesExt/fitsidi_write.jl` (478 lines), keeping only what a
fixture needs. It is now **load-bearing**: with no intermediate model there is
nothing else to assert on, so it is the only way to test CHUNK-007 without the
24 GB file. Test-only — under `test/`, never `src/`, never public.

## Watch out for
- **Never touch `hdu.data` on `UV_DATA`.** Use `getfield(hdu, :data)` — the
  extension's `_descriptor` — which returns the lazy descriptor without reading
  the data block. This session lost several minutes to a `f[9].data` that began
  streaming the whole 24.73 GB; the only way to stop it was restarting the REPL,
  since Kaimon's cancellation is cooperative.
- **The Kaimon MCP server does work**, contrary to the last two handoffs. The
  no-arg `start_session` reported "no allowed projects" while
  `~/.config/kaimon/projects.json` already listed XRadio; passing
  `project_path` explicitly started a session immediately.
- **The REPL runs with `test/` activated**, because FITSFiles is a weak
  dependency and cannot be loaded from the main project environment.
- **Q-band file ground truth**, for writing further assertions: 10 stations
  (BR FD HN KP LA MK NL OV PT SC, `NOSTA` 1:10), 4 bands x 256 channels x 4
  pols (RR LL RL LR), `REF_FREQ` 42.912 GHz, bands offset 0/128/256/384 MHz,
  channel width 0.5 MHz, 3 sources (3C273, M87, 3C279), 7245 FLAG rows, 752160
  `UV_DATA` rows, 55 baselines including autocorrelations, 0.5 s integrations,
  36 scans.
- **CHUNK-007 still owes a decision on the MSv4 `polarization` coordinate.**
  The store must carry one array-wide list to conform, and for a heterogeneous
  array any value is wrong for some baseline. Both validation files are
  homogeneously circular, so the test suite cannot catch a silent choice.
- **The UVW sign convention is still unverified** — unchanged from the last
  session, and still CHUNK-007's job. `IDIRows.uvw` holds the file's values
  verbatim, in **seconds**, so the seconds→metres conversion has not happened
  yet either.
- `read_fitsidi_rows` takes ~23 s on the Q-band file. It allocates one
  `n * 48`-byte buffer (36 MB here), not a per-row read.
