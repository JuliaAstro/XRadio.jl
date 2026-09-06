# Session Handoff — 2026-09-06

## Plan
DESIGN_REVIEW_PLAN.md — XRadio 0.1.0 (FITS-IDI → MSv4 converter)

## What was just completed
CHUNK-002, CHUNK-003, CHUNK-004, CHUNK-014 — all four ready `decide` chunks,
batched per the skill's exception. Every answer is recorded in full under the
plan's `Decisions` section; the chunk `Notes` carry only pointers. Three of the
four went with the plan's own recommendation. CHUNK-004 did not: measuring
astropy's behaviour reversed it.

## Key decisions / shim choices
- **Partitioning**: spectral window + scan, keyword `partition = (:scan,)`,
  node names `<basename>_<NN>`. 144 partitions of 270 MB at 100% fill.
- **Invented metadata**: `["UNSPECIFIED"]` for both intents,
  `"1858-11-17T00:00:00.000"` for `release_date`; all three overridable, warn on
  `release_date` only. `observer` comes from the primary header.
- **Time**: `format = "unix"`, reversing the plan's `mjd` recommendation.
  MSv4 pins `units = "s"`, and `Time(mjd_days * u.s, format="mjd")` silently
  yields 1858 rather than the intended date. Writing MJD then means
  seconds-since-MJD-0, which removes the epoch-shift advantage and makes the
  precision 4x *worse* than unix (9.5e-7 s vs 2.4e-7 s ulp).
- **Frames**: `EQUINOX "J2000"` → `fk5`, literally, erroring on anything else.
  `UVW.frame` inherits the phase center's.
- **Streaming**: deferred. Scan partitioning made it an optimization, not a
  blocker, so CHUNK-007 no longer depends on CHUNK-014.

## State of the codebase
- Files modified: `DESIGN_REVIEW_PLAN.md`, `DESIGN_REVIEW_SESSION.md` — no
  package source touched.
- Test suite: n/a (no source change; baseline from CHUNK-001 is 1532 passing)
- Ambiguity count: n/a (baseline 0)
- Staged but uncommitted: no

## Cluster status
No clusters defined in this plan. 5 of 15 chunks complete; the four `decide`
chunks are now all closed and the implementation path is unblocked.

## Next chunk
CHUNK-005: implement-fitsidi-reader — an `XRadioFITSFilesExt` extension on
`FITSFiles` parsing the FITS-IDI tables into plain Julia structures. Produces
no MSv4 types, so it is testable against the real file on its own.

## Watch out for
- **The Julia MCP server (`kaimon`) has now failed to connect for two sessions
  running.** Everything here was measured from the shell instead. CHUNK-005 is
  the first chunk that really wants Revise, so check the server is up before
  starting, and expect slow iteration if it is not.
- **FITS strings are padded, and every symptom of ignoring that is silent.**
  Card values are space-padded (`EXTNAME` is `"SOURCE  "`); table string columns
  are NUL-padded (`"3C273\0\0\0..."`). Details in CHUNK-005's Notes.
- **`hdu.data` is a `LazyStructuredData`**: `Symbol` keys, supports `keys`, does
  *not* support `haskey`.
- **`RAAPP`/`DECAPP` are verbatim copies of `RAEPO`/`DECEPO`** in this file —
  DiFX never computed apparent places. Reading them as apparent would be wrong
  and undetectable.
- The UVW **sign convention remains the one unverified question** from CHUNK-001.
  Headers cannot settle it; it needs actual `BASELINE` values against `STABXYZ`
  positions, in CHUNK-007.
- A scratch environment with the FITSFiles dev checkout is at
  `<scratchpad>/Project.toml`; the xradio source and a venv with astropy 8.0.1
  are under the previous session's scratchpad, path recorded in this session's
  transcript.
