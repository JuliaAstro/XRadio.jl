# Selecting data

Selection works at two levels, and each level has its own functions:

- **Choosing Measurement Sets** from a processing set: [`query`](@ref),
  `filter`, [`transform`](@ref) and `merge`. They keep or drop whole
  Measurement Sets and return a `ProcessingSet`, so they chain. None of them
  cuts inside a Measurement Set.
  [`splitby`](@ref) cuts each Measurement Set into one per label along
  `time`, and `groupby` gathers Measurement Sets under a key.
- **Cutting inside a Measurement Set**: indexing, `ms[...]`, with
  DimensionalData selectors. [`baseline`](@ref), [`baselines`](@ref),
  [`scan`](@ref) and [`field`](@ref) turn antenna, scan and field names into
  selectors for it.

```julia
ps = open(ProcessingSet, "track.ps.zarr")

m87 = query(ps; field_name = "M87", scan_intents = "OBSERVE_TARGET")
ms = m87[1]
vis = correlated(ms)[scan(ms, "3"), baselines(ms, "LA"; auto = false)]
```

## Choosing Measurement Sets

### `query`: match the summary

`query(ps; kw...)` keeps the Measurement Sets whose summary row matches every
keyword. The columns take the reference's names and order. The row is built
from [`partition_info`](@ref) for the `:base` data group, with three fields
renamed, plus three columns of its own:

| Column | What it holds | Form |
| --- | --- | --- |
| `name` | the Measurement Set's name in the processing set | one `String` |
| `scan_intents` | the scan intents, parsed | `Vector{Intent}` |
| `polarization` | `partition_info`'s `polarization_setup`: the polarization products, in axis order | `Vector{String}` |
| `scan_name` | the scan names, sorted and unique | `Vector{String}` |
| `spw_name` | `partition_info`'s `spectral_window_name`: the `spectral_window_name` attribute of `frequency` | one `String` |
| `spw_intents` | `partition_info`'s `spectral_window_intents` | `Vector{String}` |
| `field_name` | the field-and-source dataset's field names, sorted and unique | `Vector{String}` |
| `source_name` | its source names, sorted and unique | `Vector{String}` |
| `line_name` | its line names, empty if it names none | `Vector{String}` |
| `start_frequency` | the first channel's frequency, in Hz | one `Float64` |
| `end_frequency` | the last channel's frequency, in Hz | one `Float64` |

A Measurement Set with no `frequency` axis has an empty vector in place of each
frequency, so no frequency keyword matches it. Intents a store does not state
are `UNSPECIFIED`.

**Values.** A column that holds a collection matches when any one of its
values matches. What each keyword asks for is

| Value | Matches |
| --- | --- |
| a scalar, `field_name = "M87"` | an equal value (`isequal`) |
| a vector, tuple or set, `scan_name = ["3", "4"]` | any of its elements |
| an interval, `start_frequency = 220e9 .. 230e9` | a value inside it |
| a function, `field_name = startswith("3C")` | a value for which it returns `true` |

Keywords are combined with "and". For "all of" within one column, chain two
queries.

**Intents.** A scan intent that names no subintent matches every subintent:
`scan_intents = "CALIBRATE_BANDPASS"` keeps a Measurement Set stating
`CALIBRATE_BANDPASS#ON_SOURCE`, while `"CALIBRATE_BANDPASS#ON_SOURCE"` matches
only that pair. An intent may be written as an [`Intent`](@ref), a `ScanIntents`
member or its string form. `spw_intents` are plain strings and match
exactly.

**What throws.** Every keyword is checked against each Measurement Set's row
before any is matched, so an error is not hidden because another keyword has
already ruled a Measurement Set out.

- A keyword that is not a column throws an `ArgumentError` that lists the
  columns. The check runs per Measurement Set, so on an empty processing set a
  misspelled keyword returns an empty set instead of throwing.
- A value that can never equal the column's values throws, for example a number
  asked of a column of names. An interval is checked by its left endpoint.
  Numbers compare with numbers, strings with strings, and intents with an
  `Intent`, a `ScanIntents` member or a string. A function is not checked.
- A function that returns something other than a `Bool` throws a `TypeError`.

With no keywords, `query` returns every Measurement Set. The result keeps the
Measurement Sets' names, their order and the processing set's attributes.

### `filter`: any predicate

For anything the summary does not state, `filter(f, ps)` keeps the Measurement
Sets for which `f` returns `true`. It has no rules of its own: `f` gets each
Measurement Set.

```julia
filter(ms -> "LA" in antennas(ms), ps)
filter(ms -> ("LA", "PV") in baselines(ms), ps)
filter(ms -> length(times(ms)) > 100, ps)
```

### `transform`: change each Measurement Set

`transform(f, ps)` replaces each Measurement Set with `f(ms)` under the same
name. `f` must return a `MeasurementSet`; anything else throws an
`ArgumentError` that names the Measurement Set. To read a value out of each
Measurement Set, use `map(f, values(ps))`, which returns a `Vector` in the
processing set's order.

```julia
transform(ms -> ms[scan(ms, "3")], ps)            # one scan from each
map(frequencies, values(ps))                      # Vector of frequency axes
```

Since indexing a Measurement Set by a dimension it does not carry throws,
`transform` over a mixed processing set usually follows a `query` or `filter`
that leaves only the Measurement Sets the cut applies to.

### `merge`: combine processing sets

`merge(ps, others...)` puts the Measurement Sets of all its arguments in one
processing set, in argument order, under the union of their attributes. Two
disagreements throw an `ArgumentError`:

- a Measurement Set name that two arguments share, since the two are different
  observations and keeping one would drop the other;
- an attribute that two arguments state with different values, since the result
  carries one value and `write` puts it at the root of the store.

Processing sets taken from one store always agree, which is the case this
serves:

```julia
merge(query(ps; field_name = "M87"), query(ps; field_name = "3C 273"))
```

To combine sets that disagree, construct the result with the attributes it
should have: `ProcessingSet(sets, attrs)`.

## Splitting and grouping

A Measurement Set may hold several scans or fields, and one scan is often
spread over several Measurement Sets, one per spectral window. `splitby` and
`groupby` rearrange a processing set around such labels. Both take the
processing set first and a function of one Measurement Set second.

### `splitby`: one Measurement Set per label

`splitby(ps, f)` cuts each Measurement Set along `time` by the labels `f(ms)`
gives, one per integration. The integrations sharing a label become one
Measurement Set, named `name_label`:

```julia
splitby(ps, ms -> ms[:scan_name])                    # per scan
splitby(ps, ms -> ms[:field_name])                   # per field
splitby(ps, ms -> floor.(Int, times(ms) ./ 600))     # per ten-minute bin
```

A label keeps every integration it names, so a field observed twice is one
Measurement Set with a gap in time. Labels come in the order they first occur.
`f` must give a vector as long as `times(ms)`, or an array over `time` alone
with the same times; anything else throws an `ArgumentError`, as do two pieces
whose names coincide.

Each piece is a `view`, so splitting a lazily opened processing set reads
nothing. `read(ms)` brings one Measurement Set, a piece or a whole one, into
memory.

### `groupby`: Measurement Sets under a key

`groupby(ps, f)` returns an `OrderedDict` from each key `f(ms)` to the
processing set of the Measurement Sets that give it, in the order the keys
first occur. It neither cuts nor reads, so a key describing one scan follows a
`splitby`:

```julia
groups = groupby(splitby(ps, ms -> ms[:scan_name]),
                 ms -> (only(sources(ms)), only(scans(ms))))
for ((source, scan), group) in groups
    # every spectral window of one scan
end
```

## Cutting inside a Measurement Set

### Indexing

A `MeasurementSet` is a `DimensionalData.DimTree`, so indexing takes
DimensionalData selectors (`At`, `Near`, `..`, `Where`, and so on) on the
store's dimensions. `ms[...]` copies and `view(ms, ...)` does not; both follow
the rules below, and so does a `Dataset`.

**Naming a dimension.** A dimension is named by its type or by the store's own
name as a keyword. The two select identically:

```julia
ms[BaselineID(At(7)), Ti(t0 .. t1)]
ms[baseline_id = At(7), time = t0 .. t1]
```

A dimension the standard does not define has no type of its own; the keyword is
its only name.

**A dimension no part of the tree carries** throws an `ArgumentError` that
lists the dimensions it does carry. It does not return the tree unchanged.

**Only nodes that carry a dimension are cut.** Each selector is applied to each
layer and sub-dataset that has its dimension. The rest are kept as they are, so
selecting a time range on a Measurement Set leaves the `antenna` sub-dataset
whole.

**A single-value selector drops the dimension**, as it does on an array. A
coordinate stated over that dimension alone keeps the selected value as a
zero-dimensional array. After `ms[baseline(ms, ("BR", "LA"))]` the Measurement
Set has no `baseline_id` dimension, and `ms[:baseline_antenna1_name][]` is
`"BR"`. A one-element vector keeps the dimension: `BaselineID(At([7]))`.

On a single layer, indexing is DimensionalData's own and none of these rules
apply:

```julia
correlated(ms)[Polarization(At("XX")), Frequency(1:64)]
```

### Selectors from names

MSv4 indexes baselines by an integer `baseline_id` and states scans and fields
as a name per integration. These functions return a selector for the positions
that a name covers. They do not return data, so the selector works on the whole
Measurement Set or on one layer, and combines with other selectors:

| Call | Selector | Selects |
| --- | --- | --- |
| `baseline(ms, ("BR", "LA"))` | `BaselineID(At(id))` | one baseline, either antenna order |
| `baseline(ms, "BR-LA")` | `BaselineID(At(id))` | the same pair, as a dash label |
| `baseline(ms, 7)` | `BaselineID(At(7))` | a `baseline_id`, checked to exist |
| `baselines(ms, "LA"; auto)` | `BaselineID(At(ids))` | every baseline with `"LA"`; `auto = false` drops its autocorrelation |
| `scan(ms, "3")` | `Ti(At(times))` | the integrations of scan `"3"` |
| `field(ms, "M87")` | `Ti(At(times))` | the integrations on field `"M87"` |

A single `baseline` is a single-value selector and drops the `baseline_id`
dimension. `baselines`, `scan` and `field` select a vector and keep theirs.
`scan` and `field` select by time value, not position, so a selector taken
before a cut along another dimension, such as `baseline_id`, still applies
after it.

A name the Measurement Set does not have throws an `ArgumentError` that lists
the names it does have. A dash label must split into exactly two names, so an
antenna whose name contains a dash has to be given as a pair.

Without an argument, [`baselines`](@ref), [`scans`](@ref) and [`fields`](@ref)
list the names these accept.

### Cached lookups

Resolving a name requires a pass over the coordinates, so each Measurement Set
builds a lookup the first time it is asked and keeps it:
[`XRadio.baseline_lookup`](@ref) for baselines, and
[`XRadio.scan_lookup`](@ref) and [`XRadio.field_lookup`](@ref) for scans and
fields.

**Invalidation.** Assigning, deleting or popping a layer (`ms[:key] = A`,
`delete!`, `pop!`) drops all three lookups, and the next query rebuilds them.
Changing a coordinate's values in place, through the array `ms[:key]` returns,
does not, and the lookups then describe the old values.

**Holding a lookup.** A Measurement Set stores its layers with abstract types,
so fetching the cached lookup costs a dynamic dispatch on each call. A loop can
fetch it once and pass it to the lookup forms of the selector functions:

```julia
bl = XRadio.baseline_lookup(ms)
selectors = [baseline(bl, p) for p in pairs_wanted]

sl = XRadio.scan_lookup(ms)
t = collect(DimensionalData.val(dims(ms, Ti)))
per_scan = [correlated(ms)[scan(sl, t, s)] for s in scans(ms)]
```

## Against the reference

The reference is the Python package
[xradio](https://github.com/casangi/xradio), compared here at version 1.1.3:
`ProcessingSetXdt.query` and `MeasurementSetXdt.sel`.

| | xradio | XRadio.jl |
| --- | --- | --- |
| Choose Measurement Sets | `ps_xdt.xr_ps.query(**kw)` | `query(ps; kw...)` |
| Summary columns | `name`, `scan_intents`, `shape`, `execution_block_UID`, `polarization`, `scan_name`, `spw_name`, `spw_intents`, `field_name`, `source_name`, `line_name`, `field_coords`, `session_reference_UID`, `scheduling_block_UID`, `project_UID`, `start_frequency`, `end_frequency` | the same names, less `shape`, `execution_block_UID`, `field_coords`, `session_reference_UID`, `scheduling_block_UID` and `project_UID` |
| Any of several values | a `list`; a `tuple` is one value | a vector, tuple or set |
| A range | `slice(a, b)`, both ends included | an interval, `a .. b` |
| Substring match | `string_exact_match=False`, applied to every keyword | a function per keyword, `contains("CAL")` |
| Free-form condition | `query="..."`, a pandas query string on the summary | `filter(f, ps)`, any function of a Measurement Set |
| Scan intents | compared as whole strings | an intent without a subintent matches every subintent |
| Unknown keyword | `KeyError` from pandas | `ArgumentError` listing the columns |
| Value of the wrong type | matches nothing | `ArgumentError` |
| `data_group_name` | not matched; narrows each kept Measurement Set to that data group | not accepted |
| Change each Measurement Set | no method | `transform(f, ps)` |
| Combine processing sets | no method | `merge(ps, others...)` |
| Cut inside a Measurement Set | `ms_xdt.xr_ms.sel(**indexers)`, xarray's `sel` on the correlated dataset | `ms[dim = selector, ...]` on the whole tree |
| Unknown dimension in a cut | error from xarray | `ArgumentError` listing the dimensions |
| Data group in a cut | `sel(data_group_name=...)` drops the other groups' variables | no counterpart |
| Name-based selectors | none; select on the coordinates | `baseline`, `baselines`, `scan`, `field` |

The reference's summary also holds the correlated data's `shape`, a
`field_coords` description and the `observation_info` UIDs; XRadio's row does
not. `shape`, a bare
tuple, is the only column a tuple can match in the reference; no column of
XRadio.jl's row holds a tuple, so a tuple here is always any of its elements.

## Docstrings

```@docs; canonical=false
query
partition_info
Base.filter(::Any, ::ProcessingSet)
transform
Base.merge(::ProcessingSet, ::ProcessingSet...)
baseline
baselines
scan
field
XRadio.baseline_lookup
XRadio.scan_lookup
```
