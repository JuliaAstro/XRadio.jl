# Coming from xradio

This page maps calls in the Python package
[xradio](https://github.com/casangi/xradio) (compared at version 1.1.3) to
their XRadio.jl counterparts. Both read and write the same MSv4 Zarr stores.
Where XRadio.jl has no counterpart, the table says so.

In the Python column, `ps_xdt` is a processing set `xarray.DataTree` and
`ms_xdt` one of its Measurement Sets. In the Julia column, `ps` is a
[`ProcessingSet`](@ref) and `ms` a [`MeasurementSet`](@ref).

## Names in memory

xradio keeps the on-disk names. XRadio.jl uses the lower-case form for data
variables and drops the `_xds` suffix from sub-datasets:

| On disk and in xradio | XRadio.jl |
| --- | --- |
| `VISIBILITY`, `FLAG`, `WEIGHT`, `UVW` | `:visibility`, `:flag`, `:weight`, `:uvw` |
| `antenna_xds` | `ms.antenna` |
| `field_and_source_base_xds` | `ms.field_and_source_base` |
| dimension `time`, `baseline_id`, `frequency`, `polarization` | `Ti`, `BaselineID`, `Frequency`, `Polarization`, or the store's name as a keyword |

A data group's entries, such as `datagroup(ms)[:correlated_data]`, still name
the on-disk variable, `"VISIBILITY"`. Use [`correlated`](@ref),
[`flags`](@ref), [`weights`](@ref) and [`uvw`](@ref) to reach a group's arrays
without translating the names yourself.

Scan intents are parsed into [`Intent`](@ref)s by [`scan_intents`](@ref), and
data-group names are `Symbol`s (`:base`) rather than strings.

## Opening, loading and writing

| xradio | XRadio.jl |
| --- | --- |
| `open_processing_set(path)` | `open(ProcessingSet, path)` |
| `open_processing_set(path, scan_intents=[...])` | `query(open(ProcessingSet, path); scan_intents = [...])` |
| `load_processing_set(path)` | `read(ProcessingSet, path)`, or `open(ProcessingSet, path; lazy = false)` |
| `load_processing_set(path, sel_parms={name: {dim: slice}})` | `open` lazily, then index each Measurement Set; indexing a lazy layer reads only the selected part |
| `load_processing_set(path, include_variables=[...])` / `drop_variables=[...]` | `delete!(ms, :weight)` on each Measurement Set, for example through [`transform`](@ref) |
| `load_processing_set(path, load_sub_datasets=False)` | no counterpart; a lazily opened sub-dataset is read only when accessed |
| `load_processing_set(path, data_group_name=...)` | no counterpart |
| `ps_xdt.to_zarr(path)` | `write(path, ps)`, which checks conformance first |
| `convert_msv2_to_processing_set(ms, out)` | no counterpart |
| no counterpart | [`fitsidi2msv4`](@ref), with `FITSFiles` loaded |

## The processing set

| xradio | XRadio.jl |
| --- | --- |
| `ps_xdt["name"]` | `ps[:name]` or `ps["name"]`, or `ps.name` |
| `ps_xdt.attrs` | `metadata(ps)` |
| `ps_xdt.xr_ps.get_ms_xdt()` | `only(ps)` |
| `ps_xdt.xr_ps.summary()` | no table; [`partition_info`](@ref) gives each Measurement Set's fields: `map(partition_info, values(ps))` |
| `ps_xdt.xr_ps.query(field_name="M87")` | `query(ps; field_name = "M87")`; see [Selecting data](@ref) for the value forms and what throws |
| `ps_xdt.xr_ps.query(start_frequency=slice(a, b))` | `query(ps; start_frequency = a .. b)` |
| `ps_xdt.xr_ps.query(string_exact_match=False, field_name="M8")` | `query(ps; field_name = contains("M8"))` |
| `ps_xdt.xr_ps.query(query="...")` | `filter(f, ps)` with any function of a Measurement Set |
| `ps_xdt.xr_ps.query(data_group_name=...)` | no counterpart |
| `ps_xdt.xr_ps.get_max_dims()` | no counterpart |
| `ps_xdt.xr_ps.get_freq_axis()` | no counterpart; `frequencies.(values(ps))` gives each axis |
| `ps_xdt.xr_ps.get_combined_field_and_source_xds()` | no counterpart |
| `ps_xdt.xr_ps.get_combined_field_and_source_xds_ephemeris()` | no counterpart |
| `ps_xdt.xr_ps.get_combined_antenna_xds()` | no counterpart; [`antennas`](@ref) per Measurement Set |
| `ps_xdt.xr_ps.plot_phase_centers()`, `plot_antenna_positions()` | no counterpart |
| no counterpart | [`transform`](@ref), `merge`, and [`scans`](@ref), [`fields`](@ref), [`sources`](@ref) across the set |

## A Measurement Set

| xradio | XRadio.jl |
| --- | --- |
| `ms_xdt.VISIBILITY` | `correlated(ms)` or `ms[:visibility]` |
| `ms_xdt.FLAG`, `ms_xdt.WEIGHT`, `ms_xdt.UVW` | `flags(ms)`, `weights(ms)`, `uvw(ms)` |
| `ms_xdt.attrs["data_groups"]` | [`datagroups`](@ref)`(ms)`; one group is `datagroup(ms, :base)` |
| `ms_xdt.observation_info` | `metadata(ms)[:observation_info]` |
| `ms_xdt.frequency.values` | [`frequencies`](@ref)`(ms)` |
| `ms_xdt.time.values` | [`times`](@ref)`(ms)` |
| `ms_xdt.polarization.values` | [`polarizations`](@ref)`(ms)` |
| `ms_xdt.frequency.attrs["spectral_window_name"]` | [`spectralwindow`](@ref)`(ms)` |
| `ms_xdt.antenna_xds.antenna_name.values` | [`antennas`](@ref)`(ms)` |
| `ms_xdt.antenna_xds` | `ms.antenna` |
| `ms_xdt.xr_ms.get_field_and_source_xds("base")` | [`field_and_source`](@ref)`(ms, :base)` |
| `ms_xdt.xr_ms.get_partition_info("base")` | [`partition_info`](@ref)`(ms, :base)` |
| `ms_xdt.xr_ms.add_data_group(...)` | no counterpart |
| no counterpart | [`baselines`](@ref), [`scans`](@ref), [`fields`](@ref), [`sources`](@ref), [`mounts`](@ref), [`gaincurves`](@ref) |

## Modifying a Measurement Set

| xradio | XRadio.jl |
| --- | --- |
| `ms_xdt.copy(deep=False)` | `copy(ms)`; see [Modifying data](@ref) |
| `ms_xdt.VISIBILITY.copy(data=V)` | `rebuild(ms[:visibility], V)`; a size mismatch throws when the result is assigned |
| `ms_xdt["WEIGHT"] = da` | `ms[:weight] = A`; coordinate values that differ throw instead of aligning to NaN |

## Selecting inside a Measurement Set

xradio's `xr_ms.sel` is xarray's `sel` on the correlated dataset. In
XRadio.jl, indexing takes DimensionalData selectors on the whole tree and cuts
every layer and sub-dataset that carries the dimension. [Selecting data](@ref)
gives the rules.

| xradio | XRadio.jl |
| --- | --- |
| `ms_xdt.xr_ms.sel(time=slice(t0, t1))` | `ms[time = t0 .. t1]` |
| `ms_xdt.xr_ms.sel(polarization="XX")` | `ms[polarization = At("XX")]` |
| `ms_xdt.xr_ms.sel(time=t, method="nearest")` | `ms[time = Near(t)]` |
| `ms_xdt.xr_ms.sel(baseline_id=7)` | `ms[baseline_id = At(7)]` |
| `ms_xdt.xr_ms.sel(data_group_name="corrected")` | no counterpart |
| `ms_xdt.isel(time=slice(0, 10))` | no counterpart on the tree; on one layer, `correlated(ms)[Ti(1:10)]` (1-based) |
| select a baseline by its antennas through the `baseline_antenna1_name`/`baseline_antenna2_name` coordinates | `ms[baseline(ms, ("BR", "LA"))]`, either order |
| select one antenna's baselines the same way | `ms[baselines(ms, "LA"; auto = false)]` |
| select one scan's times through the `scan_name` coordinate | `ms[scan(ms, "3")]` |
| select one field's times through the `field_name` coordinate | `ms[field(ms, "M87")]` |

A single-value selector such as `At(7)` drops the dimension, as xarray's `sel`
with a scalar does. `At([7])` keeps it.

## Checking against the standard

| xradio | XRadio.jl |
| --- | --- |
| `xradio.schema.check.check_datatree(ps_xdt)` | `XRadio.check(ps)`, a vector of [`XRadio.SchemaIssue`](@ref)s |
| `issues.expect()` | `XRadio.validate(ps)`, which throws an [`XRadio.SchemaError`](@ref) |
| `register_dataset_type(schema)` | build a specification, from scratch with [`XRadio.DatasetSpec`](@ref) or from a standard one with [`XRadio.extend`](@ref), and pass it per call: `XRadio.check(ps; schemas = [spec])`; nothing is registered globally |
