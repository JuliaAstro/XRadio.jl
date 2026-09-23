# Modifying data

A modified Measurement Set is built from a copy, in the same three steps as
xarray's `ds.copy(deep=False)`, `da.copy(data=...)` and `ds["NAME"] = ...`:

```julia
ms2 = copy(ms)                                     # shares the arrays
ms2[:weight] = rebuild(ms[:weight], new_weights)   # same labels, new data
```

Over a processing set, the same steps go inside [`transform`](@ref):

```julia
transform(ps) do ms
    out = copy(ms)
    out[:weight] = rebuild(ms[:weight], reweight(ms))
    return out
end
```

## `copy`

`copy(ms)` shares every array with `ms` and nothing else. Each node, layer and
coordinate of the copy has its own attribute dictionary, and each sub-dataset
is a copy too, so adding, replacing or deleting a layer or attribute of the
copy leaves `ms` unchanged. Values nested inside an attribute, such as one
entry of `data_groups`, are shared, as in xarray. A lazily opened Measurement
Set stays lazy: the copy reads from the same store.

`copy` of a sub-dataset on its own, `copy(ms.antenna)`, keeps the dimensions
it reads from its Measurement Set.

## Sub-datasets

A sub-dataset is linked to its Measurement Set, as a child node is in xarray's
`DataTree`: it holds only its own dimensions, such as `antenna_name`, and reads
`time`, `frequency` and every other dimension of the Measurement Set from it.
A system calibration sampled on the correlated data's own time axis therefore
carries no time axis of its own, and `write` stores `time` once.

Assigning a sub-dataset, `ms2.system_calibration = xds` or
`ms2[:system_calibration] = xds`, links a copy of it, so a sub-dataset taken
from another Measurement Set stays linked to that one. A dimension both
state must agree in length and coordinate values, as xarray requires of a node
and its parents, or the assignment throws `DimensionMismatch`. A dimension the
sub-dataset only numbers takes its parent's coordinates.

A copy taken from another Measurement Set keeps the dimensions it read there,
so it can be assigned only to a Measurement Set with the same coordinates.
This is also what xarray does.

Opening a store links its sub-datasets the same way, and a store whose
sub-dataset disagrees with its Measurement Set throws, as `xr.open_datatree`
does. A sub-dataset placed in `branches(ms)` directly is not linked;
[`XRadio.check`](@ref) reports one that disagrees with its parent.

## Replacing a layer

`rebuild(A, data)` is DimensionalData's: it wraps `data` in `A`'s dimensions
and attributes. Assigning the result, `ms2[:key] = A`, checks it against the
tree and throws on:

- data whose size is not the lengths of its dimensions (`DimensionMismatch`);
  `rebuild` itself does not check this;
- a dimension whose length differs from the tree's (`DimensionMismatch`);
- a dimension whose coordinate values differ from the tree's
  (`ArgumentError`, naming the dimension).

The last differs from xarray, which aligns the new array to the existing
coordinates and fills the positions that do not match with NaN. Here the
coordinates of a Measurement Set change only by indexing it.

A replacement with other dimensions than the layer it replaces, such as a 2-D
array for the 4-D visibilities, is accepted, as in xarray. [`XRadio.check`](@ref)
reports it, and `write` refuses it.

Assigning or deleting a layer drops the Measurement Set's cached baseline, scan
and field lookups; see [Cached lookups](@ref).

## Data groups

A new version of the correlated data, flags or weights is a new layer, and a
data group is what says which of them belong together. [`add_datagroup!`](@ref)
adds one, taking every entry it is not given from another group, `base` unless
`shared_with` says otherwise:

```julia
ms2 = copy(ms)
ms2[:visibility_corrected] = rebuild(ms[:visibility], corrected)
add_datagroup!(ms2, :corrected; correlated_data = "VISIBILITY_CORRECTED")
```

Unlike xradio's `add_data_group`, it throws on a group name already in use and
on a name that is not a layer or sub-dataset of the Measurement Set, and a group
given no `date` or `description` gets the current time and an empty
description rather than those of the group it shares with. It replaces the
`data_groups` attribute rather than editing it, so `ms` keeps its groups.

```@docs; canonical=false
Base.copy(::Union{MeasurementSet, Dataset})
```
