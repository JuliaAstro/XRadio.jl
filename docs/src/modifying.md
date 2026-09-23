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

```@docs; canonical=false
Base.copy(::Union{MeasurementSet, Dataset})
```
