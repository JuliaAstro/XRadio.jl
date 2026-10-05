# Writing stores

`write(path, ps)` writes a processing set as a new MSv4 Zarr store, checking
it against the standard first, and `write(path, name => ms)` writes one
Measurement Set. The `mode` keyword says what happens to a store already at
`path`:

| `mode` | `path` holds a store | `path` holds nothing |
| --- | --- | --- |
| `"w-"` (default) | error | create |
| `"w"` | replace it | create |
| `"a"` | write into it | create |

These are xarray's `to_zarr` modes. To change values inside a stored array,
as xarray's `"r+"` does, open the store with `mode = "r+"` and assign into its
layers:

```julia
ps = open(ProcessingSet, path; mode = "r+")
ps[:track_0][:flag][Ti(1:10)] .= true
```

This writes only the chunks the assignment covers. Shapes, coordinates and
dimensions cannot change, and in a Zarr v2 store string-valued layers refuse
assignment. The write is not atomic, and other readers see each chunk as soon
as it is written.

## Chunk shapes

Zarr compresses and reads an array a chunk at a time, so the chunk shape
decides what reading part of a store costs. Each dataset's layers share one
set of chunk sizes, so the chunks of the visibilities, flags and weights cover
the same samples. Index coordinates are always one chunk.

By default the chunks suit reading a scan at a time: every baseline, channel
and polarization of a run of time steps. Time is cut so that a chunk of the
dataset's largest layer holds about `chunk_bytes` (16 MiB) uncompressed. If one
time step alone exceeds that, frequency is cut too, and then baseline.
Polarization is not cut: it has at most four entries, it varies fastest on
disk, and cross-hand work reads all of it.

`chunks` fixes the sizes of the dimensions it names, by their MSv4 names, and
the rest are fitted to `chunk_bytes` as above:

```julia
write(path, ps)                                             # a scan at a time
write(path, ps; chunks = (; baseline_id = 1))               # a baseline at a time
write(path, ps; chunks = (; frequency = 16))                # a narrow band at a time
write(path, ps; chunks = (; polarization = 1))              # parallel hands alone
write(path, ps; chunks = (time = 64, frequency = 128), chunk_bytes = typemax(Int))
```

A store has one chunk shape, and a shape that makes one kind of read cheap
makes another dear: with one baseline per chunk, reading a scan touches a
chunk per baseline. Changing the shape means rewriting the arrays. A name in
`chunks` that is not a dimension of the data is an error.

## Writing into an existing store

```julia
ps = open(ProcessingSet, path)
ms = copy(ps[:track_0])
ms[:visibility_corrected] = rebuild(ms[:visibility], corrected)
add_datagroup!(ms, :corrected; correlated_data = "VISIBILITY_CORRECTED")
write(path, :track_0 => ms; mode = "a")
```

This writes what the store's `track_0` does not hold, here the new layer, and
replaces the attributes of the Measurement Set and its sub-datasets, which
carry the new data group. The stored visibilities are not rewritten. The
store's other Measurement Sets and its own attributes are left as they are;
writing a processing set with `mode = "a"` writes each of its Measurement Sets
this way and replaces the store's attributes with the set's. The store keeps
its Zarr version. This is xarray's `to_zarr(mode = "a")`, with these
differences:

- A layer the store already holds is left alone only when `ms` holds it as
  read lazily from this store. A layer holding anything else is an error unless
  `replace = true`, which rewrites it; xarray overwrites it. String layers,
  which are read into memory even lazily, are compared by value.
- A coordinate that differs from the stored one is an error. One `ms` adds is
  written.
- `replace = true` is accepted only with `mode = "a"`.

Layers the store holds and `ms` does not stay in the store, and the attributes
of layers left alone are not rewritten. A Measurement Set the store does not
hold is written whole.

```@docs; canonical=false
Base.write(::AbstractString, ::Pair{Symbol, MeasurementSet})
```

## Writing a template, then filling it

`write(path, ps; data = false)` writes a store shaped like `ps` without its
data: coordinates, attributes and sub-datasets are written as usual, and each
Measurement Set's data variables are created empty. Results are then written
into it region by region, through the layers of a store opened with
`mode = "r+"`, without holding the whole result in memory:

```julia
write(out, ps; data = false)
ms = open(ProcessingSet, out; mode = "r+")[:track_0]
ms[:visibility][Ti(1:10)] .= vis      # the first ten integrations
ms[:flag][Ti(1:10)] .= flags
```

An unwritten sample reads as flagged: `FLAG` is `true` there, a real-valued
layer `NaN`, and a complex one zero. A half-filled store therefore conforms and
reads correctly, and the empty layers take no space on disk beyond their
metadata. `data = false` also works with `mode = "a"`, creating the layers a
Measurement Set adds empty, but cannot be combined with `replace = true`.
