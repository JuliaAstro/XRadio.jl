# Building data

[`XRadio.measurement_set`](@ref) and [`XRadio.subdataset`](@ref) build MSv4
datasets from plain arrays. The schema supplies the dimensions, the axis order
and the attributes the standard fixes; the caller supplies values and the
attributes the standard requires but cannot know, such as a spectral window's
name or a direction's frame. Arrays are in Julia order, fastest-varying
dimension first, which is the reverse of the order the standard states.

```@docs; canonical=false
XRadio.measurement_set
```
