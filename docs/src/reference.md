# Reference

```@autodocs
Modules = [XRadio]
Private = false
```

## Methods of Base functions

```@docs
Base.open(::Type{ProcessingSet}, ::Any)
Base.read(::Type{ProcessingSet}, ::Any)
Base.write(::AbstractString, ::ProcessingSet)
Base.isvalid(::Union{Dataset, MeasurementSet, ProcessingSet})
Base.copy(::Union{MeasurementSet, Dataset})
```

## Intents

```@autodocs
Modules = [XRadio.ScanIntents, XRadio.SubscanIntents]
```

## Internals named above

```@docs
XRadio.relaxed_eltype
XRadio.GAIN_VARIABLES
XRadio.GAIN_QUANTITIES
```
