"""
The MSv4 dimension vocabulary.

MSv4 names dimensions with lowercase `snake_case` strings; `DimensionalData`
names them with singleton `Dimension` types. Every dimension the schema uses
gets a type here, and [`MSV4_DIMS`](@ref) is the bijection between the two
spellings. `time` maps to `DimensionalData`'s built-in `Ti` so that its time
selectors work on MSv4 data.

The vocabulary is open at the edges. A store may carry data the standard does
not describe, and a dataset added to one brings its own dimensions with it, so
a name outside `MSV4_DIMS` is labelled with `DimensionalData`'s generic
`Dim{name}` rather than rejected. That type carries the name itself, which is
what lets such a dimension survive a read and a write without anything having
to be told about it.
"""

@dim BaselineID "Baseline identifier"
@dim Frequency "Frequency (Hz)"
@dim Polarization "Polarization product"
@dim UVWLabel "UVW component"

@dim AntennaName "Antenna name"
@dim StationName "Station name"
@dim FieldName "Field name"
@dim ReceptorLabel "Receptor (feed) label"
@dim ElementID "Phased-array element identifier"

"""
    AntennaName

The `antenna_name` dimension: unique name for each antenna(_station).
"""
AntennaName

"""
    StationName

The `station_name` dimension of the weather sub-dataset: station name.
"""
StationName

@dim CartesianPosLabel "Cartesian position component"
@dim CartesianPosLabelLocal "Cartesian position component, antenna-local frame"
@dim EllipsoidDirLabel "Ellipsoidal direction component"
@dim EllipsoidDisLabel "Ellipsoidal distance component"
@dim SkyDirLabel "Sky direction component"
@dim SkyDisLabel "Sky distance component"
@dim LocalSkyDirLabel "Local (azimuth/elevation) sky direction component"

@dim LineLabel "Spectral line label"
@dim ToneLabel "Phase-calibration tone label"
@dim PolyTerm "Polynomial term"
@dim NPolynomial "Pointing polynomial term"

@dim TimeEphemeris "Ephemeris sample time"
@dim TimePointing "Pointing sample time"
@dim TimeSystemCal "System-calibration sample time"
@dim TimeWeather "Weather sample time"
@dim TimePhaseCal "Phase-calibration sample time"
@dim FrequencySystemCal "System-calibration frequency"

"""
    MSV4_DIMS

Maps each MSv4 dimension name to the `Dimension` type representing it. This is
the complete set of dimensions the MSv4 schema defines; a store naming anything
else does not conform.

This is the complete set of dimensions the MSv4 schema defines. It is not the
complete set a store may name: anything else is labelled `Dim{name}`.

See also [`dimtype`](@ref) and [`msv4name`](@ref).
"""
const MSV4_DIMS = (
    time = Ti,
    baseline_id = BaselineID,
    frequency = Frequency,
    polarization = Polarization,
    uvw_label = UVWLabel,
    antenna_name = AntennaName,
    station_name = StationName,
    field_name = FieldName,
    receptor_label = ReceptorLabel,
    element_id = ElementID,
    cartesian_pos_label = CartesianPosLabel,
    cartesian_pos_label_local = CartesianPosLabelLocal,
    ellipsoid_dir_label = EllipsoidDirLabel,
    ellipsoid_dis_label = EllipsoidDisLabel,
    sky_dir_label = SkyDirLabel,
    sky_dis_label = SkyDisLabel,
    local_sky_dir_label = LocalSkyDirLabel,
    line_label = LineLabel,
    tone_label = ToneLabel,
    poly_term = PolyTerm,
    n_polynomial = NPolynomial,
    time_ephemeris = TimeEphemeris,
    time_pointing = TimePointing,
    time_system_cal = TimeSystemCal,
    time_weather = TimeWeather,
    time_phase_cal = TimePhaseCal,
    frequency_system_cal = FrequencySystemCal,
)

# Reverse lookup, built once from MSV4_DIMS so the two spellings cannot drift.
const _DIM_TO_NAME = Dict{UnionAll, Symbol}(
    v => k for (k, v) in pairs(MSV4_DIMS)
)

"""
    dimtype(name) -> Type

The `Dimension` type for MSv4 dimension `name`, given as a `Symbol` or
`AbstractString`.

A name the standard does not define gets `Dim{name}`, which labels the
dimension and states its own name back. To work with such a dimension under a
type of your own, `set` it after reading and set it back before writing:

```julia
DimensionalData.@dim TonePhase "Phase-cal tone phase"
mine = set(array, Dim{:tone_phase} => TonePhase)
theirs = set(mine, TonePhase => Dim{:tone_phase})
```

```jldoctest
julia> XRadio.dimtype("baseline_id")
BaselineID

julia> XRadio.dimtype(:tone_phase)
Dim{:tone_phase}
```
"""
dimtype(name::Symbol) = get(MSV4_DIMS, name, Dim{name})
dimtype(name::AbstractString) = dimtype(Symbol(name))

"""
    msv4name(D) -> Symbol

The MSv4 dimension name for a `Dimension` type or instance — the inverse of
[`dimtype`](@ref).

Throws for a dimension type of your own, which carries no MSv4 name: writing it
would put a name of the writer's invention into the store. `Dim{name}` is the
spelling a dimension outside the standard has in a store, so set one back to it
before writing.

```jldoctest
julia> XRadio.msv4name(BaselineID)
:baseline_id

julia> XRadio.msv4name(Dim{:tone_phase})
:tone_phase
```
"""
msv4name(::Type{D}) where {D <: Dimension} = _dimname(Base.typename(D).wrapper)
msv4name(::Type{<:Dim{S}}) where {S} = S::Symbol
msv4name(d::Dimension) = msv4name(typeof(d))

function _dimname(W::UnionAll)
    return get(_DIM_TO_NAME, W) do
        throw(
            ArgumentError(
                "`$W` is not an MSv4 dimension and carries no MSv4 name. A                  dimension the standard does not define is spelled                  `Dim{name}` in a store; `set` this one back to that before                  writing."
            )
        )
    end
end
