"""
The MSv4 dimension vocabulary.

MSv4 names dimensions with lowercase `snake_case` strings; `DimensionalData`
names them with singleton `Dimension` types. Every dimension the schema uses
gets a type here, and [`MSV4_DIMS`](@ref) is the bijection between the two
spellings. `time` maps to `DimensionalData`'s built-in `Ti` so that its time
selectors work on MSv4 data.
"""

@dim BaselineId "Baseline identifier"
@dim Frequency "Frequency (Hz)"
@dim Polarization "Polarization product"
@dim UvwLabel "UVW component"

@dim AntennaName "Antenna name"
@dim StationName "Station name"
@dim FieldName "Field name"
@dim ReceptorLabel "Receptor (feed) label"
@dim ElementId "Phased-array element identifier"

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

See also [`dimtype`](@ref) and [`msv4name`](@ref).
"""
const MSV4_DIMS = (
    time = Ti,
    baseline_id = BaselineId,
    frequency = Frequency,
    polarization = Polarization,
    uvw_label = UvwLabel,
    antenna_name = AntennaName,
    station_name = StationName,
    field_name = FieldName,
    receptor_label = ReceptorLabel,
    element_id = ElementId,
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
    dimtype(name) -> UnionAll

The `Dimension` type for MSv4 dimension `name`, given as a `Symbol` or
`AbstractString`. Throws for a name outside the MSv4 vocabulary.

```jldoctest
julia> XRadio.dimtype("baseline_id")
BaselineId
```
"""
dimtype(name::Symbol) = get(MSV4_DIMS, name) do
    throw(ArgumentError("`$name` is not an MSv4 dimension. Known: $(join(sort(collect(string.(keys(MSV4_DIMS)))), ", "))"))
end
dimtype(name::AbstractString) = dimtype(Symbol(name))

"""
    msv4name(D) -> Symbol

The MSv4 dimension name for a `Dimension` type or instance — the inverse of
[`dimtype`](@ref).

```jldoctest
julia> XRadio.msv4name(BaselineId)
:baseline_id
```
"""
msv4name(::Type{D}) where {D <: Dimension} = _dimname(Base.typename(D).wrapper)
msv4name(d::Dimension) = msv4name(typeof(d))

function _dimname(W::UnionAll)
    return get(_DIM_TO_NAME, W) do
        throw(ArgumentError("`$W` is not an MSv4 dimension"))
    end
end
