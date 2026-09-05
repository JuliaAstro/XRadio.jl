"""
Measure attribute bundles.

Every MSv4 array that carries a physical quantity declares what kind of measure
it is in its `type` attribute, and the kind fixes which further attributes must
accompany it — a `sky_coord` needs a `frame`, a `time` needs a `scale` and a
`format`, and so on. The helpers here build those bundles so a dataset table
stays one readable line per array.

The permitted values are the standard's, transcribed from the MSv4 schema.
"""

const TIME_SCALES = ["tai", "tcb", "tcg", "tdb", "tt", "ut1", "utc"]
const TIME_FORMATS = ["unix", "mjd", "cxcsec", "gps"]

const SKY_FRAMES = [
    "altaz", "barycentricmeanecliptic", "barycentrictrueecliptic", "cirs",
    "custombarycentricecliptic", "fk4", "fk4noterms", "fk5", "galactic",
    "galacticlsr", "galactocentric", "gcrs", "geocentricmeanecliptic",
    "geocentrictrueecliptic", "hadec", "hcrs", "heliocentriceclipticiau76",
    "heliocentricmeanecliptic", "heliocentrictrueecliptic", "icrs", "itrs",
    "lsr", "lsrd", "lsrk", "precessedgeocentric", "supergalactic", "teme", "tete",
]

const UVW_FRAMES = ["icrs", "fk5", "APP"]

const SPECTRAL_OBSERVERS = [
    "TOPO", "BARY", "REST", "gcrs", "hcrs", "icrs", "lsr", "lsrd", "lsrk",
]

const DOPPLER_TYPES = [
    "radio", "z", "ratio", "beta", "gamma", "optical", "true", "relativistic",
]

const COORDINATE_SYSTEMS = [
    "geocentric", "geodetic", "orbital", "planetcentric", "planetodetic", "topocentric",
]

const ELLIPSOIDS = ["GRS80", "WGS72", "WGS84"]

const QUANTITY_UNITS = ["s", "Hz", "m", "m/s", "rad", "K", "K/Jy", "Pa", "/m^2"]

_type_attr(v::AbstractString) = AttrSpec(:type, :str; values = [v], default = v)

"""
    quantity(units) -> Vector{AttrSpec}

Attributes of a plain physical quantity: its `type` and its `units`.
"""
quantity(units::AbstractString) = [
    _type_attr("quantity"),
    AttrSpec(:units, :str; values = [units], default = units, doc = "Units of the quantity"),
]

"""
    sky_coord(; frame = "icrs") -> Vector{AttrSpec}

Attributes of a direction on the sky. `frame` names the astropy reference frame
the direction is expressed in. The units are radians for a direction and metres
for a distance along the line of sight, which is why both are permitted.
"""
sky_coord(; frame::AbstractString = "icrs") = [
    _type_attr("sky_coord"),
    AttrSpec(:units, :str; values = ["m", "rad"], default = "rad"),
    AttrSpec(:frame, :str; values = SKY_FRAMES, default = frame, doc = "Astropy reference frame of the direction"),
]

"""
    time_measure(; type = "time", integration_time = false) -> Vector{AttrSpec}

Attributes of a time axis: the astropy `scale` and `format` needed to interpret
the stored numbers as absolute times. The main dataset's `time` coordinate also
carries its nominal sampling interval as an `integration_time` quantity.

Sub-datasets sampled on their own clock — weather, pointing, ephemerides —
name their measure after that axis rather than plain `"time"`, so that a
dataset carrying both cannot confuse the two.
"""
function time_measure(; type::AbstractString = "time", integration_time::Bool = false)
    attrs = [
        _type_attr(type),
        AttrSpec(:units, :str; values = ["s"], default = "s"),
        AttrSpec(:scale, :str; values = TIME_SCALES, default = "utc", doc = "Astropy time scale"),
        AttrSpec(:format, :str; values = TIME_FORMATS, default = "unix", doc = "Astropy time format"),
    ]
    if integration_time
        push!(
            attrs,
            AttrSpec(
                :integration_time, :dataarray;
                doc = "Nominal sampling interval, in seconds"
            )
        )
    end
    return attrs
end

"""
    location(; coordinate_system = "geocentric") -> Vector{AttrSpec}

Attributes of an absolute position in space.

`coordinate_system` says how the position is parameterized, which in turn
decides which of the label coordinates it is indexed by and whether its units
are metres or radians. A geodetic or planetodetic position also needs the
reference `ellipsoid` its latitude is measured against.

See [`offset`](@ref) for a position measured relative to another.
"""
location(; coordinate_system::AbstractString = "geocentric") = [
    _type_attr("location"),
    AttrSpec(:units, :str; values = ["m", "rad"], default = "m"),
    AttrSpec(:frame, :str; values = ["ITRS", "Undefined"], default = "ITRS"),
    AttrSpec(:coordinate_system, :str; values = COORDINATE_SYSTEMS, default = coordinate_system),
    AttrSpec(:origin_object_name, :str; default = "earth", doc = "Body the position is measured on"),
    AttrSpec(
        :ellipsoid, :str; values = ELLIPSOIDS, default = "GRS80", optional = true,
        doc = "Reference ellipsoid, required for a geodetic or planetodetic position"
    ),
]

"""
    offset(origin) -> Vector{AttrSpec}

Attributes of a position measured relative to another array rather than to the
centre of a body. `origin` names that array — for a phased array's elements,
the station's own `ANTENNA_POSITION`.
"""
offset(origin::AbstractString) = [
    _type_attr("location"),
    AttrSpec(:units, :str; values = ["m"], default = "m"),
    AttrSpec(:coordinate_system, :str; values = ["topocentric"], default = "topocentric"),
    AttrSpec(:origin, :str; values = [origin], default = origin, doc = "Array the offset is measured from"),
]

"""
    spectral_coord(; observer = "icrs") -> Vector{AttrSpec}

Attributes of a frequency: its units and the reference frame (`observer`) it is
expressed in. This is all a frequency-valued *variable* needs; the frequency
*axis* of a dataset carries more, for which see [`frequency_axis`](@ref).
"""
spectral_coord(; observer::AbstractString = "icrs") = [
    _type_attr("spectral_coord"),
    AttrSpec(:units, :str; values = ["Hz"], default = "Hz"),
    AttrSpec(
        :observer, :str; values = SPECTRAL_OBSERVERS, default = observer,
        doc = "Reference frame the frequencies are expressed in"
    ),
]

"""
    frequency_axis(; observer = "icrs") -> Vector{AttrSpec}

Attributes of a dataset's frequency axis.

Beyond the [`spectral_coord`](@ref) attributes, the axis describes the single
spectral window it came from — its name and intents, the nominal channel width,
and the reference frequency and doppler shift needed to move the frequencies
between frames. MSv4 keeps this here rather than in a separate table because
one dataset covers exactly one spectral window.
"""
frequency_axis(; observer::AbstractString = "icrs") = vcat(
    spectral_coord(; observer),
    AttrSpec[
        long_name("Frequency"),
        AttrSpec(
            :spectral_window_name, :str;
            doc = "Name of the spectral window, which replaces MSv2's data description index"
        ),
        AttrSpec(
            :spectral_window_intents, Symbol("list[str]");
            doc = "What the spectral window was set up to observe"
        ),
        AttrSpec(:frequency_group_name, :str; optional = true, doc = "Group the spectral window belongs to"),
        AttrSpec(:reference_frequency, :dataarray; doc = "Frequency the doppler shift is defined against"),
        AttrSpec(:channel_width, :dataarray; doc = "Nominal channel width, in hertz"),
        AttrSpec(:doppler, :dataarray; optional = true, doc = "Doppler shift applied to the frequencies"),
    ],
)

"""
    doppler(; units = "m/s", doppler_type = "radio") -> Vector{AttrSpec}

Attributes of a doppler shift. `doppler_type` names the velocity convention.
"""
doppler(; units::AbstractString = "m/s", doppler_type::AbstractString = "radio") = [
    _type_attr("doppler"),
    AttrSpec(:units, :str; values = ["m/s", "ratio"], default = units),
    AttrSpec(:doppler_type, :str; values = DOPPLER_TYPES, default = doppler_type),
]

"""
    uvw_measure() -> Vector{AttrSpec}

Attributes of the `UVW` baseline coordinates.
"""
uvw_measure() = [
    _type_attr("uvw"),
    AttrSpec(:units, :str; values = ["m"], default = "m"),
    AttrSpec(:frame, :str; values = UVW_FRAMES, default = "icrs"),
    AttrSpec(:long_name, :str; default = "Baseline coordinates", optional = true),
    AttrSpec(:allow_multiple_versions, :bool; default = true, optional = true),
]

"""
    rotation_matrix() -> Vector{AttrSpec}

Attributes of a dimensionless rotation matrix.
"""
rotation_matrix() = [
    _type_attr("rotation_matrix"),
    AttrSpec(:units, :str; values = ["dimensionless"], default = "dimensionless"),
]

"""
    long_name(default) -> AttrSpec

The optional human-readable label MSv4 attaches to most arrays.
"""
long_name(default::AbstractString) = AttrSpec(:long_name, :str; default, optional = true)

"""
    versioned() -> AttrSpec

Marks an array that data groups may carry several versions of, such as
`VISIBILITY` alongside `VISIBILITY_CORRECTED`.
"""
versioned() = AttrSpec(:allow_multiple_versions, :bool; default = true, optional = true)

"""
    time_sampling(name) -> Vector{AttrSpec}

Attributes of a per-sample time such as `TIME_CENTROID`: enough to interpret it
as an absolute time, but without the nominal interval a time *axis* carries.
"""
time_sampling(name::AbstractString = "Time sampling data") = [
    AttrSpec(:units, :str; values = ["s"], default = "s"),
    AttrSpec(:scale, :str; values = TIME_SCALES, default = "utc"),
    AttrSpec(:format, :str; values = TIME_FORMATS, default = "unix"),
    long_name(name),
]

"""
    frequency_sampling(name) -> Vector{AttrSpec}

Attributes of a per-sample frequency such as `FREQUENCY_CENTROID` or
`EFFECTIVE_CHANNEL_WIDTH`.
"""
frequency_sampling(name::AbstractString = "Frequency sampling data") = [
    AttrSpec(:units, :str; values = ["Hz"], default = "Hz"),
    AttrSpec(:observer, :str; values = SPECTRAL_OBSERVERS, default = "icrs"),
    long_name(name),
]

# The bit meanings MSv4 assigns to an unsigned FLAG array, in bit order.
const FLAG_BITS = [
    "UNSPECIFIED_BIT", "STATIC_BIT", "CAM_BIT", "DATA_LOST_BIT",
    "INGEST_RFI_BIT", "PREDICTED_RFI_BIT", "CAL_RFI_BIT", "POSTPROC_BIT",
]
