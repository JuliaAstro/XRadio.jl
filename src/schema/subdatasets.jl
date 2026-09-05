"""
The sub-datasets of an MSv4.

Everything needed to interpret the correlated data that is not the correlated
data itself: the antennas, the field and source, and whatever calibration and
environment measurements the telescope recorded. Each is optional except
`antenna`, and each hangs off the Measurement Set as a named child.

Several of these describe the same antennas, so they repeat the same identity
coordinates. `_antenna_identity` builds that block once.
"""

# antenna_name and the per-antenna labels that travel with it. `mount` and
# `telescope_name` are absent from the system calibration dataset.
function _antenna_identity(; mount::Bool = true, telescope::Bool = true)
    out = ArraySpec[
        ArraySpec(
            "antenna_name", (:antenna_name,), String;
            attrs = [long_name("Antenna name")], doc = "Name of each antenna."
        ),
        ArraySpec(
            "station_name", (:antenna_name,), String;
            doc = "Name of the station each antenna occupies. Antennas move between stations; \
                   the pair identifies a physical location."
        ),
    ]
    mount && push!(
        out,
        ArraySpec(
            "mount", (:antenna_name,), String;
            doc = "Mount type, which fixes how parallactic angle enters the feed response: \
                   `ALT-AZ`, `EQUATORIAL`, `X-Y`, `ALT-AZ+NASMYTH-R`, `ALT-AZ+NASMYTH-L`."
        )
    )
    telescope && push!(
        out,
        ArraySpec(
            "telescope_name", (:antenna_name,), String;
            doc = "Telescope each antenna belongs to. Heterogeneous VLBI arrays mix several."
        )
    )
    append!(
        out,
        ArraySpec[
            ArraySpec(
                "receptor_label", (:receptor_label,), String;
                doc = "Labels of the feed's receptors."
            ),
            ArraySpec(
                "polarization_type", (:antenna_name, :receptor_label), String;
                doc = "Polarization each receptor responds to, such as `R`/`L` or `X`/`Y`."
            ),
        ],
    )
    return out
end

# The label coordinates that say how a position variable is parameterized.
# A position is written in whichever of these its `coordinate_system` calls for.
const _POSITION_LABELS = ArraySpec[
    ArraySpec(
        "cartesian_pos_label", (:cartesian_pos_label,), String; optional = true,
        doc = "Labels of a Cartesian position: `x`, `y`, `z`."
    ),
    ArraySpec(
        "ellipsoid_dir_label", (:ellipsoid_dir_label,), String; optional = true,
        doc = "Labels of an ellipsoidal direction: `lon`, `lat`."
    ),
    ArraySpec(
        "ellipsoid_dis_label", (:ellipsoid_dis_label,), String; optional = true,
        doc = "Label of an ellipsoidal distance: `dist`."
    ),
]

# The dimension tuples a position variable may span, one per parameterization.
_position_dims(lead::Symbol...) = (
    (lead..., :ellipsoid_dir_label),
    (lead..., :ellipsoid_dis_label),
    (lead..., :cartesian_pos_label),
)

"""
    ANTENNA_SCHEMA

Where the antennas are and what they are. Combines MSv2's `ANTENNA` and `FEED`
tables, which MSv4 can merge because one MSv4 admits only a single feed type
per antenna.
"""
const ANTENNA_SCHEMA = DatasetSpec(
    :antenna, "antenna";
    coords = vcat(_antenna_identity(), _POSITION_LABELS),
    vars = ArraySpec[
        ArraySpec(
            "ANTENNA_POSITION", _position_dims(:antenna_name), Float64;
            attrs = location(),
            doc = "Position of each antenna, in the parameterization its `coordinate_system` names."
        ),
        ArraySpec(
            "ANTENNA_DISH_DIAMETER", (:antenna_name,), Float64;
            attrs = quantity("m"), optional = true,
            doc = "Physical dish diameter."
        ),
        ArraySpec(
            "ANTENNA_EFFECTIVE_DISH_DIAMETER", (:antenna_name,), Float64;
            attrs = quantity("m"), optional = true,
            doc = "Diameter of an unblocked dish of the same collecting area."
        ),
        ArraySpec(
            "ANTENNA_BLOCKAGE", (:antenna_name,), Float64;
            attrs = quantity("m"), optional = true,
            doc = "Diameter blocked by the feed support structure."
        ),
        ArraySpec(
            "ANTENNA_RECEPTOR_ANGLE", (:antenna_name, :receptor_label), Float64;
            attrs = quantity("rad"), optional = true,
            doc = "Rotation of each receptor about the feed axis."
        ),
        ArraySpec(
            "ANTENNA_FOCUS_LENGTH", (:antenna_name,), Float64;
            attrs = quantity("m"), optional = true,
            doc = "Focal length of the dish."
        ),
    ],
    attrs = AttrSpec[
        AttrSpec(
            :overall_telescope_name, :str; optional = true,
            doc = "Name for the array as a whole, where the individual telescopes differ."
        ),
        AttrSpec(
            :relocatable_antennas, :bool; optional = true,
            doc = "Whether antennas move between stations, as at the VLA."
        ),
        AttrSpec(:type, :str; values = ["antenna"], default = "antenna"),
    ],
    doc = "Antenna positions and feed properties. Replaces MSv2's ANTENNA and FEED tables."
)

"""
    POINTING_SCHEMA

Where each antenna was actually pointed, as a function of time.

The direction may be sampled on the main `time` axis or on the telescope's own
`time_pointing` clock, and may be given either directly or as polynomial
coefficients over `n_polynomial`.
"""
const POINTING_SCHEMA = DatasetSpec(
    :pointing, "pointing";
    coords = ArraySpec[
        ArraySpec(
            "antenna_name", (:antenna_name,), String;
            attrs = [long_name("Antenna name")], doc = "Name of each antenna."
        ),
        ArraySpec(
            "local_sky_dir_label", (:local_sky_dir_label,), String;
            doc = "Labels of a local sky direction: `az`, `alt`."
        ),
        ArraySpec(
            "time", (:time,), Float64;
            attrs = time_measure(), optional = true,
            doc = "Pointing sampled on the correlated data's own time axis."
        ),
        ArraySpec(
            "time_pointing", (:time_pointing,), Float64;
            attrs = time_measure(; type = "time_pointing"), optional = true,
            doc = "Pointing sampled on the telescope's own, usually faster, clock."
        ),
        ArraySpec(
            "n_polynomial", (:n_polynomial,), Int64; optional = true,
            doc = "Polynomial term index, when the pointing is stored as coefficients."
        ),
    ],
    vars = ArraySpec[
        ArraySpec(
            "POINTING_BEAM",
            (
                (:time, :antenna_name, :local_sky_dir_label),
                (:time_pointing, :antenna_name, :local_sky_dir_label),
                (:time, :antenna_name, :local_sky_dir_label, :n_polynomial),
                (:time_pointing, :antenna_name, :local_sky_dir_label, :n_polynomial),
            ),
            Float64;
            attrs = sky_coord(),
            doc = "Commanded beam direction of each antenna."
        ),
        ArraySpec(
            "POINTING_DISH_MEASURED",
            (
                (:time, :antenna_name, :local_sky_dir_label),
                (:time_pointing, :antenna_name, :local_sky_dir_label),
            ),
            Float64;
            attrs = sky_coord(), optional = true,
            doc = "Direction the dish encoders actually reported."
        ),
        ArraySpec(
            "POINTING_OVER_THE_TOP",
            ((:time, :antenna_name), (:time_pointing, :antenna_name)), Bool;
            optional = true,
            doc = "True where the antenna was pointed past the zenith."
        ),
    ],
    attrs = [AttrSpec(:type, :str; values = ["pointing"], default = "pointing")],
    doc = "Antenna pointing directions over time. Replaces MSv2's POINTING table."
)

"""
    WEATHER_SCHEMA

Ground weather and ionospheric conditions at each station, used to model
tropospheric and ionospheric delay.
"""
const WEATHER_SCHEMA = let
    tdims = ((:station_name, :time), (:station_name, :time_weather))
    DatasetSpec(
        :weather, "weather";
        coords = vcat(
            ArraySpec[
                ArraySpec(
                    "station_name", (:station_name,), String;
                    doc = "Station the measurements were taken at."
                ),
                ArraySpec(
                    "time", (:time,), Float64;
                    attrs = time_measure(), optional = true,
                    doc = "Weather sampled on the correlated data's own time axis."
                ),
                ArraySpec(
                    "time_weather", (:time_weather,), Float64;
                    attrs = time_measure(; type = "time_weather"), optional = true,
                    doc = "Weather sampled on the weather station's own clock."
                ),
            ],
            _POSITION_LABELS,
        ),
        vars = ArraySpec[
            ArraySpec(
                "STATION_POSITION", _position_dims(:station_name), Float64;
                attrs = location(),
                doc = "Position of each weather station."
            ),
            ArraySpec(
                "H2O", tdims, Float64; attrs = quantity("/m^2"), optional = true,
                doc = "Precipitable water vapour column density."
            ),
            ArraySpec(
                "IONOS_ELECTRON", tdims, Float64; attrs = quantity("/m^2"), optional = true,
                doc = "Total electron content of the ionosphere."
            ),
            ArraySpec(
                "PRESSURE", tdims, Float64; attrs = quantity("Pa"), optional = true,
                doc = "Ambient atmospheric pressure."
            ),
            ArraySpec(
                "REL_HUMIDITY", tdims, Float64; optional = true,
                doc = "Ambient relative humidity, as a fraction."
            ),
            ArraySpec(
                "TEMPERATURE", tdims, Float64; attrs = quantity("K"), optional = true,
                doc = "Ambient air temperature."
            ),
            ArraySpec(
                "DEW_POINT", tdims, Float64; attrs = quantity("K"), optional = true,
                doc = "Temperature at which the air would saturate."
            ),
            ArraySpec(
                "WIND_DIRECTION", tdims, Float64; attrs = quantity("rad"), optional = true,
                doc = "Direction the wind blows from, north through east."
            ),
            ArraySpec(
                "WIND_SPEED", tdims, Float64; attrs = quantity("m/s"), optional = true,
                doc = "Wind speed."
            ),
        ],
        attrs = [AttrSpec(:type, :str; values = ["weather"], default = "weather")],
        doc = "Weather and ionospheric conditions. Replaces MSv2's WEATHER table."
    )
end

"""
    SYSTEM_CALIBRATION_SCHEMA

Receiver and sky temperatures measured by the telescope, from which the
absolute flux scale is derived. Every temperature may be sampled on either time
axis and at either frequency resolution, or not resolved in frequency at all.
"""
const SYSTEM_CALIBRATION_SCHEMA = let
    # Every temperature admits the same six shapes: two time axes crossed with
    # full frequency resolution, the coarser system-calibration frequency
    # channels, or no frequency dependence.
    tdims = [
        (t..., f...) for t in ((:antenna_name, :time_system_cal), (:antenna_name, :time))
            for f in ((:receptor_label, :frequency_system_cal), (:receptor_label, :frequency), (:receptor_label,))
    ]
    temperature(name, doc) = ArraySpec(
        name, tdims, Float64; attrs = quantity("K"), optional = true, doc
    )
    DatasetSpec(
        :system_calibration, "system_calibration";
        coords = vcat(
            _antenna_identity(; mount = false, telescope = false),
            ArraySpec[
                ArraySpec(
                    "time", (:time,), Float64;
                    attrs = time_measure(), optional = true,
                    doc = "Measurements sampled on the correlated data's own time axis."
                ),
                ArraySpec(
                    "time_system_cal", (:time_system_cal,), Float64;
                    attrs = time_measure(), optional = true,
                    doc = "Measurements sampled on the calibration system's own clock."
                ),
                ArraySpec(
                    "frequency", (:frequency,), Float64;
                    attrs = frequency_axis(), optional = true,
                    doc = "Measurements resolved at the correlated data's channel resolution."
                ),
                ArraySpec(
                    "frequency_system_cal", (:frequency_system_cal,), Int64; optional = true,
                    doc = "Measurements resolved at the calibration system's own, coarser, resolution."
                ),
            ],
        ),
        vars = ArraySpec[
            ArraySpec(
                "PHASE_DIFFERENCE",
                ((:antenna_name, :time_system_cal), (:antenna_name, :time)), Float64;
                attrs = quantity("rad"), optional = true,
                doc = "Phase difference between the two receptors of the feed."
            ),
            temperature("TCAL", "Temperature of the calibration noise source."),
            temperature("TRX", "Receiver noise temperature."),
            temperature("TSKY", "Sky brightness temperature."),
            temperature("TSYS", "Total system temperature."),
            temperature("TANT", "Antenna temperature."),
            temperature("TANT_SYS", "Ratio of antenna to system temperature."),
        ],
        attrs = [AttrSpec(:type, :str; values = ["system_calibration"], default = "system_calibration")],
        doc = "Receiver and sky temperatures. Replaces MSv2's SYSCAL table."
    )
end

"""
    GAIN_CURVE_SCHEMA

Antenna gain as a polynomial in elevation, together with the sensitivity needed
to convert it to a flux scale. Used mainly for VLBI, where no single calibrator
is visible to the whole array.
"""
const GAIN_CURVE_SCHEMA = DatasetSpec(
    :gain_curve, "gain_curve";
    coords = vcat(
        _antenna_identity(),
        ArraySpec[
            ArraySpec(
                "gain_curve_type", (:antenna_name,), String; optional = true,
                doc = "Whether the polynomial is in elevation or in zenith angle."
            ),
            ArraySpec(
                "poly_term", (:poly_term,), Int64;
                doc = "Index of each polynomial coefficient."
            ),
        ],
    ),
    vars = ArraySpec[
        ArraySpec(
            "GAIN_CURVE", (:antenna_name, :poly_term, :receptor_label), Float64;
            doc = "Polynomial coefficients of the gain against elevation."
        ),
        ArraySpec(
            "GAIN_CURVE_INTERVAL", (:antenna_name,), Float64;
            attrs = quantity("s"),
            doc = "Time over which the gain curve is taken to hold."
        ),
        ArraySpec(
            "GAIN_CURVE_SENSITIVITY", (:antenna_name, :receptor_label), Float64;
            attrs = quantity("K/Jy"),
            doc = "Degrees per jansky, converting the gain curve to a flux scale."
        ),
    ],
    attrs = AttrSpec[
        AttrSpec(:measured_date, :str; doc = "ISO 8601 date the gain curve was measured."),
        AttrSpec(:type, :str; values = ["gain_curve"], default = "gain_curve"),
    ],
    doc = "Antenna gain against elevation. Replaces MSv2's GAIN_CURVE table."
)

"""
    PHASE_CALIBRATION_SCHEMA

Injected phase-calibration tones and the measured cable delay, used in VLBI to
track the instrumental phase and delay of each station.
"""
const PHASE_CALIBRATION_SCHEMA = let
    tone_dims = (
        (:antenna_name, :time, :receptor_label, :tone_label),
        (:antenna_name, :time_phase_cal, :receptor_label, :tone_label),
    )
    scalar_dims = ((:antenna_name, :time), (:antenna_name, :time_phase_cal))
    DatasetSpec(
        :phase_calibration, "phase_calibration";
        coords = vcat(
            _antenna_identity(),
            ArraySpec[
                ArraySpec(
                    "time", (:time,), Float64;
                    attrs = time_measure(), optional = true,
                    doc = "Tones sampled on the correlated data's own time axis."
                ),
                ArraySpec(
                    "time_phase_cal", (:time_phase_cal,), Float64; optional = true,
                    doc = "Tones sampled on the phase-calibration system's own clock."
                ),
                ArraySpec(
                    "tone_label", (:tone_label,), String; optional = true,
                    doc = "Label of each injected tone."
                ),
            ],
        ),
        vars = ArraySpec[
            ArraySpec(
                "PHASE_CAL", tone_dims, Union{ComplexF32, ComplexF64};
                doc = "Complex amplitude measured for each injected tone."
            ),
            ArraySpec(
                "PHASE_CAL_CABLE_CAL", scalar_dims, Float64; attrs = quantity("s"),
                doc = "Measured round-trip delay of the cable feeding the station."
            ),
            ArraySpec(
                "PHASE_CAL_INTERVAL", scalar_dims, Float64; attrs = quantity("s"),
                doc = "Time over which each tone measurement was integrated."
            ),
            ArraySpec(
                "PHASE_CAL_TONE_FREQUENCY", tone_dims, Float64; attrs = quantity("Hz"),
                doc = "Frequency of each injected tone."
            ),
        ],
        attrs = [AttrSpec(:type, :str; values = ["phase_calibration"], default = "phase_calibration")],
        doc = "Injected phase-calibration tones and cable delay, for VLBI."
    )
end

"""
    PHASED_ARRAY_SCHEMA

Layout of a phased-array station: where each element sits, in the station's own
local frame, and how that frame relates to the geocentric one.
"""
const PHASED_ARRAY_SCHEMA = DatasetSpec(
    :phased_array, "phased_array";
    coords = ArraySpec[
        ArraySpec(
            "antenna_name", (:antenna_name,), String;
            attrs = [long_name("Antenna name")], doc = "Name of each station."
        ),
        ArraySpec(
            "element_id", (:element_id,), Union{Int32, Int64};
            doc = "Identifier of each element within its station."
        ),
        ArraySpec("receptor_label", (:receptor_label,), String; doc = "Labels of the feed's receptors."),
        ArraySpec(
            "polarization_type", (:antenna_name, :receptor_label), String;
            doc = "Polarization each receptor responds to."
        ),
        ArraySpec(
            "cartesian_pos_label", (:cartesian_pos_label,), String;
            doc = "Labels of a geocentric Cartesian position: `x`, `y`, `z`."
        ),
        ArraySpec(
            "cartesian_pos_label_local", (:cartesian_pos_label_local,), String;
            doc = "Labels of the station's local Cartesian frame: `p`, `q`, `r`."
        ),
    ],
    vars = ArraySpec[
        ArraySpec(
            "PHASED_ARRAY_COORDINATE_AXES",
            (:antenna_name, :cartesian_pos_label_local, :cartesian_pos_label), Float64;
            attrs = rotation_matrix(),
            doc = "Rotation taking the station's local frame to the geocentric one."
        ),
        ArraySpec(
            "PHASED_ARRAY_ELEMENT_OFFSET",
            (:antenna_name, :cartesian_pos_label_local, :element_id), Float64;
            attrs = offset("ANTENNA_POSITION"),
            doc = "Offset of each element from the station position, in the local frame."
        ),
        ArraySpec(
            "PHASED_ARRAY_ELEMENT_FLAG",
            (:antenna_name, :receptor_label, :element_id), Bool;
            doc = "True where an element was not contributing."
        ),
    ],
    attrs = [AttrSpec(:type, :str; values = ["phased_array"], default = "phased_array")],
    doc = "Element layout of a phased-array station, such as LOFAR or SKA-Low."
)

# Attributes shared by both field-and-source datasets.
_field_source_attrs(type, ephemeris::Bool) = begin
    attrs = AttrSpec[
        AttrSpec(
            :doppler_shift_velocity, :str; values = ["ratio", "m/s"], optional = true,
            doc = "Units of the velocity the frequency axis was doppler-shifted by."
        ),
        AttrSpec(:source_model_url, :str; optional = true, doc = "Where the source model came from."),
    ]
    ephemeris && push!(
        attrs,
        AttrSpec(:ephemeris_name, :str; optional = true, doc = "Name of the ephemeris used.")
    )
    push!(attrs, AttrSpec(:type, :str; values = [type], default = type))
    return attrs
end

"""
    FIELD_SOURCE_SCHEMA

What was being observed and where it is. Merges MSv2's `FIELD` and `SOURCE`
tables, which MSv4 can combine because a data group points at exactly one of
these datasets.
"""
const FIELD_SOURCE_SCHEMA = DatasetSpec(
    :field_and_source, "field_and_source";
    coords = ArraySpec[
        ArraySpec(
            "source_name", (:field_name,), String;
            doc = "Astronomical object observed in each field."
        ),
        ArraySpec("field_name", (:field_name,), String; doc = "Name of each field."),
        ArraySpec(
            "sky_dir_label", (:sky_dir_label,), String;
            doc = "Labels of a sky direction: `ra`, `dec`."
        ),
        ArraySpec(
            "cartesian_pos_label", (:cartesian_pos_label,), String; optional = true,
            doc = "Labels of a Cartesian position: `x`, `y`, `z`."
        ),
        ArraySpec(
            "line_label", (:line_label,), String; optional = true,
            doc = "Label of each spectral line of interest."
        ),
        ArraySpec(
            "line_names", (:field_name, :line_label), String; optional = true,
            doc = "Name of each spectral line."
        ),
    ],
    vars = ArraySpec[
        ArraySpec(
            "FIELD_PHASE_CENTER_DIRECTION", (:field_name, :sky_dir_label), Float64;
            attrs = sky_coord(), optional = true,
            doc = "Direction the visibilities are phased to."
        ),
        ArraySpec(
            "FIELD_REFERENCE_CENTER_DIRECTION", (:field_name, :sky_dir_label), Float64;
            attrs = sky_coord(), optional = true,
            doc = "Direction the field's coordinate grid is referenced to."
        ),
        ArraySpec(
            "SOURCE_DIRECTION", (:field_name, :sky_dir_label), Float64;
            attrs = sky_coord(), optional = true,
            doc = "Direction of the source itself, which need not be the phase center."
        ),
        ArraySpec(
            "LINE_REST_FREQUENCY", (:field_name, :line_label), Float64;
            attrs = spectral_coord(), optional = true,
            doc = "Rest frequency of each spectral line."
        ),
        ArraySpec(
            "LINE_SYSTEMIC_VELOCITY", (:field_name, :line_label), Float64;
            attrs = quantity("m/s"), optional = true,
            doc = "Systemic velocity of the source at each line."
        ),
        ArraySpec(
            "OBSERVER_POSITION", (:cartesian_pos_label,), Float64;
            attrs = location(), optional = true,
            doc = "Position the observation is referenced to."
        ),
    ],
    attrs = _field_source_attrs("field_and_source", false),
    doc = "Field and source being observed. Replaces MSv2's FIELD and SOURCE tables."
)

"""
    FIELD_SOURCE_EPHEMERIS_SCHEMA

The field and source of a moving target, whose position is tabulated against
time rather than fixed. Built from JPL Horizons ephemerides rather than from
MSv2's own ephemeris tables.
"""
const FIELD_SOURCE_EPHEMERIS_SCHEMA = let
    # A direction on the sky, sampled on either clock, given either as a
    # direction or as a distance.
    sky_dims = ((:time, :sky_dir_label), (:time, :sky_dis_label))
    sky_dims_both = (
        (:time, :sky_dir_label), (:time, :sky_dis_label),
        (:time_ephemeris, :sky_dir_label), (:time_ephemeris, :sky_dis_label),
    )
    scalar_dims = ((:time,), (:time_ephemeris,))
    body_dims = (_position_dims(:time)..., _position_dims(:time_ephemeris)...)
    DatasetSpec(
        :field_and_source_ephemeris, "field_and_source_ephemeris";
        coords = ArraySpec[
            ArraySpec("source_name", (:time,), String; doc = "Object observed at each time."),
            ArraySpec("field_name", (:time,), String; doc = "Field observed at each time."),
            ArraySpec(
                "time", (:time,), Float64; attrs = time_measure(),
                doc = "The correlated data's own time axis."
            ),
            ArraySpec(
                "sky_dir_label", (:sky_dir_label,), String; optional = true,
                doc = "Labels of a sky direction: `ra`, `dec`."
            ),
            ArraySpec(
                "sky_dis_label", (:sky_dis_label,), String; optional = true,
                doc = "Label of a distance along the line of sight: `dist`."
            ),
            ArraySpec(
                "ellipsoid_dir_label", (:ellipsoid_dir_label,), String; optional = true,
                doc = "Labels of an ellipsoidal direction on the body: `lon`, `lat`."
            ),
            ArraySpec(
                "ellipsoid_dis_label", (:ellipsoid_dis_label,), String; optional = true,
                doc = "Label of an ellipsoidal distance: `dist`."
            ),
            ArraySpec(
                "cartesian_pos_label", (:cartesian_pos_label,), String; optional = true,
                doc = "Labels of a Cartesian position: `x`, `y`, `z`."
            ),
            ArraySpec(
                "line_label", (:line_label,), String; optional = true,
                doc = "Label of each spectral line of interest."
            ),
            ArraySpec(
                "line_names", (:field_name, :line_label), String; optional = true,
                doc = "Name of each spectral line."
            ),
            ArraySpec(
                "time_ephemeris", (:time_ephemeris,), Float64;
                attrs = time_measure(; type = "time_ephemeris"), optional = true,
                doc = "The ephemeris table's own, usually coarser, time axis."
            ),
        ],
        vars = ArraySpec[
            ArraySpec(
                "FIELD_PHASE_CENTER_DIRECTION", sky_dims, Float64;
                attrs = sky_coord(), optional = true,
                doc = "Direction the visibilities are phased to."
            ),
            ArraySpec(
                "FIELD_PHASE_CENTER_DISTANCE", sky_dims, Float64;
                attrs = sky_coord(), optional = true,
                doc = "Distance to the phase center, for a nearby body."
            ),
            ArraySpec(
                "FIELD_REFERENCE_CENTER_DIRECTION", sky_dims, Float64;
                attrs = sky_coord(), optional = true,
                doc = "Direction the field's coordinate grid is referenced to."
            ),
            ArraySpec(
                "FIELD_REFERENCE_CENTER_DISTANCE", sky_dims, Float64;
                attrs = sky_coord(), optional = true,
                doc = "Distance to the reference center."
            ),
            ArraySpec(
                "LINE_REST_FREQUENCY", (:time, :line_label), Float64;
                attrs = spectral_coord(), optional = true,
                doc = "Rest frequency of each spectral line."
            ),
            ArraySpec(
                "LINE_SYSTEMIC_VELOCITY", (:time, :line_label), Float64;
                attrs = quantity("m/s"), optional = true,
                doc = "Systemic velocity of the source at each line."
            ),
            ArraySpec(
                "SOURCE_LOCATION", sky_dims_both, Float64;
                attrs = sky_coord(), optional = true,
                doc = "Tabulated position of the moving body."
            ),
            ArraySpec(
                "SOURCE_RADIAL_VELOCITY", scalar_dims, Float64;
                attrs = quantity("m/s"), optional = true,
                doc = "Radial velocity of the body relative to the observer."
            ),
            ArraySpec(
                "NORTH_POLE_POSITION_ANGLE", scalar_dims, Float64;
                attrs = quantity("rad"), optional = true,
                doc = "Position angle of the body's north pole, east of north."
            ),
            ArraySpec(
                "NORTH_POLE_ANGULAR_DISTANCE", scalar_dims, Float64;
                attrs = quantity("rad"), optional = true,
                doc = "Angular distance from the body's north pole to its disc center."
            ),
            ArraySpec(
                "SUB_OBSERVER_DIRECTION", body_dims, Float64;
                attrs = location(), optional = true,
                doc = "Point on the body directly beneath the observer."
            ),
            ArraySpec(
                "SUB_SOLAR_POSITION", body_dims, Float64;
                attrs = location(), optional = true,
                doc = "Point on the body directly beneath the sun."
            ),
            ArraySpec(
                "HELIOCENTRIC_RADIAL_VELOCITY", scalar_dims, Float64;
                attrs = quantity("m/s"), optional = true,
                doc = "Radial velocity of the body relative to the sun."
            ),
            ArraySpec(
                "OBSERVER_PHASE_ANGLE", scalar_dims, Float64;
                attrs = quantity("rad"), optional = true,
                doc = "Sun–body–observer angle."
            ),
            ArraySpec(
                "OBSERVER_POSITION", _position_dims(), Float64;
                attrs = location(), optional = true,
                doc = "Position the observation is referenced to."
            ),
        ],
        attrs = _field_source_attrs("field_and_source_ephemeris", true),
        doc = "Field and source of a moving target, tabulated against time."
    )
end
