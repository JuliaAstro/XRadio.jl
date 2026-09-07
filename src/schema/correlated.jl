"""
The correlated (main) datasets.

An MSv4 holds one of these at its root: [`VISIBILITY_SCHEMA`](@ref) for an
interferometer, [`SPECTRUM_SCHEMA`](@ref) for a single dish. The two differ
only in their science array and in whether the second dimension counts
baselines or antennas, so most of their content is built here once.

A single dataset covers one observation, spectral window, polarization setup,
observation mode, processor and beam per antenna. That restriction is what lets
the data be dense rectangular arrays whose shape does not change with time;
baselines that drop out are padded with `NaN` rather than removed.
"""

# Attributes every correlated dataset carries.
_correlated_attrs(types) = [
    AttrSpec(:observation_info, :dict; doc = "Observation-level metadata; see [`OBSERVATION_INFO`](@ref)"),
    AttrSpec(:processor_info, :dict; doc = "Correlator metadata; see [`PROCESSOR_INFO`](@ref)"),
    AttrSpec(
        :data_groups, :dict;
        doc = "Groupings of correlated data with their flags and weights; see [`DATA_GROUP`](@ref). \
               A group named `base` is mandatory"
    ),
    AttrSpec(:schema_version, :str; doc = "MSv4 schema version the dataset was written against"),
    AttrSpec(:creator, :dict; doc = "Software that wrote the dataset; see [`CREATOR`](@ref)"),
    AttrSpec(:creation_date, :str; doc = "ISO 8601 time at which the dataset was written"),
    AttrSpec(:type, :str; values = types, default = first(types)),
]

# Coordinates shared by both correlated datasets. `second` is the dimension
# that indexes the array elements being correlated: baselines for an
# interferometer, antennas for a single dish.
function _correlated_coords(second::Symbol)
    return ArraySpec[
        ArraySpec(
            "time", (:time,), Float64;
            attrs = time_measure(; integration_time = true),
            doc = "Mid-point of each nominal sampling interval, whose length is the \
                   `integration_time` attribute."
        ),
        ArraySpec(
            "frequency", (:frequency,), Float64;
            attrs = frequency_axis(),
            doc = "Channel centre frequencies of the dataset's single spectral window."
        ),
        ArraySpec(
            "polarization", (:polarization,), String;
            attrs = [long_name("Polarization")],
            doc = "Polarization products, named in the Stokes convention: `XX`, `RL`, `I`, and so on."
        ),
        ArraySpec(
            "field_name", (:time,), String;
            doc = "Field observed at each time step, tying the data to `field_and_source`."
        ),
        ArraySpec(
            "scan_name", (:time,), String;
            attrs = [
                AttrSpec(
                    :scan_intents, Symbol("list[str]");
                    doc = "What each scan was observed for, such as `OBSERVE_TARGET#ON_SOURCE`"
                ),
                long_name("Scan name"),
            ],
            doc = "Scan each time step belongs to. Scans are named, not numbered, so that \
                   selecting a subset never requires renumbering."
        ),
        ArraySpec(
            "polarization_mixed", (second, :polarization), String;
            optional = true,
            doc = "Per-element polarization products, for arrays whose antennas do not share \
                   a polarization setup."
        ),
    ]
end

# Data variables shared by both correlated datasets, given the element
# dimension and the shape of the science arrays.
function _correlated_vars(second::Symbol, shape)
    return ArraySpec[
        ArraySpec(
            "FLAG", shape, Union{Bool, UInt8, UInt16, UInt32, UInt64};
            attrs = [
                AttrSpec(
                    :flag_bits, Symbol("list[str]"); default = FLAG_BITS, optional = true,
                    doc = "Meaning of each bit, when the flags are stored as a bit field \
                           rather than as booleans"
                ),
                long_name("Visibility flags"),
                versioned(),
            ],
            doc = "True where the corresponding datum is bad. Unsigned element types carry \
                   one flagging reason per bit, described by `flag_bits`."
        ),
        ArraySpec(
            "WEIGHT", shape, Union{Float16, Float32, Float64};
            attrs = [long_name("Visibility weights"), versioned()],
            doc = "Statistical weight of each datum, equal to `1/σ²`. MSv4 weights are \
                   per channel, so this has the same shape as the science array."
        ),
        ArraySpec(
            "EFFECTIVE_INTEGRATION_TIME",
            ((:time, second), (:time, second, :frequency, :polarization)), Float64;
            attrs = quantity("s"), optional = true,
            doc = "Time actually integrated, once flagged data is excluded."
        ),
        ArraySpec(
            "TIME_CENTROID", (:time, second), Float64;
            attrs = time_sampling(), optional = true,
            doc = "Flux-weighted centroid of the integration, which need not be its mid-point."
        ),
        ArraySpec(
            "TIME_CENTROID_EXTRA_PRECISION", (:time, second), Float64;
            attrs = time_sampling(), optional = true,
            doc = "Correction added to `TIME_CENTROID` where one float is not precise enough."
        ),
        ArraySpec(
            "EFFECTIVE_CHANNEL_WIDTH", shape, Float64;
            attrs = frequency_sampling(), optional = true,
            doc = "Bandwidth actually integrated per channel, once flagged data is excluded."
        ),
        ArraySpec(
            "FREQUENCY_CENTROID", (:frequency,), Float64;
            attrs = frequency_sampling(), optional = true,
            doc = "Flux-weighted centroid frequency of each channel."
        ),
    ]
end

const _VIS_SHAPE = (:time, :baseline_id, :frequency, :polarization)

"""
    VISIBILITY_SCHEMA

The main dataset of an interferometric MSv4.

Its science array is `VISIBILITY`, of shape `(time, baseline_id, frequency,
polarization)`. MSv2 stored one row per time-baseline pair; MSv4 splits that
row index into separate `time` and `baseline_id` dimensions, which is what lets
the data be indexed and sliced by label.

A `radiometer` dataset uses this same layout.
"""
const VISIBILITY_SCHEMA = DatasetSpec(
    :visibility, ["visibility", "radiometer"];
    coords = vcat(
        _correlated_coords(:baseline_id),
        ArraySpec[
            ArraySpec(
                "baseline_id", (:baseline_id,), Union{Int32, Int64};
                attrs = [long_name("Baseline ID")],
                doc = "Identifier of each antenna pair."
            ),
            ArraySpec(
                "baseline_antenna1_name", (:baseline_id,), String;
                attrs = [long_name("Antenna name by baseline_id")],
                doc = "Name of the first antenna of each baseline."
            ),
            ArraySpec(
                "baseline_antenna2_name", (:baseline_id,), String;
                attrs = [long_name("Antenna name by baseline_id")],
                doc = "Name of the second antenna of each baseline."
            ),
            ArraySpec(
                "uvw_label", (:uvw_label,), String; optional = true,
                attrs = [long_name("U/V/W label")],
                doc = "Labels of the `UVW` components: `u`, `v`, `w`."
            ),
        ],
    ),
    vars = vcat(
        ArraySpec[
            ArraySpec(
                "VISIBILITY", _VIS_SHAPE, Union{ComplexF32, ComplexF64};
                attrs = [
                    AttrSpec(:units, :str; default = "Jy"),
                    long_name("Visibility values"),
                    versioned(),
                ],
                doc = "Complex visibilities, in janskys."
            ),
            ArraySpec(
                "UVW", (:time, :baseline_id, :uvw_label), Union{Float16, Float32, Float64};
                attrs = uvw_measure(),
                doc = "Baseline vector in the `uvw` frame, in metres: the first antenna's \
                       position minus the second's, with `w` towards the source, `v` in the \
                       plane through the source and the pole, and `u` towards increasing \
                       longitude. The standard states the direction and the difference as \
                       opposites of each other — \"from antenna 2 to antenna 1\" against \
                       `POSITION2 - POSITION1` — and the direction is the clause its own \
                       stores follow, which the test suite measures from antenna positions \
                       rather than assumes. The sign of the visibility phase follows from it."
            ),
        ],
        _correlated_vars(:baseline_id, _VIS_SHAPE),
    ),
    attrs = _correlated_attrs(["visibility", "radiometer"]),
    doc = "Main dataset of an interferometric Measurement Set v4."
)

const _SPEC_SHAPE = (:time, :antenna_name, :frequency, :polarization)

"""
    SPECTRUM_SCHEMA

The main dataset of a single-dish MSv4.

Its science array is `SPECTRUM`, of shape `(time, antenna_name, frequency,
polarization)`. There are no baselines and so no `UVW`.
"""
const SPECTRUM_SCHEMA = DatasetSpec(
    :spectrum, "spectrum";
    coords = vcat(
        _correlated_coords(:antenna_name),
        ArraySpec[
            ArraySpec(
                "antenna_name", (:antenna_name,), String;
                attrs = [long_name("Antenna name")],
                doc = "Name of each antenna contributing a spectrum."
            ),
        ],
    ),
    vars = vcat(
        ArraySpec[
            ArraySpec(
                "SPECTRUM", _SPEC_SHAPE, Union{Float16, Float32, Float64};
                attrs = [
                    AttrSpec(:units, :str; default = "Jy"),
                    long_name("Spectrum values"),
                ],
                doc = "Single-dish spectra, in janskys."
            ),
        ],
        _correlated_vars(:antenna_name, _SPEC_SHAPE),
    ),
    attrs = _correlated_attrs(["spectrum"]),
    doc = "Main dataset of a single-dish Measurement Set v4."
)
