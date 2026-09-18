"""
MSv4's metadata dictionaries.

Where MSv2 kept small, mostly-scalar tables — `OBSERVATION`, `PROCESSOR` — MSv4
folds them into plain dictionaries on the main dataset's attributes. The
specifications here say what those dictionaries must contain.
"""

"""
    OBSERVATION_INFO

Who observed, under which project, and when the data may become public.
Corresponds to MSv2's `OBSERVATION` table.
"""
const OBSERVATION_INFO = DictSpec(
    :observation_info,
    [
        AttrSpec(:observer, Symbol("list[str]"); doc = "Names of the observers"),
        AttrSpec(:release_date, :str; doc = "ISO 8601 date on which the data may become public"),
        AttrSpec(:project_UID, :str; doc = "Project code, or the ASDM ExecBlock `projectUID`"),
        AttrSpec(:execution_block_UID, :str; optional = true, doc = "Archive UID of the execution block"),
        AttrSpec(:session_reference_UID, :str; optional = true, doc = "Observing session reference"),
        AttrSpec(:observing_log, :str; optional = true, doc = "Observing log as supplied by the telescope"),
        AttrSpec(:scheduling_block_UID, :str; optional = true, doc = "Archive UID of the scheduling block"),
    ];
    doc = "Observation-level metadata, from MSv2's OBSERVATION table."
)

"""
    PROCESSOR_INFO

What produced the correlated data. Corresponds to MSv2's `PROCESSOR` table.
"""
const PROCESSOR_INFO = DictSpec(
    :processor_info,
    [
        AttrSpec(
            :type, :str;
            doc = "Processor type; reserved values are CORRELATOR (interferometric), \
                   SPECTROMETER (single dish), and RADIOMETER (generic detector)"
        ),
        AttrSpec(:sub_type, :str; doc = "Processor sub-type, such as `GBT` or `JIVE`"),
    ];
    doc = "Correlator or spectrometer metadata, from MSv2's PROCESSOR table."
)

"""
    CREATOR

The software that wrote the Measurement Set.
"""
const CREATOR = DictSpec(
    :creator,
    [
        AttrSpec(:software_name, :str; doc = "Software that wrote the Measurement Set"),
        AttrSpec(:version, :str; doc = "Version of that software"),
    ];
    doc = "Provenance of the Measurement Set."
)

"""
    DATA_GROUP

One entry of the `data_groups` dictionary: a set of correlated data, flags,
weights and baseline coordinates that belong together.

A dataset may hold several versions of its science arrays — `VISIBILITY`
alongside `VISIBILITY_CORRECTED`, say — and a data group is what records which
flags and weights go with which. Variables may be shared between groups. The
group named `base` is mandatory.
"""
const DATA_GROUP = DictSpec(
    :data_group,
    [
        AttrSpec(:correlated_data, :str; doc = "Name of the correlated data variable, such as `VISIBILITY`"),
        AttrSpec(:flag, :str; doc = "Name of the group's flag variable"),
        AttrSpec(:weight, :str; doc = "Name of the group's weight variable"),
        AttrSpec(:uvw, :str; optional = true, doc = "Name of the group's UVW variable; absent for single-dish data"),
        AttrSpec(:field_and_source, :str; doc = "Name of the group's field and source sub-dataset"),
        AttrSpec(:description, :str; doc = "What distinguishes this group"),
        AttrSpec(:date, :str; doc = "ISO 8601 creation time of the group"),
    ];
    doc = "One correlated-data / flag / weight / uvw grouping."
)

"""
    DATA_GROUPS

The `data_groups` attribute of a correlated dataset: one [`DATA_GROUP`](@ref)
per grouping, keyed by the group's name. `base` is the only name the standard
fixes; a dataset may add whatever others it needs, and their contents are its
own business.
"""
const DATA_GROUPS = DictSpec(
    :data_groups,
    [AttrSpec(:base, :dict; nested = DATA_GROUP, doc = "The grouping every dataset carries")];
    doc = "The correlated-data groupings a dataset holds."
)

"""
    INFO_DICTS

Every metadata dictionary the standard defines, by name.
"""
const INFO_DICTS = Dict{Symbol, DictSpec}(
    d.name => d for d in (OBSERVATION_INFO, PROCESSOR_INFO, CREATOR, DATA_GROUP, DATA_GROUPS)
)

"""
    MSV4_SCHEMA_VERSION

The version of the MSv4 schema this package implements.
"""
const MSV4_SCHEMA_VERSION = v"4.0.0"
