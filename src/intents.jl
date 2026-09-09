"""
Scan intents, as the ASDM enumerations define them.

MSv4 states a scan intent as `"INTENT#SUBINTENT"` — one member of ASDM's
`ScanIntent` enumeration, optionally followed by one of its `SubscanIntent`
enumeration. Both vocabularies are closed, so [`Intent`](@ref) holds two
enumeration members rather than a string, and a value outside them is an error
rather than a name that quietly matches nothing.

`ScanIntent` and `SubscanIntent` each define `TEST` and `UNSPECIFIED`, so the
members are scoped to a module apiece, as ASDM scopes them to a namespace.

Spectral window intents are a different vocabulary — "continuum", "spectral
line" — which neither MSv4 nor ASDM enumerates, so
[`spectral_window_intents`](@ref) stays strings.

The member lists are transcribed from ASDM's generated `CScanIntent.h` and
`CSubscanIntent.h`. Take any addition from those files rather than by hand.
"""

"""
    ScanIntents

What a scan was observed for. ASDM's `ScanIntent` enumeration.

The enumeration type itself is `ScanIntents.T`.
"""
@enumx ScanIntents begin
    """
    Amplitude calibration scan.
    """
    CALIBRATE_AMPLI

    """
    Atmosphere calibration scan.
    """
    CALIBRATE_ATMOSPHERE

    """
    Bandpass calibration scan.
    """
    CALIBRATE_BANDPASS

    """
    Delay calibration scan.
    """
    CALIBRATE_DELAY

    """
    Flux measurement scan.
    """
    CALIBRATE_FLUX

    """
    Focus calibration scan. Z coordinate to be derived.
    """
    CALIBRATE_FOCUS

    """
    Focus calibration scan; X focus coordinate to be derived.
    """
    CALIBRATE_FOCUS_X

    """
    Focus calibration scan; Y focus coordinate to be derived.
    """
    CALIBRATE_FOCUS_Y

    """
    Phase calibration scan.
    """
    CALIBRATE_PHASE

    """
    Pointing calibration scan.
    """
    CALIBRATE_POINTING

    """
    Polarization calibration scan.
    """
    CALIBRATE_POLARIZATION

    """
    Measure relative gains of sidebands.
    """
    CALIBRATE_SIDEBAND_RATIO

    """
    Data from the water vapor radiometers (and correlation data) are used to derive their calibration parameters.
    """
    CALIBRATE_WVR

    """
    Skydip calibration scan.
    """
    DO_SKYDIP

    """
    Holography calibration scan.
    """
    MAP_ANTENNA_SURFACE

    """
    Data on a celestial calibration source are used to derive a map of the primary beam.
    """
    MAP_PRIMARY_BEAM

    """
    Target source scan.
    """
    OBSERVE_TARGET

    CALIBRATE_POL_LEAKAGE

    CALIBRATE_POL_ANGLE

    """
    Used for development.
    """
    TEST

    """
    Unspecified scan intent.
    """
    UNSPECIFIED

    """
    Requested by EVLA.
    """
    CALIBRATE_ANTENNA_POSITION

    """
    Requested by EVLA.
    """
    CALIBRATE_ANTENNA_PHASE

    """
    Requested by EVLA.
    """
    MEASURE_RFI

    """
    Requested by EVLA.
    """
    CALIBRATE_ANTENNA_POINTING_MODEL

    """
    Requested by EVLA.
    """
    SYSTEM_CONFIGURATION

    """
    Calculate and apply phasing solutions. Applicable at ALMA.
    """
    CALIBRATE_APPPHASE_ACTIVE

    """
    Apply previously obtained phasing solutions. Applicable at ALMA.
    """
    CALIBRATE_APPPHASE_PASSIVE

    OBSERVE_CHECK_SOURCE

    """
    Enable a gain differential target type.
    """
    CALIBRATE_DIFFGAIN
end

"""
    SubscanIntents

What a subscan within a scan was for. ASDM's `SubscanIntent` enumeration.

The enumeration type itself is `SubscanIntents.T`.
"""
@enumx SubscanIntents begin
    """
    On-source measurement.
    """
    ON_SOURCE

    """
    Off-source measurement.
    """
    OFF_SOURCE

    """
    Pointing measurement, some antennas are on -ource, some off-source.
    """
    MIXED

    """
    Reference measurement (used for boresight in holography).
    """
    REFERENCE

    """
    Antennas are scanning.
    """
    SCANNING

    """
    Hot load measurement.
    """
    HOT

    """
    Ambient load measurement.
    """
    AMBIENT

    """
    Signal sideband measurement.
    """
    SIGNAL

    """
    Image sideband measurement.
    """
    IMAGE

    """
    Reserved for development.
    """
    TEST

    """
    Unspecified.
    """
    UNSPECIFIED
end

"""
    Intent(intent::ScanIntents.T, subintent = nothing)
    Intent(s::AbstractString)

One scan intent: what a scan was observed for, and what a subscan within it
was for.

```julia
Intent(ScanIntents.OBSERVE_TARGET, SubscanIntents.ON_SOURCE)
Intent("OBSERVE_TARGET#ON_SOURCE")     # the same, parsed
```

Parsing takes MSv4's `"INTENT#SUBINTENT"` form, and the subintent may be
omitted — converted stores write a bare `"UNSPECIFIED"`. An omitted subintent
is `nothing`, which is not the same as `SubscanIntents.UNSPECIFIED`: the first
says the store names no subintent, the second says it names that one. So
`string(intent)` reproduces exactly what the store held.

Both halves come from closed ASDM enumerations, so a value outside them throws
rather than becoming a name that matches nothing.
"""
struct Intent
    intent::ScanIntents.T
    subintent::Union{SubscanIntents.T, Nothing}
end

Intent(intent::ScanIntents.T) = Intent(intent, nothing)

function Intent(s::AbstractString)
    parts = split(s, '#')
    length(parts) <= 2 || throw(
        ArgumentError(
            "`$s` is not a scan intent; MSv4 states one as \"INTENT#SUBINTENT\", " *
                "or as an intent alone"
        )
    )
    intent = _enum_member(ScanIntents.T, parts[1], "scan intent")
    length(parts) == 1 && return Intent(intent)
    return Intent(intent, _enum_member(SubscanIntents.T, parts[2], "subscan intent"))
end

Intent(i::Intent) = i

Base.string(i::Intent) =
    i.subintent === nothing ? string(i.intent) : string(i.intent, '#', i.subintent)
Base.print(io::IO, i::Intent) = print(io, string(i))
Base.show(io::IO, i::Intent) = print(io, "Intent(\"", string(i), "\")")

# The ASDM enumerations spell their members exactly as MSv4 stores them, so a
# member is found by name rather than through a table.
function _enum_member(::Type{T}, name::AbstractString, what::AbstractString) where {T}
    for member in instances(T)
        String(Symbol(member)) == name && return member
    end
    throw(
        ArgumentError(
            "`$name` is not an ASDM $what; the $(length(instances(T))) it defines are " *
                join((String(Symbol(m)) for m in instances(T)), ", ")
        )
    )
end
