"""
How an antenna carries its feed.

The mount fixes how the feed rotates against the sky as the antenna tracks, and
so how parallactic angle enters the antenna's polarization response. MSv4
stores it as a bare name in the `antenna` dataset's `mount`; the types here
carry the rotation law itself, which the name only stands for.
"""

"""
    AbstractMount

An antenna mount.

MSv4 names five: `ALT-AZ`, `EQUATORIAL`, `X-Y`, `ALT-AZ+NASMYTH-R` and
`ALT-AZ+NASMYTH-L`. Four of those rotate the feed by a fixed combination of the
parallactic angle and the elevation, and [`Mount`](@ref) holds that combination
as coefficients; the rest are separate types, because their rotation is a
different function of where the antenna is pointing and cannot be written as
such a combination.

A package with a mount none of these describe — an orbiting station whose
position moves, an array whose feed rotation is tabulated rather than computed
— defines its own subtype and a [`msv4_mount`](@ref) method giving the name to
write. That is the whole of the interface; [`axis_offset`](@ref) is defined for
mounts that have one.
"""
abstract type AbstractMount end

# The axis offset of a mount whose two axes intersect, which is most of them.
const _NO_OFFSET = (0.0, 0.0, 0.0)

"""
    Mount(parallactic, elevation, offset = (0.0, 0.0, 0.0))

A mount whose feed rotates by `parallactic · χ + elevation · E`, where `χ` is
the parallactic angle and `E` the elevation.

The four coefficient pairs MSv4 names have constructors of their own —
[`MountAltAz`](@ref), [`MountEquatorial`](@ref), [`MountNasmythR`](@ref) and
[`MountNasmythL`](@ref) — and a mount those do not anticipate is still
expressible as a direct `Mount(1//2, 0)`. `offset` is the mount's axis offsets;
see [`axis_offset`](@ref).

```jldoctest
julia> XRadio.msv4_mount(MountAltAz())
"ALT-AZ"

julia> XRadio.axis_offset(MountNasmythL((0.0, 0.0, 3.2)))
(0.0, 0.0, 3.2)
```
"""
struct Mount{A, B, C} <: AbstractMount
    parallactic::A
    elevation::B
    offset::C
end

Mount(parallactic, elevation) = Mount(parallactic, elevation, _NO_OFFSET)

"""
    MountXY(offset = (0.0, 0.0, 0.0))

A mount whose primary axis is horizontal rather than vertical.

Its keyhole — the direction the primary axis points, where tracking a source
demands an unbounded slew rate — sits on the horizon instead of at the zenith,
which is why satellite-tracking antennas are built this way. The feed rotation
follows the position angle of that fixed horizon point rather than of the
zenith, so it is not the parallactic angle scaled or offset by anything, and no
[`Mount`](@ref) describes it.
"""
struct MountXY{C} <: AbstractMount
    offset::C
end

"""
    MountOrbiting(offset = (0.0, 0.0, 0.0))

The mount of an antenna in orbit, whose feed rotation depends on the
spacecraft's attitude and position rather than on a station's horizon.

Nothing here computes that rotation. The type records which antennas are in
orbit so a package that does can dispatch on them.
"""
struct MountOrbiting{C} <: AbstractMount
    offset::C
end

"""
    MountOther(offset = (0.0, 0.0, 0.0))

A mount the writer of a file declined to classify.

FITS-IDI's `MNTSTA = 6` says only that the mount is none of the ones AIPS
names, so nothing about the feed rotation is known.
"""
struct MountOther{C} <: AbstractMount
    offset::C
end

for T in (:MountXY, :MountOrbiting, :MountOther)
    @eval $T() = $T(_NO_OFFSET)
end

"""
    MountAltAz(offset = (0.0, 0.0, 0.0))

An alt-azimuth mount. Its feed rotates by exactly the parallactic angle.
"""
MountAltAz(offset = _NO_OFFSET) = Mount(1, 0, offset)

"""
    MountEquatorial(offset = (0.0, 0.0, 0.0))

An equatorial mount. Its primary axis points at the celestial pole, so the feed
does not rotate against the sky at all.
"""
MountEquatorial(offset = _NO_OFFSET) = Mount(0, 0, offset)

"""
    MountNasmythR(offset = (0.0, 0.0, 0.0))

A right Nasmyth mount: alt-azimuth, with the beam taken out along the elevation
bearing by a tertiary mirror that adds the elevation to the feed rotation.
"""
MountNasmythR(offset = _NO_OFFSET) = Mount(1, 1, offset)

"""
    MountNasmythL(offset = (0.0, 0.0, 0.0))

A left Nasmyth mount, which subtracts the elevation from the feed rotation
where [`MountNasmythR`](@ref) adds it.
"""
MountNasmythL(offset = _NO_OFFSET) = Mount(1, -1, offset)

"""
    axis_offset(m::AbstractMount)

Displacement between the mount's two rotation axes, in metres, along each axis
of the station's own frame. It is zero on a mount whose axes intersect, and is
what FITS-IDI's `STAXOF` records.

MSv4 has no field for it, so it survives only in memory.
"""
axis_offset(m::AbstractMount) = m.offset

"""
    msv4_mount(m::AbstractMount) -> String

The name MSv4's `mount` field takes for `m`.

A [`Mount`](@ref) whose coefficients are none of the four MSv4 names is an
error rather than a default: the store records a name, and there is none to
record.
"""
function msv4_mount(m::Mount)
    p, e = m.parallactic, m.elevation
    p == 0 && e == 0 && return "EQUATORIAL"
    p == 1 && e == 0 && return "ALT-AZ"
    p == 1 && e == 1 && return "ALT-AZ+NASMYTH-R"
    p == 1 && e == -1 && return "ALT-AZ+NASMYTH-L"
    throw(ArgumentError(
        "Mount($p, $e) is not one of the four coefficient pairs MSv4 names, " *
            "so it has no `mount` string to write"
    ))
end

msv4_mount(::MountXY) = "X-Y"
msv4_mount(::MountOrbiting) = "OTHER"
msv4_mount(::MountOther) = "OTHER"
