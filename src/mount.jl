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
write. To read it back from a store, it adds a method to
[`XRadio.mount_from_msv4`](@ref) for that name:

```julia
struct MountBWGR{C} <: XRadio.AbstractMount
    offset::C
end
XRadio.msv4_mount(::MountBWGR) = "ALT-AZ+BWG-R"
XRadio.mount_from_msv4(::Val{Symbol("ALT-AZ+BWG-R")}, parallactic, elevation, offset) =
    MountBWGR(offset)
```

The name must be one XRadio does not already read, or the store gives back
XRadio's own mount. Storing a mount also stores its [`axis_offset`](@ref), which
reads an `offset` field; a type without one defines `axis_offset` itself.
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

The named constructors, and [`mounts`](@ref), give `Float64` coefficients, so
an array of the common mounts has one concrete element type.

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

Base.:(==)(a::Mount, b::Mount) =
    a.parallactic == b.parallactic && a.elevation == b.elevation && a.offset == b.offset
Base.hash(m::Mount, h::UInt) = hash((m.parallactic, m.elevation, m.offset), hash(:Mount, h))

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
orbit so a package that does can dispatch on them. It is written as
`"ORBITING"`, which MSv4 does not reserve but accepts, as it accepts any name.
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
MountAltAz(offset = _NO_OFFSET) = Mount(1.0, 0.0, offset)

"""
    MountEquatorial(offset = (0.0, 0.0, 0.0))

An equatorial mount. Its primary axis points at the celestial pole, so the feed
does not rotate against the sky at all.
"""
MountEquatorial(offset = _NO_OFFSET) = Mount(0.0, 0.0, offset)

"""
    MountNasmythR(offset = (0.0, 0.0, 0.0))

A right Nasmyth mount: alt-azimuth, with the beam taken out along the elevation
bearing by a tertiary mirror that adds the elevation to the feed rotation.
"""
MountNasmythR(offset = _NO_OFFSET) = Mount(1.0, 1.0, offset)

"""
    MountNasmythL(offset = (0.0, 0.0, 0.0))

A left Nasmyth mount, which subtracts the elevation from the feed rotation
where [`MountNasmythR`](@ref) adds it.
"""
MountNasmythL(offset = _NO_OFFSET) = Mount(1.0, -1.0, offset)

"""
    axis_offset(m::AbstractMount)

Displacement between the mount's two rotation axes, in metres, along each axis
of the station's own frame. It is zero on a mount whose axes intersect, and is
what FITS-IDI's `STAXOF` records.

MSv4 has no field for it; [`XRadio.set_mounts!`](@ref) stores it as the
`ANTENNA_AXIS_OFFSET` extension variable.
"""
axis_offset(m::AbstractMount) = m.offset

"""
    msv4_mount(m::AbstractMount) -> String

The name MSv4's `mount` field takes for `m`.

A [`Mount`](@ref) whose coefficients are none of the four MSv4 names is an
error rather than a default: the name is all this returns, and there is none
that says what the mount is. [`XRadio.set_mounts!`](@ref) stores such a mount
faithfully, beside its coefficients.
"""
function msv4_mount(m::Mount)
    name = _coefficient_name(m.parallactic, m.elevation)
    name === nothing && throw(ArgumentError(
        "Mount($(m.parallactic), $(m.elevation)) is not one of the four " *
            "coefficient pairs MSv4 names, so it has no `mount` string to write"
    ))
    return name
end

msv4_mount(::MountXY) = "X-Y"
msv4_mount(::MountOrbiting) = "ORBITING"
msv4_mount(::MountOther) = "OTHER"

# The coefficient pairs MSv4 names, and the constructor each name stands for.
const _COEFFICIENT_MOUNTS = (
    ("EQUATORIAL", 0, 0, MountEquatorial),
    ("ALT-AZ", 1, 0, MountAltAz),
    ("ALT-AZ+NASMYTH-R", 1, 1, MountNasmythR),
    ("ALT-AZ+NASMYTH-L", 1, -1, MountNasmythL),
)

# The mounts that are types rather than coefficient pairs, by the name each
# writes.
const _TYPED_MOUNTS = Dict("X-Y" => MountXY, "ORBITING" => MountOrbiting, "OTHER" => MountOther)

function _coefficient_name(p, e)
    for (name, cp, ce, _) in _COEFFICIENT_MOUNTS
        p == cp && e == ce && return name
    end
    return nothing
end

"""
    mount_from_msv4(name, parallactic = NaN, elevation = NaN, offset = (0.0, 0.0, 0.0))
    mount_from_msv4(::Val{name}, parallactic, elevation, offset)

The mount a store describes, from its `mount` name and, where the store carries
them, the `MOUNT_PARALLACTIC` and `MOUNT_ELEVATION` coefficients and the
`ANTENNA_AXIS_OFFSET`.

Finite coefficients describe a [`Mount`](@ref) and are authoritative: a pair
MSv4 names comes back as its own constructor, and any other pair, written under
`"OTHER"`, as a direct `Mount`. `NaN` coefficients — or none — leave the name to
say what the mount is. A name that contradicts the coefficients is an error.

A name XRadio does not read itself is passed on as `Val(Symbol(name))`, which
is where a package reads its own mount type back; see
[`AbstractMount`](@ref). A name nothing reads is an error.
"""
function mount_from_msv4(
        name::AbstractString, parallactic = NaN, elevation = NaN, offset = _NO_OFFSET
    )
    stated = !isnan(parallactic) && !isnan(elevation)
    if isnan(parallactic) != isnan(elevation)
        throw(ArgumentError(
            "mount `$name` states one coefficient and not the other " *
                "(parallactic $parallactic, elevation $elevation)"
        ))
    end
    if stated
        for (n, p, e, make) in _COEFFICIENT_MOUNTS
            (parallactic == p && elevation == e) || continue
            name == n || throw(ArgumentError(
                "mount `$name` states the coefficients ($parallactic, $elevation), " *
                    "which are those of `$n`"
            ))
            return make(offset)
        end
        name == "OTHER" || throw(ArgumentError(
            "mount `$name` states the coefficients ($parallactic, $elevation), " *
                "which MSv4 has no name for; a store writes such a mount as `OTHER`"
        ))
        return Mount(Float64(parallactic), Float64(elevation), offset)
    end
    haskey(_TYPED_MOUNTS, name) && return _TYPED_MOUNTS[name](offset)
    for (n, _, _, make) in _COEFFICIENT_MOUNTS
        n == name && return make(offset)
    end
    return mount_from_msv4(Val(Symbol(name)), parallactic, elevation, offset)
end

mount_from_msv4(::Val{N}, parallactic, elevation, offset) where {N} = throw(ArgumentError(
    "XRadio has no mount type for `$N`; it reads " *
        join(repr.(vcat(collect(first.(_COEFFICIENT_MOUNTS)), sort!(collect(keys(_TYPED_MOUNTS))))), ", ") *
        ". A package reads its own with a method " *
        "`XRadio.mount_from_msv4(::Val{Symbol(\"$N\")}, parallactic, elevation, offset)`"
))

# The coefficients stored for a mount, `NaN` where the mount is not a
# coefficient pair.
_mount_coefficients(m::Mount) = (Float64(m.parallactic), Float64(m.elevation))
_mount_coefficients(::AbstractMount) = (NaN, NaN)

# The name stored beside the coefficients. A pair MSv4 does not name is written
# as `OTHER`, which the coefficients beside it make unambiguous.
_stored_mount_name(m::Mount) = something(_coefficient_name(m.parallactic, m.elevation), "OTHER")
_stored_mount_name(m::AbstractMount) = msv4_mount(m)
