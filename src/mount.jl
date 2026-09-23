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

The reference implementation names nine: `EQUATORIAL`, `ALT-AZ`,
`ALT-AZ+NASMYTH-R`, `ALT-AZ+NASMYTH-L`, `ALT-AZ+BWG-R`, `ALT-AZ+BWG-L`,
`ALT-AZ+ROTATOR`, `X-Y` and `SPACE-HALCA`. The first six rotate the feed by a
fixed combination of the parallactic angle, the elevation and the azimuth, and
[`Mount`](@ref) holds that combination as coefficients; the rest are separate
types, because their rotation is a different function of where the antenna is
pointing, or of how it is driven, and cannot be written as such a combination.

A package with a mount none of these describe — an orbiting station whose
position moves, an array whose feed rotation is tabulated rather than computed
— defines its own subtype and a [`msv4_mount`](@ref) method giving the name to
write. To read it back from a store, it adds a method to
[`XRadio.mount_from_msv4`](@ref) for that name:

```julia
struct MountTabulated{C} <: XRadio.AbstractMount
    offset::C
end
XRadio.msv4_mount(::MountTabulated) = "TABULATED"
XRadio.mount_from_msv4(::Val{:TABULATED}, parallactic, elevation, azimuth, offset) =
    MountTabulated(offset)
```

The name must be one XRadio does not already read, or the store gives back
XRadio's own mount. Storing a mount also stores its [`axis_offset`](@ref), which
reads an `offset` field; a type without one defines `axis_offset` itself.
"""
abstract type AbstractMount end

"""
    Mount(parallactic, elevation, azimuth = 0; offset = (0.0, 0.0, 0.0))

A mount whose feed rotates by `parallactic · χ + elevation · E + azimuth · A`,
where `χ` is the parallactic angle, `E` the elevation and `A` the azimuth.

The six coefficient triples MSv4 names have constructors of their own —
[`MountAltAz`](@ref), [`MountEquatorial`](@ref), [`MountNasmythR`](@ref),
[`MountNasmythL`](@ref), [`MountBWGR`](@ref) and [`MountBWGL`](@ref) — and a
mount those do not anticipate is still expressible as a direct
`Mount(1//2, 0)`. `offset` is the mount's axis offsets; see
[`axis_offset`](@ref).

The named constructors, and [`mounts`](@ref), give `Float64` coefficients, so
an array of the common mounts has one concrete element type.

```jldoctest
julia> XRadio.msv4_mount(MountAltAz())
"ALT-AZ"

julia> XRadio.axis_offset(MountNasmythL((0.0, 0.0, 3.2)))
(0.0, 0.0, 3.2)
```
"""
struct Mount{A, B, C, O} <: AbstractMount
    parallactic::A
    elevation::B
    azimuth::C
    offset::O
end

Mount(parallactic, elevation, azimuth = zero(elevation); offset = ntuple(Returns(0.0), 3)) =
    Mount(parallactic, elevation, azimuth, offset)

_coefficients(m::Mount) = (m.parallactic, m.elevation, m.azimuth)

Base.:(==)(a::Mount, b::Mount) = _coefficients(a) == _coefficients(b) && a.offset == b.offset
Base.hash(m::Mount, h::UInt) = hash((_coefficients(m)..., m.offset), hash(:Mount, h))

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

A mount the writer of a file declined to classify, so nothing about its feed
rotation is known. It is written as `"OTHER"`, and AIPS has no `MNTSTA` code
for it.
"""
struct MountOther{C} <: AbstractMount
    offset::C
end

"""
    MountRotator(offset = (0.0, 0.0, 0.0))

An alt-azimuth mount with a feed rotator, a third axis that turns the feed to
hold it at an angle on the sky; ASKAP's antennas are built this way. The angle
depends on how the rotator is driven, so nothing here computes it. It is
written as `"ALT-AZ+ROTATOR"`.
"""
struct MountRotator{C} <: AbstractMount
    offset::C
end

"""
    MountHALCA(offset = (0.0, 0.0, 0.0))

The mount of HALCA, the space-VLBI antenna of VSOP, whose orientation follows a
model of its own. Nothing here computes it. It is written as `"SPACE-HALCA"`;
[`MountOrbiting`](@ref) is any other antenna in orbit.
"""
struct MountHALCA{C} <: AbstractMount
    offset::C
end

for T in (:MountXY, :MountOrbiting, :MountOther, :MountRotator, :MountHALCA)
    @eval $T() = $T(ntuple(Returns(0.0), 3))
end

"""
    MountAltAz(offset = (0.0, 0.0, 0.0))

An alt-azimuth mount. Its feed rotates by exactly the parallactic angle.
"""
MountAltAz(offset = ntuple(Returns(0.0), 3)) = Mount(1.0, 0.0, 0.0, offset)

"""
    MountEquatorial(offset = (0.0, 0.0, 0.0))

An equatorial mount. Its primary axis points at the celestial pole, so the feed
does not rotate against the sky at all.
"""
MountEquatorial(offset = ntuple(Returns(0.0), 3)) = Mount(0.0, 0.0, 0.0, offset)

"""
    MountNasmythR(offset = (0.0, 0.0, 0.0))

A right Nasmyth mount: alt-azimuth, with the beam taken out along the elevation
bearing by a tertiary mirror that adds the elevation to the feed rotation.
"""
MountNasmythR(offset = ntuple(Returns(0.0), 3)) = Mount(1.0, 1.0, 0.0, offset)

"""
    MountNasmythL(offset = (0.0, 0.0, 0.0))

A left Nasmyth mount, which subtracts the elevation from the feed rotation
where [`MountNasmythR`](@ref) adds it.
"""
MountNasmythL(offset = ntuple(Returns(0.0), 3)) = Mount(1.0, -1.0, 0.0, offset)

"""
    MountBWGR(offset = (0.0, 0.0, 0.0))

A right beam-waveguide mount: alt-azimuth, with the beam brought down to the
pedestal by a guide whose second rotating mirror subtracts the azimuth where
[`MountNasmythR`](@ref) has only the elevation added, so the feed rotates by
`χ + E − A` (Dodson & Rioja, arXiv:2210.13381). A guide that does not bring the
image back to the central axis reverses the sign of `A`, which is a direct
[`Mount`](@ref).
"""
MountBWGR(offset = ntuple(Returns(0.0), 3)) = Mount(1.0, 1.0, -1.0, offset)

"""
    MountBWGL(offset = (0.0, 0.0, 0.0))

A left beam-waveguide mount, whose feed rotates by `χ − E + A`: the signs of
[`MountBWGR`](@ref)'s elevation and azimuth terms reversed.
"""
MountBWGL(offset = ntuple(Returns(0.0), 3)) = Mount(1.0, -1.0, 1.0, offset)

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

A [`Mount`](@ref) whose coefficients are none of the six MSv4 names is an
error rather than a default: the name is all this returns, and there is none
that says what the mount is. [`XRadio.set_mounts!`](@ref) stores such a mount
faithfully, beside its coefficients.
"""
function msv4_mount(m::Mount)
    name = _coefficient_name(_coefficients(m)...)
    name === nothing && throw(ArgumentError(
        "Mount$(_coefficients(m)) is not one of the six coefficient triples " *
            "MSv4 names, so it has no `mount` string to write"
    ))
    return name
end

msv4_mount(::MountXY) = "X-Y"
msv4_mount(::MountOrbiting) = "ORBITING"
msv4_mount(::MountOther) = "OTHER"
msv4_mount(::MountRotator) = "ALT-AZ+ROTATOR"
msv4_mount(::MountHALCA) = "SPACE-HALCA"

# The coefficient triples MSv4 names, and the constructor each name stands for.
const _COEFFICIENT_MOUNTS = (
    ("EQUATORIAL", 0, 0, 0, MountEquatorial),
    ("ALT-AZ", 1, 0, 0, MountAltAz),
    ("ALT-AZ+NASMYTH-R", 1, 1, 0, MountNasmythR),
    ("ALT-AZ+NASMYTH-L", 1, -1, 0, MountNasmythL),
    ("ALT-AZ+BWG-R", 1, 1, -1, MountBWGR),
    ("ALT-AZ+BWG-L", 1, -1, 1, MountBWGL),
)

# The mounts that are types rather than coefficient triples, by the name each
# writes.
const _TYPED_MOUNTS = Dict(
    "X-Y" => MountXY, "ORBITING" => MountOrbiting, "OTHER" => MountOther,
    "ALT-AZ+ROTATOR" => MountRotator, "SPACE-HALCA" => MountHALCA,
)

function _coefficient_name(p, e, a)
    for (name, cp, ce, ca, _) in _COEFFICIENT_MOUNTS
        p == cp && e == ce && a == ca && return name
    end
    return nothing
end

"""
    mount_from_msv4(name, parallactic = NaN, elevation = NaN, azimuth = NaN,
                    offset = (0.0, 0.0, 0.0))
    mount_from_msv4(::Val{name}, parallactic, elevation, azimuth, offset)

The mount a store describes, from its `mount` name and, where the store carries
them, the `MOUNT_PARALLACTIC`, `MOUNT_ELEVATION` and `MOUNT_AZIMUTH`
coefficients and the `ANTENNA_AXIS_OFFSET`.

Finite coefficients describe a [`Mount`](@ref) and are authoritative: a triple
MSv4 names comes back as its own constructor, and any other triple, written
under `"OTHER"`, as a direct `Mount`. `NaN` coefficients — or none — leave the
name to say what the mount is. A name that contradicts the coefficients is an
error, as is a triple stated in part.

A name XRadio does not read itself is passed on as `Val(Symbol(name))`, which
is where a package reads its own mount type back; see
[`AbstractMount`](@ref). A name nothing reads is an error.
"""
function mount_from_msv4(
        name::AbstractString, parallactic = NaN, elevation = NaN, azimuth = NaN,
        offset = ntuple(Returns(0.0), 3)
    )
    triple = (parallactic, elevation, azimuth)
    if !(all(isnan, triple) || !any(isnan, triple))
        throw(ArgumentError(
            "mount `$name` states some of its coefficients and not others " *
                "(parallactic $parallactic, elevation $elevation, azimuth $azimuth)"
        ))
    end
    if !any(isnan, triple)
        for (n, p, e, a, make) in _COEFFICIENT_MOUNTS
            triple == (p, e, a) || continue
            name == n || throw(ArgumentError(
                "mount `$name` states the coefficients $triple, which are those of `$n`"
            ))
            return make(offset)
        end
        name == "OTHER" || throw(ArgumentError(
            "mount `$name` states the coefficients $triple, which MSv4 has no " *
                "name for; a store writes such a mount as `OTHER`"
        ))
        return Mount(map(Float64, triple)..., offset)
    end
    haskey(_TYPED_MOUNTS, name) && return _TYPED_MOUNTS[name](offset)
    for (n, _, _, _, make) in _COEFFICIENT_MOUNTS
        n == name && return make(offset)
    end
    return mount_from_msv4(Val(Symbol(name)), parallactic, elevation, azimuth, offset)
end

mount_from_msv4(::Val{N}, parallactic, elevation, azimuth, offset) where {N} = throw(ArgumentError(
    "XRadio has no mount type for `$N`; it reads " *
        join(repr.(vcat(collect(first.(_COEFFICIENT_MOUNTS)), sort!(collect(keys(_TYPED_MOUNTS))))), ", ") *
        ". A package reads its own with a method " *
        "`XRadio.mount_from_msv4(::Val{Symbol(\"$N\")}, parallactic, elevation, azimuth, offset)`"
))

# The coefficients stored for a mount, `NaN` where the mount is not a
# coefficient triple.
_mount_coefficients(m::Mount) = map(Float64, _coefficients(m))
_mount_coefficients(::AbstractMount) = (NaN, NaN, NaN)

# The name stored beside the coefficients. A triple MSv4 does not name is
# written as `OTHER`, which the coefficients beside it make unambiguous.
_stored_mount_name(m::Mount) = something(_coefficient_name(_coefficients(m)...), "OTHER")
_stored_mount_name(m::AbstractMount) = msv4_mount(m)
