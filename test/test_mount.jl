# A mount an outside package might add: its feed rotation is measured rather
# than computed, so no coefficient pair describes it, and naming it is the whole
# of what XRadio needs.
struct TabulatedMount <: XRadio.AbstractMount
    offset::NTuple{3, Float64}
end
TabulatedMount() = TabulatedMount((0.0, 0.0, 0.0))
XRadio.msv4_mount(::TabulatedMount) = "ALT-AZ"

# The same, under a name of its own, so a store gives it back.
struct MeasuredMount{C} <: XRadio.AbstractMount
    offset::C
end
XRadio.msv4_mount(::MeasuredMount) = "TEST-MEASURED"
XRadio.mount_from_msv4(::Val{Symbol("TEST-MEASURED")}, parallactic, elevation, offset) =
    MeasuredMount(offset)

@testset "antenna mounts" begin
    @testset "the coefficient mounts name themselves" begin
        @test XRadio.msv4_mount(MountAltAz()) == "ALT-AZ"
        @test XRadio.msv4_mount(MountEquatorial()) == "EQUATORIAL"
        @test XRadio.msv4_mount(MountNasmythR()) == "ALT-AZ+NASMYTH-R"
        @test XRadio.msv4_mount(MountNasmythL()) == "ALT-AZ+NASMYTH-L"

        @test MountAltAz() == Mount(1, 0)
        @test (MountNasmythR().elevation, MountNasmythL().elevation) == (1, -1)
        @test Mount(1 // 2, 0) == Mount(0.5, 0.0)
        @test hash(Mount(1 // 2, 0)) == hash(Mount(0.5, 0.0))
    end

    # The feed rotation of these is not a combination of parallactic angle and
    # elevation, so they are types rather than coefficient pairs.
    @testset "mounts outside that family are their own types" begin
        @test XRadio.msv4_mount(MountXY()) == "X-Y"
        @test XRadio.msv4_mount(MountOrbiting()) == "ORBITING"
        @test XRadio.msv4_mount(MountOther()) == "OTHER"
        @test !any(m -> m isa Mount, (MountXY(), MountOrbiting(), MountOther()))
        @test all(m -> m isa XRadio.AbstractMount, (MountXY(), MountOrbiting(), MountOther()))
    end

    @testset "a coefficient pair MSv4 does not name has no name to write" begin
        @test_throws "not one of the four coefficient pairs MSv4 names" XRadio.msv4_mount(Mount(1 // 2, 0))
    end

    @testset "the axis offset survives in memory" begin
        @test XRadio.axis_offset(MountAltAz((1.0, 2.0, 3.0))) == (1.0, 2.0, 3.0)
        @test XRadio.axis_offset(MountXY((0.0, 0.0, 4.5))) == (0.0, 0.0, 4.5)
        @test XRadio.axis_offset(MountEquatorial()) == (0.0, 0.0, 0.0)
    end

    @testset "an outside package needs only to name its mount" begin
        @test XRadio.msv4_mount(TabulatedMount()) == "ALT-AZ"
    end

    @testset "a name reads back as its mount" begin
        @test XRadio.mount_from_msv4("ALT-AZ") === MountAltAz()
        @test XRadio.mount_from_msv4("ALT-AZ+NASMYTH-L") === MountNasmythL()
        @test XRadio.mount_from_msv4("X-Y", NaN, NaN, (0.0, 0.0, 2.0)) === MountXY((0.0, 0.0, 2.0))
        @test XRadio.mount_from_msv4("ORBITING") === MountOrbiting()
        @test XRadio.mount_from_msv4("OTHER") === MountOther()
        @test_throws "XRadio has no mount type for `ALT-AZ+BWG-R`" XRadio.mount_from_msv4("ALT-AZ+BWG-R")
    end

    @testset "coefficients are authoritative where stated" begin
        @test XRadio.mount_from_msv4("ALT-AZ", 1.0, 0.0) === MountAltAz()
        @test XRadio.mount_from_msv4("OTHER", 0.5, 0.0) === Mount(0.5, 0.0)
        @test_throws "which are those of `ALT-AZ`" XRadio.mount_from_msv4("EQUATORIAL", 1.0, 0.0)
        @test_throws "a store writes such a mount as `OTHER`" XRadio.mount_from_msv4("ALT-AZ", 0.5, 0.0)
        @test_throws "states one coefficient and not the other" XRadio.mount_from_msv4("ALT-AZ", 1.0, NaN)
    end

    @testset "every mount round-trips through an antenna dataset" begin
        xds = synthetic_antenna(["AA", "BB", "CC", "DD"])
        stored = AbstractMount[
            MountNasmythR((0.0, 0.0, 1.5)), Mount(1 // 2, 0, (0.5, 0.0, 0.0)),
            MountOrbiting(), MountXY((0.0, 2.0, 0.0)),
        ]
        XRadio.set_mounts!(xds, stored)
        @test isempty(check(xds))
        @test collect(xds[:mount]) == ["ALT-AZ+NASMYTH-R", "OTHER", "ORBITING", "X-Y"]
        @test collect(mounts(xds)) == stored
        @test lookup(mounts(xds), AntennaName) == ["AA", "BB", "CC", "DD"]
        @test eltype(mounts(xds)) == Union{
            Mount{Float64, Float64, NTuple{3, Float64}},
            MountOrbiting{NTuple{3, Float64}}, MountXY{NTuple{3, Float64}},
        }
        @test_throws DimensionMismatch XRadio.set_mounts!(xds, stored[1:2])
    end

    @testset "a store without the coefficients reads the name alone" begin
        xds = synthetic_antenna(["AA", "BB"])
        @test collect(mounts(xds)) == fill(MountAltAz(), 2)
        @test isconcretetype(eltype(mounts(xds)))
    end

    @testset "per-antenna accessors on a Measurement Set" begin
        ms = synthetic_set(1)[1]
        ant = branches(ms)[:antenna]
        @test antenna_positions(ms) == ant[:antenna_position]
        @test polarization_types(ms) == ant[:polarization_type]
        @test collect(mounts(ms)) == fill(MountAltAz(), length(antennas(ms)))

        # MSv4 makes the receptor angle optional, and an angle nobody stated is
        # unknown rather than zero.
        angles = receptor_angles(ms)
        @test all(isnan, angles)
        @test lookup(angles, AntennaName) == antennas(ms)
        @test lookup(angles, ReceptorLabel) == ["X", "Y"]
        ant[:antenna_receptor_angle] = DimArray(
            zeros(size(angles)), (dims(angles, ReceptorLabel), dims(angles, AntennaName))
        )
        @test all(iszero, receptor_angles(ms))
    end

    @testset "a Measurement Set without antennas describes none" begin
        ms = synthetic_set(1)[1]
        delete!(branches(ms), :antenna)
        @test isempty(mounts(ms))
        @test size(antenna_positions(ms)) == (3, 0)
        @test isempty(receptor_angles(ms))
        @test isempty(polarization_types(ms))
    end

    @testset "an outside mount with a name of its own reads back" begin
        xds = synthetic_antenna(["AA", "BB"])
        stored = [MountAltAz(), MeasuredMount((0.0, 0.0, 2.0))]
        XRadio.set_mounts!(xds, stored)
        @test collect(xds[:mount]) == ["ALT-AZ", "TEST-MEASURED"]
        @test collect(mounts(xds)) == stored
        @test eltype(mounts(xds)) ==
            Union{typeof(MountAltAz()), MeasuredMount{NTuple{3, Float64}}}

        # A standard name stays XRadio's, whatever a package defines.
        XRadio.set_mounts!(xds, [TabulatedMount(), MountAltAz()])
        @test collect(mounts(xds)) == [MountAltAz(), MountAltAz()]

        @test_throws "A package reads its own with a method" XRadio.mount_from_msv4("NO-SUCH-MOUNT")
    end
end
