# A mount an outside package might add: its feed rotation is measured rather
# than computed, so no coefficient pair describes it, and naming it is the whole
# of what XRadio needs.
struct TabulatedMount <: XRadio.AbstractMount end
XRadio.msv4_mount(::TabulatedMount) = "ALT-AZ"

@testset "antenna mounts" begin
    @testset "the coefficient mounts name themselves" begin
        @test XRadio.msv4_mount(MountAltAz()) == "ALT-AZ"
        @test XRadio.msv4_mount(MountEquatorial()) == "EQUATORIAL"
        @test XRadio.msv4_mount(MountNasmythR()) == "ALT-AZ+NASMYTH-R"
        @test XRadio.msv4_mount(MountNasmythL()) == "ALT-AZ+NASMYTH-L"

        @test MountAltAz() == Mount(1, 0)
        @test (MountNasmythR().elevation, MountNasmythL().elevation) == (1, -1)
    end

    # The feed rotation of these is not a combination of parallactic angle and
    # elevation, so they are types rather than coefficient pairs.
    @testset "mounts outside that family are their own types" begin
        @test XRadio.msv4_mount(MountXY()) == "X-Y"
        @test XRadio.msv4_mount(MountOrbiting()) == "OTHER"
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
end
