using OffsetArrays: OffsetVector

using XRadio: gain_variable, gain_quantity, msv4_gain_curve, rescale_gain_angle
using XRadio: gain_curve_from_msv4

# A curve of a shape XRadio does not define, which is the whole of what an
# outside package has to write: the gain at the angle, and the angle's name.
struct SquaredCosineGain <: AbstractGainCurve end
XRadio.gain_variable(::SquaredCosineGain) = :zenith_angle
XRadio.gain_quantity(::SquaredCosineGain) = :power
XRadio.gain(::SquaredCosineGain, x::Real) = cos(x)^2

# A curve that implements none of the interface, which is how a curve defined
# outside XRadio arrives before it says anything about itself.
struct BareGain <: AbstractGainCurve end

@testset "gain curves" begin
    @testset "a polynomial is evaluated where it says it is" begin
        # A VLBA Q-band curve: 0.683 + 0.0165 za - 0.000214 za^2 in degrees,
        # which is what the file holds and what this holds in radians.
        degrees = [0.683, 0.0165, -0.000214]
        curve = GainPolynomial(
            rescale_gain_angle(degrees, 180 / π); variable = :zenith_angle
        )
        @test gain_variable(curve) == :zenith_angle
        @test gain_quantity(curve) == :power

        # The same number by either route, which is the point of rescaling.
        @test gain(curve, deg2rad(45)) ≈ evalpoly(45, degrees)
        @test gain(curve, (; elevation = deg2rad(45))) ≈ evalpoly(45, degrees)
        # A curve against elevation is a different curve, not a spelling.
        elevation = GainPolynomial(curve.coefficients)
        @test gain(elevation, (; elevation = deg2rad(45))) ≈ evalpoly(45, degrees)
        @test gain(elevation, (; elevation = deg2rad(20))) ≈ evalpoly(20, degrees)
        @test gain(curve, (; elevation = deg2rad(20))) ≈ evalpoly(70, degrees)

        # Coefficients are held lowest order first, and a constant curve is a
        # curve.
        @test gain(GainPolynomial([2.0]), 1.7) == 2.0
        @test gain(GainPolynomial([0.0, 1.0]), 1.7) == 1.7
    end

    @testset "degrees on disk, radians in memory" begin
        degrees = [0.683, 0.0165, -0.000214]
        radians = rescale_gain_angle(degrees, 180 / π)
        @test radians[1] == degrees[1]
        @test radians[2] ≈ degrees[2] * 180 / π
        @test radians[3] ≈ degrees[3] * (180 / π)^2
        # The round trip is what a file boundary does in both directions.
        @test rescale_gain_angle(radians, π / 180) ≈ degrees
    end

    @testset "what MSv4 can hold, and what it cannot" begin
        coefficients = [1.0, 0.1]
        for (variable, quantity, type) in (
                (:elevation, :power, "POWER(EL)"),
                (:zenith_angle, :power, "POWER(ZA)"),
                (:elevation, :voltage, "VOLTAGE(EL)"),
                (:zenith_angle, :voltage, "VOLTAGE(ZA)"),
            )
            written = msv4_gain_curve(GainPolynomial(coefficients; variable, quantity))
            @test written.type == type
            @test written.coefficients ≈ rescale_gain_angle(coefficients, π / 180)
        end

        # MSv4 names no curve against hour angle or declination, so there is
        # nothing to write rather than something to approximate.
        @test msv4_gain_curve(
            GainPolynomial(coefficients; variable = :hour_angle)
        ) === nothing
        @test msv4_gain_curve(TabulatedGain([0.1, 0.2], [0.9, 1.0])) === nothing
        @test msv4_gain_curve(SphericalHarmonicGain(zeros(10))) === nothing
        # A curve from outside says nothing to MSv4 unless it says otherwise.
        @test msv4_gain_curve(SquaredCosineGain()) === nothing
    end

    @testset "what a store holds comes back as the curve it was" begin
        coefficients = [1.0, 0.1]
        for (variable, quantity, type) in (
                (:elevation, :power, "POWER(EL)"),
                (:zenith_angle, :power, "POWER(ZA)"),
                (:elevation, :voltage, "VOLTAGE(EL)"),
                (:zenith_angle, :voltage, "VOLTAGE(ZA)"),
            )
            curve = GainPolynomial(coefficients; variable, quantity)
            back = gain_curve_from_msv4(msv4_gain_curve(curve)...)
            @test gain_variable(back) === variable
            @test gain_quantity(back) === quantity
            @test back.coefficients ≈ coefficients
        end

        # A store pads every antenna's polynomial out to the longest one in its
        # dataset and records no term count. A zero top-order coefficient
        # cannot change a gain, so the shortest polynomial saying the same
        # thing is what comes back.
        padded = gain_curve_from_msv4("POWER(EL)", [0.99, 0.01, 0.0, 0.0])
        @test length(padded.coefficients) == 2
        @test gain(padded, 0.5) ≈ gain(
            gain_curve_from_msv4("POWER(EL)", [0.99, 0.01]), 0.5
        )

        # A gain polynomial has at least a constant term, so one that is zero
        # throughout keeps it rather than becoming empty.
        @test gain_curve_from_msv4("POWER(EL)", zeros(4)).coefficients == [0.0]

        # The coefficients come back in the precision they were held in: the
        # degrees-to-radians rescaling takes their type rather than widening
        # every store to `Float64`.
        single = GainPolynomial(Float32[0.9, 0.01])
        written = msv4_gain_curve(single)
        @test eltype(written.coefficients) === Float32
        @test eltype(gain_curve_from_msv4(written...).coefficients) === Float32
        @test gain_curve_from_msv4("POWER(EL)", Float32[NaN]).coefficients[1] === NaN32

        @test_throws "is not a gain curve MSv4 names" gain_curve_from_msv4(
            "POWER(HA)", [1.0]
        )
    end

    @testset "a whole layer of curves is evaluated by broadcasting" begin
        # A store states no curve for some antennas, and each of those holds a
        # curve whose coefficients are `NaN` rather than a marker standing in
        # for its absence. The layer is one concrete type throughout, so
        # evaluating it gives an array of numbers with gaps.
        over = (AntennaName(["AA", "BB"]), ReceptorLabel(["pol_0", "pol_1"]))
        layer = DimArray(
            [
                GainPolynomial([1.0, 0.1])  GainPolynomial([NaN])
                GainPolynomial([2.0, 0.2])  GainPolynomial([3.0, 0.3])
            ],
            over,
        )
        @test eltype(layer) === GainPolynomial{Vector{Float64}}

        got = gain.(layer, 0.5)
        @test dims(got) == over
        @test got[At("AA"), At("pol_0")] ≈ 1.05
        @test isnan(got[At("AA"), At("pol_1")])
        @test got[At("BB"), At("pol_1")] ≈ 3.15

        # Nothing in the layer needs a type of its own, so the result is
        # concrete whatever precision the angle is asked in.
        @test eltype(got) === Float64
        @test eltype(gain.(layer, 0.5f0)) === Float64
        @test eltype(gain.(layer, 1)) === Float64
        @test eltype(gain.(layer, Ref((; elevation = 0.5)))) === Float64

        # The precision follows the coefficients rather than being `Float64`,
        # and a gap in single precision is `NaN32`.
        small = DimArray(
            [GainPolynomial(Float32[1, 0.1]), GainPolynomial(Float32[NaN])],
            over[1],
        )
        @test eltype(gain.(small, 0.5f0)) === Float32
        @test gain.(small, 0.5f0)[At("BB")] === NaN32
    end

    @testset "a store that names no angle is not a curve" begin
        # MSv4 makes `gain_curve_type` optional, but coefficients that name no
        # angle could be against elevation or against zenith angle, and the two
        # disagree everywhere but the horizon.
        xds = XRadio.subdataset(
            :gain_curve;
            antenna_name = ["AA", "BB"], receptor_label = ["pol_0", "pol_1"],
            poly_term = Int64[0, 1],
            gain_curve = reshape(Float64[1, 2, 3, 4, 5, 6, 7, 8], 2, 2, 2),
            gain_curve_sensitivity = fill(0.1, 2, 2),
            gain_curve_interval = fill(60.0, 2),
        )
        @test_throws "carries no `gain_curve_type`" gaincurves(xds)
    end

    @testset "a tabulated curve interpolates and does not extrapolate" begin
        angles = deg2rad.([10.0, 20.0, 40.0])
        curve = TabulatedGain(angles, [0.90, 0.95, 1.00])

        @test gain(curve, angles[1]) == 0.90
        @test gain(curve, angles[end]) == 1.00
        @test gain(curve, deg2rad(15)) ≈ 0.925
        # Unequal spacing is interpolated on its own interval.
        @test gain(curve, deg2rad(30)) ≈ 0.975
        @test gain(curve, (; elevation = deg2rad(15))) ≈ 0.925

        # A gain curve says nothing about angles the telescope was not
        # pointed at, and saying so is better than extrapolating a fit.
        @test_throws DomainError gain(curve, deg2rad(5))
        @test_throws DomainError gain(curve, deg2rad(41))

        @test_throws "one gain for each angle" TabulatedGain(angles, [0.9, 1.0])
        @test_throws "angles increase" TabulatedGain(reverse(angles), [0.9, 0.95, 1.0])
        @test_throws "two angles or more" TabulatedGain([0.1], [0.9])
    end

    @testset "a spherical-harmonic curve" begin
        sky = (; hour_angle = 0.3, declination = 0.9)

        # The constant term is the gain when nothing else is present.
        @test gain(SphericalHarmonicGain([2.5; zeros(9)]), sky) == 2.5

        # A00 and A10 alone give A00 + A10 cos(co-declination), and the
        # co-declination is the complement of the declination the caller gave.
        c = [1.0, 0.5, zeros(8)...]
        @test gain(SphericalHarmonicGain(c), sky) ≈ 1.0 + 0.5 * cos(π / 2 - 0.9)

        # The odd azimuthal terms turn with hour angle where the even ones do
        # not, which is what tells them apart.
        even = [0.0, 0.0, 1.0, zeros(7)...]
        odd = [0.0, 0.0, 0.0, 1.0, zeros(6)...]
        s = sin(π / 2 - 0.9)
        @test gain(SphericalHarmonicGain(even), sky) ≈ s * cos(0.3)
        @test gain(SphericalHarmonicGain(odd), sky) ≈ s * sin(0.3)

        @test_throws "has ten coefficients" SphericalHarmonicGain(zeros(9))
        # It is a function of two angles, so one is not enough.
        @test_throws "which the position does not give" gain(
            SphericalHarmonicGain(zeros(10)), (; hour_angle = 0.3)
        )
    end

    @testset "a curve measured at one declination carries it" begin
        inner = GainPolynomial([1.0, 0.2]; variable = :hour_angle)
        curve = AtDeclination(inner, deg2rad(30))

        @test curve.declination ≈ deg2rad(30)
        @test gain_variable(curve) == :hour_angle
        @test gain_quantity(curve) == :power
        @test gain(curve, 0.4) == gain(inner, 0.4)
        @test gain(curve, (; hour_angle = 0.4)) == gain(inner, 0.4)

        # The declination it carries is the one the wrapped curve sees, so
        # wrapping a curve that needs a declination makes it evaluable at a
        # position that gives only the hour angle.
        harmonic = SphericalHarmonicGain([1.0, 0.5, zeros(8)...])
        at = AtDeclination(harmonic, deg2rad(30))
        @test gain(at, (; hour_angle = 0.3)) ==
            gain(harmonic, (; hour_angle = 0.3, declination = deg2rad(30)))
        # A caller who says where the source is keeps their own declination.
        @test gain(at, (; hour_angle = 0.3, declination = 0.9)) ==
            gain(harmonic, (; hour_angle = 0.3, declination = 0.9))
        # Wrapping does not make a curve writable, and does not make one
        # unwritable either.
        @test msv4_gain_curve(curve) === nothing
        @test msv4_gain_curve(AtDeclination(GainPolynomial([1.0]), 0.5)).type ==
            "POWER(EL)"
    end

    @testset "a curve defined outside XRadio" begin
        curve = SquaredCosineGain()
        @test gain(curve, 0.0) == 1.0
        # The generic path resolves the angle the curve asks for from the
        # position, which is all an outside curve has to say.
        @test gain(curve, (; elevation = π / 2)) == 1.0
        @test gain(curve, (; elevation = 0.0)) ≈ 0.0 atol = 1.0e-15
    end

    @testset "a curve that says neither what it is nor what of" begin
        # Both halves of the interface an outside curve implements fail the
        # same way when it implements neither: a sentence naming what is
        # missing.
        @test_throws "defines neither `XRadio.gain_variable`" gain_variable(BareGain())
        @test_throws "does not define `XRadio.gain_quantity`" gain_quantity(BareGain())
    end

    @testset "what a curve refuses" begin
        @test_throws "not an angle a gain curve may depend on" GainPolynomial(
            [1.0]; variable = :azimuth
        )
        @test_throws "is not a gain" GainPolynomial([1.0]; quantity = :temperature)
        @test_throws "at least a constant term" GainPolynomial(Float64[])

        # A position that does not carry the angle names what it does carry.
        @test_throws "is a function of `elevation`" gain(
            GainPolynomial([1.0]), (; hour_angle = 0.2)
        )
    end

    @testset "coefficients in any vector" begin
        # Horner runs over the coefficients' values, so a polynomial whose
        # terms are held in an array indexed from zero evaluates the same.
        coefficients = [1.0, 0.5, 0.25]
        offset = OffsetVector(coefficients, 0:2)
        @test gain(GainPolynomial(offset), 0.4) == gain(GainPolynomial(coefficients), 0.4)
        @test msv4_gain_curve(GainPolynomial(offset)).coefficients ≈
            msv4_gain_curve(GainPolynomial(coefficients)).coefficients

        angles = OffsetVector(deg2rad.([10.0, 20.0, 40.0]), 0:2)
        gains = OffsetVector([0.90, 0.95, 1.00], 0:2)
        @test gain(TabulatedGain(angles, gains), deg2rad(15)) ≈ 0.925
    end
end
