using OffsetArrays: OffsetVector

using XRadio: gain_variable, gain_quantity, msv4_gain_curve, rescale_gain_angle

# A curve of a shape XRadio does not define, which is the whole of what an
# outside package has to write: the gain at the angle, and the angle's name.
struct SquaredCosineGain <: AbstractGainCurve end
XRadio.gain_variable(::SquaredCosineGain) = :zenith_angle
XRadio.gain_quantity(::SquaredCosineGain) = :power
XRadio.gain(::SquaredCosineGain, x::Real) = cos(x)^2

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
