"""
How an antenna's gain varies with where it points.

A gain curve is dimensionless — the gain relative to its own peak — and the
sensitivity in kelvin per jansky that turns it into a flux scale is stated
beside it rather than within it, which is how both FITS-IDI and MSv4 state it.

MSv4 models one shape of curve, a polynomial in elevation or in zenith angle,
and that is all a conforming store can carry. Telescopes have recorded more
than that for as long as there have been gain curves: AIPS Memo 114r admits
tabulated curves and a spherical-harmonic expansion in hour angle and
co-declination as well. The types here model what telescopes measure rather
than what the store can hold, and [`XRadio.msv4_gain_curve`](@ref) is the one
place the difference shows.
"""

"""
    AbstractGainCurve

An antenna's gain as a function of where it points.

The hierarchy is open. A curve none of these describes — a spline, a
measured beam map, a model with a frequency dependence — is a subtype defined
from outside, and everything here works on it.

# Interface

A new curve implements one of

- `gain(curve, x::Real)`, its gain at the single angle it depends on, together
  with `XRadio.gain_variable(curve)` naming that angle; or
- `gain(curve, sky::NamedTuple)` directly, where no single angle determines the
  gain.

together with `XRadio.gain_quantity(curve)`, saying whether it is a power gain
or a voltage gain, and optionally `XRadio.msv4_gain_curve(curve)`, which says
how to write the curve to a conforming store and defaults to `nothing` — the
honest answer for every shape but the polynomial.

Angles are radians throughout, as everywhere else in XRadio. FITS-IDI and
casacore both write these curves against degrees, so a polynomial's
coefficients differ from the ones on disk by a power of `180/π` for each term;
that conversion happens at the file boundary and nowhere else.
"""
abstract type AbstractGainCurve end

"""
    XRadio.GAIN_VARIABLES

The angles a gain curve may be a function of, as
[`XRadio.gain_variable`](@ref) names them. `:zenith_angle` and
`:codeclination` are complements of `:elevation` and `:declination` and are
derived from them rather than asked for separately.
"""
const GAIN_VARIABLES = (:elevation, :zenith_angle, :hour_angle, :declination, :codeclination)

"""
    XRadio.GAIN_QUANTITIES

What a gain curve is a curve *of*. A power gain is the square of the voltage
gain, so a curve is meaningless without saying which it is; MSv4 names both.
"""
const GAIN_QUANTITIES = (:power, :voltage)

function _check_gain_variable(variable)
    variable in GAIN_VARIABLES || throw(ArgumentError(
        "$(repr(variable)) is not an angle a gain curve may depend on; " *
            "XRadio knows $(join(GAIN_VARIABLES, ", "))"
    ))
    return variable
end

function _check_gain_quantity(quantity)
    quantity in GAIN_QUANTITIES || throw(ArgumentError(
        "$(repr(quantity)) is not a gain; a curve is one of " *
            join(GAIN_QUANTITIES, ", ")
    ))
    return quantity
end

"""
    gain(curve, x::Real) -> Real
    gain(curve, sky::NamedTuple) -> Real

The gain of `curve`, either at the angle it depends on or at a position on the
sky.

`sky` gives where the source is, as seen from one antenna at one time, in
radians: any of `elevation`, `hour_angle` and `declination`. A curve takes the
angle it needs and says so when the position does not carry it, so a caller
supplies only what its curve is a function of.

```jldoctest
julia> curve = GainPolynomial([1.0, 0.5]);

julia> gain(curve, 0.4)
1.2

julia> gain(curve, (; elevation = 0.4))
1.2
```

A whole array of curves is evaluated by broadcasting, which needs nothing of
its own:

```julia
gain.(gaincurves(ms)[:curve], deg2rad(45))
```

keeps the antenna and receptor axes, and is `NaN` for an antenna whose gain the
store does not state.
"""
function gain end

gain(curve::AbstractGainCurve, sky::NamedTuple) =
    gain(curve, _sky_angle(sky, gain_variable(curve)))

"""
    XRadio.gain_variable(curve) -> Symbol

The angle `curve` is a function of, one of [`XRadio.GAIN_VARIABLES`](@ref).

A curve that no single angle determines — a spherical-harmonic expansion, a
beam map — defines `gain(curve, ::NamedTuple)` instead of this.
"""
function gain_variable(curve::AbstractGainCurve)
    throw(ArgumentError(
        "$(typeof(curve)) defines neither `XRadio.gain_variable` nor " *
            "`gain(curve, ::NamedTuple)`, so there is no saying what to " *
            "evaluate it at"
    ))
end

"""
    XRadio.gain_quantity(curve) -> Symbol

Whether `curve` is a power gain or a voltage gain, one of
[`XRadio.GAIN_QUANTITIES`](@ref).
"""
function gain_quantity(curve::AbstractGainCurve)
    throw(ArgumentError(
        "$(typeof(curve)) does not define `XRadio.gain_quantity`, so there " *
            "is no saying whether it is a power gain or a voltage gain"
    ))
end

# The angle a curve asks for, from the position a caller gave. Elevation and
# declination are the ones a position states; their complements are arithmetic
# on them and are not asked for twice.
function _sky_angle(sky::NamedTuple, which::Symbol)
    which === :elevation && return _sky_needs(sky, :elevation)
    which === :zenith_angle && return _complement(_sky_needs(sky, :elevation))
    which === :hour_angle && return _sky_needs(sky, :hour_angle)
    which === :declination && return _sky_needs(sky, :declination)
    which === :codeclination && return _complement(_sky_needs(sky, :declination))
    return _check_gain_variable(which)
end

# The angle to the pole from the angle to the horizon, keeping the caller's
# own precision rather than widening it to the precision of `π`.
_complement(x) = oftype(float(x), π / 2) - x

_sky_needs(sky::NamedTuple, name::Symbol) = haskey(sky, name) ? sky[name] :
    throw(ArgumentError(
        "this gain curve is a function of `$name`, which the position does " *
            "not give; it has $(join(keys(sky), ", "))"
    ))

"""
    GainPolynomial(coefficients; variable = :elevation, quantity = :power)
    GainPolynomial(variable, quantity, coefficients)

A gain curve that is a polynomial in one angle, lowest order first.

This is the one shape MSv4 stores, and the shape almost every observatory
publishes: VLBI gain curves are quadratics in elevation or zenith angle, whose
coefficients are fitted once and hold for a season.

```jldoctest
julia> curve = GainPolynomial([0.683, 0.0165, -0.000214]; variable = :zenith_angle);

julia> XRadio.msv4_gain_curve(curve).type
"POWER(ZA)"
```
"""
struct GainPolynomial{C} <: AbstractGainCurve
    variable::Symbol
    quantity::Symbol
    coefficients::C

    function GainPolynomial{C}(variable, quantity, coefficients) where {C}
        _check_gain_variable(variable)
        _check_gain_quantity(quantity)
        isempty(coefficients) && throw(ArgumentError(
            "a gain polynomial has at least a constant term"
        ))
        return new{C}(variable, quantity, coefficients)
    end
end

GainPolynomial(variable, quantity, coefficients) =
    GainPolynomial{typeof(coefficients)}(variable, quantity, coefficients)
GainPolynomial(coefficients; variable::Symbol = :elevation, quantity::Symbol = :power) =
    GainPolynomial(variable, quantity, coefficients)

gain_variable(curve::GainPolynomial) = curve.variable
gain_quantity(curve::GainPolynomial) = curve.quantity

function gain(curve::GainPolynomial, x::Real)
    # Horner over the values rather than the indices, so a polynomial whose
    # coefficients are held in any vector evaluates the same way.
    total = zero(x * zero(eltype(curve.coefficients)))
    for k in Iterators.reverse(curve.coefficients)
        total = total * x + k
    end
    return total
end

"""
    TabulatedGain(angles, gains; variable = :elevation, quantity = :power)
    TabulatedGain(variable, quantity, angles, gains)

A gain curve measured at a list of angles, linearly interpolated between them.

`angles` is in radians and must increase. Asking for the gain outside the range
that was measured is a [`DomainError`](@ref) rather than an extrapolation: a
gain curve says nothing about elevations the telescope was never pointed at.

MSv4 has no field for a curve of this shape, so
[`XRadio.msv4_gain_curve`](@ref) returns `nothing` for it and a conforming
store cannot hold it.
"""
struct TabulatedGain{X, G} <: AbstractGainCurve
    variable::Symbol
    quantity::Symbol
    angles::X
    gains::G

    function TabulatedGain{X, G}(variable, quantity, angles, gains) where {X, G}
        _check_gain_variable(variable)
        _check_gain_quantity(quantity)
        axes(angles) == axes(gains) || throw(DimensionMismatch(
            "a tabulated gain curve has one gain for each angle: " *
                "$(axes(angles)) against $(axes(gains))"
        ))
        length(angles) >= 2 || throw(ArgumentError(
            "a tabulated gain curve is measured at two angles or more"
        ))
        issorted(angles) || throw(ArgumentError(
            "a tabulated gain curve's angles increase; these do not"
        ))
        return new{X, G}(variable, quantity, angles, gains)
    end
end

TabulatedGain(variable, quantity, angles, gains) =
    TabulatedGain{typeof(angles), typeof(gains)}(variable, quantity, angles, gains)
TabulatedGain(angles, gains; variable::Symbol = :elevation, quantity::Symbol = :power) =
    TabulatedGain(variable, quantity, angles, gains)

gain_variable(curve::TabulatedGain) = curve.variable
gain_quantity(curve::TabulatedGain) = curve.quantity

function gain(curve::TabulatedGain, x::Real)
    angles, gains = curve.angles, curve.gains
    first(angles) <= x <= last(angles) || throw(DomainError(
        x,
        "the gain curve was measured from $(first(angles)) to $(last(angles)) " *
            "radians and says nothing outside that"
    ))
    i = searchsortedlast(angles, x)
    i == lastindex(angles) && return gains[i]
    j = nextind(angles, i)
    return gains[i] + (gains[j] - gains[i]) * (x - angles[i]) / (angles[j] - angles[i])
end

"""
    SphericalHarmonicGain(coefficients; quantity = :power)

A gain curve expanded in spherical harmonics of hour angle and co-declination.

`coefficients` are the expansion's ten terms, in increasing degree and order
with each even azimuthal term before its odd one: `A00`, `A10`, `A11E`, `A11O`,
`A20`, `A21E`, `A21O`, `A22E`, `A22O`, `A30`. AIPS Memo 114r fixes that order,
and FITS-IDI files are written in it. The gain at co-declination `θ` and hour
angle `φ` is

    Σ Aₗ₀ Pₗ₀(cos θ) + AₗₘE Pₗₘ(cos θ) cos(m φ) + AₗₘO Pₗₘ(cos θ) sin(m φ)

with the unnormalized associated Legendre functions taken without the
Condon–Shortley phase: `P₁₀ = c`, `P₁₁ = s`, `P₂₀ = (3c² − 1)/2`, `P₂₁ = 3cs`,
`P₂₂ = 3s²`, `P₃₀ = (5c³ − 3c)/2`, for `c = cos θ` and `s = sin θ`. The
convention is written out because AIPS Memo 114r does not state one and a
different normalization would change the gain without changing anything a
reader could check.

The form comes from the Green Bank 140-foot telescope and is rare. MSv4 has no
field for it.
"""
struct SphericalHarmonicGain{C} <: AbstractGainCurve
    quantity::Symbol
    coefficients::C

    function SphericalHarmonicGain{C}(quantity, coefficients) where {C}
        _check_gain_quantity(quantity)
        length(coefficients) == 10 || throw(ArgumentError(
            "a spherical-harmonic gain curve has ten coefficients, one for " *
                "each term of the expansion; this has $(length(coefficients))"
        ))
        return new{C}(quantity, coefficients)
    end
end

SphericalHarmonicGain(quantity, coefficients) =
    SphericalHarmonicGain{typeof(coefficients)}(quantity, coefficients)
SphericalHarmonicGain(coefficients; quantity::Symbol = :power) =
    SphericalHarmonicGain(quantity, coefficients)

gain_quantity(curve::SphericalHarmonicGain) = curve.quantity

function gain(curve::SphericalHarmonicGain, sky::NamedTuple)
    θ = _sky_angle(sky, :codeclination)
    φ = _sky_angle(sky, :hour_angle)
    c, s = cos(θ), sin(θ)
    a00, a10, a11e, a11o, a20, a21e, a21o, a22e, a22o, a30 = curve.coefficients
    p10, p11 = c, s
    p20, p21, p22 = (3c^2 - 1) / 2, 3c * s, 3s^2
    p30 = (5c^3 - 3c) / 2
    return a00 +
        a10 * p10 + p11 * (a11e * cos(φ) + a11o * sin(φ)) +
        a20 * p20 + p21 * (a21e * cos(φ) + a21o * sin(φ)) +
        p22 * (a22e * cos(2φ) + a22o * sin(2φ)) +
        a30 * p30
end

"""
    AtDeclination(curve, declination)

A gain curve that holds at one declination, which a curve against hour angle
does.

The gain of an antenna tracking a source depends on where the source is, and
hour angle alone does not say: the same hour angle at a different declination
is a different elevation. FITS-IDI records the declination beside such a curve
and this carries it, so that a curve read from a file states the conditions it
was measured under rather than losing them.

Evaluating this is evaluating the curve it wraps, at a position that carries
`declination`: a caller who names a declination of their own keeps it, and one
who does not gets this one. Nothing here checks that the position asked about
is at `declination` — an antenna's gain is not undefined half a degree away,
and where the tolerance lies is the caller's to judge.
"""
struct AtDeclination{C <: AbstractGainCurve, T} <: AbstractGainCurve
    curve::C
    declination::T
end

gain(w::AtDeclination, x::Real) = gain(w.curve, x)
gain(w::AtDeclination, sky::NamedTuple) =
    gain(w.curve, merge((; declination = w.declination), sky))
gain_variable(w::AtDeclination) = gain_variable(w.curve)
gain_quantity(w::AtDeclination) = gain_quantity(w.curve)

# The four curves MSv4 names, by what they are a gain of and what they are a
# function of.
const _MSV4_GAIN_CURVE_TYPES = Dict(
    (:power, :elevation) => "POWER(EL)",
    (:power, :zenith_angle) => "POWER(ZA)",
    (:voltage, :elevation) => "VOLTAGE(EL)",
    (:voltage, :zenith_angle) => "VOLTAGE(ZA)",
)

"""
    XRadio.msv4_gain_curve(curve) -> (; type, coefficients) or nothing

How to write `curve` to a conforming store, or `nothing` where MSv4 cannot
express it.

`type` is the name MSv4's `gain_curve_type` takes and `coefficients` are the
polynomial's, **against the angle in degrees**, which is the convention
casacore and AIPS write and MSv4 inherits without stating. In memory the same
polynomial is against radians, so each term is scaled by `π/180` to its own
order on the way out.

Returning `nothing` is the honest answer for every shape but the polynomial,
and for a polynomial in an angle MSv4 does not name. It is what makes a curve
XRadio can model but the standard cannot hold visible rather than silent: a
converter drops such a curve from the store and says so, instead of writing
coefficients that do not mean what the field says they mean.
"""
msv4_gain_curve(::AbstractGainCurve) = nothing

function msv4_gain_curve(curve::GainPolynomial)
    type = get(_MSV4_GAIN_CURVE_TYPES, (curve.quantity, curve.variable), nothing)
    type === nothing && return nothing
    return (; type, coefficients = to_degree_angle(curve.coefficients))
end

"""
    XRadio.to_radian_angle(coefficients) -> Vector
    XRadio.to_degree_angle(coefficients) -> Vector

The same polynomial against the other angle unit: FITS-IDI and casacore write a
gain curve against degrees, XRadio evaluates against radians, and the term of
order `k` differs by the ratio to the `k`th power.

The ratio takes the coefficients' own precision. As a bare `Float64` it would
widen a single-precision curve at every file boundary, so a file or a store
holding `Float32` could not be read back as one.
"""
to_radian_angle(coefficients) =
    rescale_gain_angle(coefficients, _angle_factor(coefficients, 180 / π))

@doc (@doc to_radian_angle)
to_degree_angle(coefficients) =
    rescale_gain_angle(coefficients, _angle_factor(coefficients, π / 180))

_angle_factor(coefficients, factor) = oftype(float(first(coefficients)), factor)

msv4_gain_curve(w::AtDeclination) = msv4_gain_curve(w.curve)

# The inverse of `_MSV4_GAIN_CURVE_TYPES`, so a store's `gain_curve_type` says
# what the coefficients beside it are a gain of and a function of.
const _GAIN_CURVE_FROM_MSV4 = Dict(v => k for (k, v) in _MSV4_GAIN_CURVE_TYPES)

"""
    XRadio.gain_curve_from_msv4(type, coefficients) -> GainPolynomial

The curve a conforming store's `gain_curve_type` and `GAIN_CURVE` row describe
— the inverse of [`XRadio.msv4_gain_curve`](@ref).

`coefficients` are against the angle in **degrees**, which is how MSv4 writes
them, and come back against radians in the precision they were held in.

Trailing zeros are dropped. MSv4 pads every antenna's polynomial out to the
longest one in its dataset and records no term count, so a trailing zero could
be padding or a real coefficient and nothing in the store tells them apart. A
zero top-order coefficient contributes nothing to the gain, so the shortest
polynomial giving the same gain is taken; a curve that is zero throughout keeps
its constant term.

```jldoctest
julia> curve = XRadio.gain_curve_from_msv4("POWER(ZA)", [0.99, 0.0, 0.0]);

julia> XRadio.gain_variable(curve), length(curve.coefficients)
(:zenith_angle, 1)
```
"""
function gain_curve_from_msv4(type::AbstractString, coefficients)
    known = get(_GAIN_CURVE_FROM_MSV4, String(type), nothing)
    known === nothing && throw(ArgumentError(
        "`$type` is not a gain curve MSv4 names; it holds one of " *
            join(sort(collect(keys(_GAIN_CURVE_FROM_MSV4))), ", ")
    ))
    quantity, variable = known
    return GainPolynomial(
        variable, quantity, to_radian_angle(_drop_trailing_zeros(coefficients))
    )
end

# Everything from the last non-zero coefficient down, and never less than the
# constant term — a gain polynomial has at least one.
function _drop_trailing_zeros(coefficients)
    stop = lastindex(coefficients)
    while stop > firstindex(coefficients) && iszero(coefficients[stop])
        stop -= 1
    end
    return coefficients[firstindex(coefficients):stop]
end

"""
    XRadio.rescale_gain_angle(coefficients, factor) -> Vector

The coefficients of the same polynomial in an angle scaled by `factor`.

A polynomial in radians becomes one in degrees with `factor = π/180` and back
with `factor = 180/π`, since the term of order `k` picks up `factor^k`. Both
directions are needed at a file boundary, because gain curves are written
against degrees and evaluated here against radians.
"""
rescale_gain_angle(coefficients, factor) =
    [k * factor^(i - 1) for (i, k) in enumerate(coefficients)]
