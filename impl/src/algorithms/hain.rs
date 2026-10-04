/// The algorithm implemented here is based on:
///     "Fast, precise flattening of cubic Bézier path and offset curves"
///     http://sibgrapi.sid.inpe.br/col/sid.inpe.br/banon/2004/08.13.18.12/doc/BezierOffsetRendering.pdf
///
/// The basic premise is that for a small t the third order term in the
/// equation of a cubic bezier curve is insignificantly small. This can
/// then be approximated by a quadratic equation for which the maximum
/// difference from a linear approximation can be much more easily determined.
use crate::{CubicBezierSegment, LineSegment};

const EPSILON: f32 = 1e-4;

#[inline(always)]
fn in_range(t: f32) -> bool {
    t >= 0.0 && t < 1.0
}

pub fn flatten_cubic<F: FnMut(&LineSegment)>(
    bezier: &CubicBezierSegment,
    tolerance: f32,
    call_back: &mut F,
) {
    let mut inflection_1 = -1.0;
    let mut inflection_2 = 2.0;
    let inflection_count =
        find_cubic_bezier_inflection_points(&bezier, &mut inflection_1, &mut inflection_2);

    if inflection_count == 0 {
        flatten_cubic_no_inflection(bezier, 0.0, 1.0, tolerance, call_back);
        return;
    }

    // From now on we know that we have one or two inflection points.

    let mut t1_min = inflection_1;
    let mut t1_max = inflection_1;
    inflection_approximation_range(bezier, inflection_1, tolerance, &mut t1_min, &mut t1_max);

    let mut t2_min = inflection_2;
    let mut t2_max = inflection_2;
    if inflection_count == 2 {
        inflection_approximation_range(bezier, inflection_2, tolerance, &mut t2_min, &mut t2_max);
    }

    if inflection_count == 1 && t1_min < 0.0 && t1_max > 1.0 {
        // The first inflection range covers the entire curve.
        call_back(&bezier.baseline());
        return;
    }

    let mut from = bezier.from;

    if t1_min > 0.0 {
        from = flatten_cubic_no_inflection(bezier, 0.0, t1_min, tolerance, call_back);
    }

    if in_range(t1_max) && (inflection_count == 1 || t2_min > t1_max) {
        // There is no second inflection point or the second inflection point's range is in the first's.
        // Add a line to the end of the first approximation range;
        let after = bezier.after_split(t1_max);
        if from != after.from {
            call_back(&LineSegment {
                from,
                to: after.from,
            });
        }
        from = after.from;

        if inflection_count == 1 || t2_min > 1.0 {
            flatten_cubic_no_inflection(bezier, t1_max, 1.0, tolerance, call_back);
            from = after.to;
        }
    } else if inflection_count == 1 && t1_max >= 1.0 {
        // The first inflection range extends to the end of the curve:
        // approximate everything up to the endpoint with a line segment.
        call_back(&LineSegment {
            from,
            to: bezier.to,
        });
        return;
    } else if inflection_count == 2 && t2_min > 1.0 {
        call_back(&LineSegment {
            from,
            to: bezier.to,
        });
        return;
    }

    if inflection_count == 2 && t2_min < 1.0 && t2_max > 0.0 {
        if t2_min > 0.0 && t2_min < t1_max {
            // Skip t2_min since it is in the first approximation range.
            // Clamp to 1.0: an approximation range may extend past the end of
            // the curve, and sample() extrapolates beyond it.
            let to = bezier.sample(t1_max.min(1.0));
            call_back(&LineSegment { from, to });
            from = to;
        } else if t2_min > 0.0 && t1_max > 0.0 {
            from = flatten_cubic_no_inflection(bezier, t1_max, t2_min, tolerance, call_back);
        } else if t2_min > 0.0 {
            from = flatten_cubic_no_inflection(bezier, 0.0, t2_min, tolerance, call_back);
        }

        if t2_max < 1.0 {
            let after = bezier.after_split(t2_max);
            if from != after.from {
                call_back(&LineSegment {
                    from,
                    to: after.from,
                });
            }
            flatten_cubic_no_inflection(bezier, t2_max, 1.0, tolerance, call_back);
        } else if from != bezier.to {
            call_back(&LineSegment {
                from,
                to: bezier.to,
            });
        }
    }
}

fn flatten_cubic_no_inflection<F: FnMut(&LineSegment)>(
    original: &CubicBezierSegment,
    lo: f32,
    hi: f32,
    tolerance: f32,
    call_back: &mut F,
) -> crate::Point {
    // Flatten the [lo, hi] portion of `original` into chords.
    //
    // Each chord endpoint comes from a fresh split of the original curve at
    // the step boundary instead of from chaining the remainder sub-curve:
    // chaining adds one rounding error per split, and on large curves
    // flattened with hundreds of chords the accumulated drift pushes the
    // polyline away from the curve by more than the tolerance. A fresh split
    // of the original costs the same single de Casteljau pass per chord.
    let end = if hi >= 1.0 {
        original.to
    } else if lo <= 0.0 {
        original.before_split(hi).to
    } else {
        original.split_range(lo..hi).to
    };

    let hi = hi as f64;
    let mut ta = lo as f64;
    let mut sub = if lo <= 0.0 {
        *original
    } else {
        original.after_split(lo)
    };

    loop {
        let step = no_inflection_flattening_step(&sub, tolerance);
        let tb = ta + step as f64 * (1.0 - ta);

        if tb >= hi {
            call_back(&LineSegment {
                from: sub.from,
                to: end,
            });
            return end;
        }

        let next = original.after_split(tb as f32);
        call_back(&LineSegment {
            from: sub.from,
            to: next.from,
        });
        ta = tb;
        sub = next;
    }
}

fn no_inflection_flattening_step(bezier: &CubicBezierSegment, tolerance: f32) -> f32 {
    let v1 = bezier.ctrl1 - bezier.from;
    let v2 = bezier.ctrl2 - bezier.from;
    let v3 = bezier.to - bezier.from;

    let a = (v1.x * v1.x + v1.y * v1.y).sqrt();
    if a == 0.0 {
        // The curve starts with zero velocity: the frame used below is
        // undefined. Only degenerate sub-curves end up here.
        return 1.0;
    }

    // In the frame at the start of the curve (origin at bezier.from, x axis
    // along v1), the curve is exactly
    //     x(u) = 3 a u + (3 x2 - 6 a) u^2 + (3 a - 3 x2 + x3) u^3
    //     s(u) = 3 y2 u^2 + (y3 - 3 y2) u^3
    // with x2 = dot(v2, v1) / a and x3 = dot(v3, v1) / a the longitudinal
    // components, and y2 = cross(v1, v2) / a and y3 = cross(v1, v3) / a the
    // transverse ones.
    //
    // The paper's step size is derived from the transverse part alone,
    // assuming x(u) ~ 3 a u. That under-estimates the chord distance wherever
    // x(u) departs from 3 a u: on curves whose start tangent is not aligned
    // with the direction of travel (for example near a cusp) the true distance
    // can be twice the approximation. Instead, use the exact distance from the
    // curve point at u to the chord over [0, t]:
    //     d(u) = u * (t - u) * |A + B * u| / |B(t)|
    // with A = 3 a t (s2 + s3 t), B = s3 x(t) - q s(t),
    //      s2 = 3 y2, s3 = y3 - 3 y2, p = 3 x2 - 6 a, q = 3 a - 3 x2 + x3,
    //      x(t) = ((q t + p) t + 3 a) t and s(t) = (s2 + s3 t) t^2.
    // This is exact: the numerator is a cubic in u that vanishes at u = 0 and
    // u = t, so dividing it by u (t - u) leaves a linear function, and the
    // maximum of u (t - u) |A + B u| over [0, t] is found by solving a
    // quadratic (Rolle's theorem guarantees a critical point inside [0, t]).
    //
    // Take t where the maximum of d(u) over [0, t] meets the tolerance budget
    // (0.9 * tolerance minus an f32 rounding allowance, see below), using the
    // fixed-point iteration t <- t * (budget / d_max(t))^(1/3). The exponent
    // matches the fastest growth of d_max (between quadratic and cubic in t),
    // so the iteration converges from either side, and a final verification
    // loop shrinks any step whose certified error still exceeds the budget.
    let x2 = v2.dot(v1) / a;
    let x3 = v3.dot(v1) / a;
    let y2 = v1.cross(v2) / a;
    let y3 = v1.cross(v3) / a;
    let s2 = 3.0 * y2;
    let s3 = y3 - 3.0 * y2;
    if s2 == 0.0 && s3 == 0.0 {
        // The remaining curve lies on its start tangent: its baseline chord
        // is an exact approximation.
        return 1.0;
    }
    let p = 3.0 * x2 - 6.0 * a;
    let q = 3.0 * a - 3.0 * x2 + x3;

    // Reserve part of the budget for the f32 rounding of the split itself:
    // chord endpoints and sub-curve control points are rounded to f32, so the
    // polyline deviates from the original curve by a couple of ulp of the
    // coordinate magnitude even when the certified chord distance is zero.
    let max_coord = bezier
        .from
        .x
        .abs()
        .max(bezier.from.y.abs())
        .max(bezier.ctrl1.x.abs())
        .max(bezier.ctrl1.y.abs())
        .max(bezier.ctrl2.x.abs())
        .max(bezier.ctrl2.y.abs())
        .max(bezier.to.x.abs())
        .max(bezier.to.y.abs());
    let budget = (0.9 * tolerance - 2.0 * f32::EPSILON * max_coord).max(0.1 * tolerance);

    // Maximum of d(u) over u in [0, t].
    //
    // Solving the critical-point equation directly,
    //     u = (tB - A +/- sqrt(A^2 + A t B + t^2 B^2)) / (3B),
    // loses the root near u = t/2 to cancellation whenever the cubic term is
    // negligible (|B t| << |A|): both subtracted terms are ~|A| while their
    // difference is ~tB/2, so f32 rounding can move the root anywhere in
    // [0, t]. A root that lands near 0 or t makes the maximum look ~0 and the
    // fixed point below then overshoots the safe step size by orders of
    // magnitude. Instead substitute v = u/t and z = v - 1/2, which turns the
    // critical points into the roots of
    //     3 b z^2 + (2 A + b) z - b/4 = 0
    // with A = 3 a t (s2 + s3 t) and b = t B. Their product is the constant
    // -1/12: compute the well-conditioned root (numerator terms of the same
    // sign) and obtain the other by Vieta, which needs no subtraction at all.
    // Whatever rounding is left in b only slides the root along the flat top
    // of the maximum, where |f| is insensitive to it first-order. The value at
    // v = 1/2 (the exact maximizer when b = 0) is always included as a floor.
    let d_max = |t: f32| -> f32 {
        let w = s2 + s3 * t;
        let xt = ((q * t + p) * t + 3.0 * a) * t;
        let st = w * t * t;
        let bt2 = xt * xt + st * st;
        if bt2 == 0.0 {
            // Degenerate chord: the curve comes back to its start point.
            return f32::INFINITY;
        }
        let big_a = 3.0 * a * t * w;
        let b = t * (s3 * xt - q * st);
        let mut m: f32 = 0.25 * (big_a + 0.5 * b).abs();
        if b != 0.0 {
            let x = 2.0 * big_a + b;
            let s = (x * x + 3.0 * b * b).sqrt();
            let z1 = if x >= 0.0 { -(x + s) } else { s - x } / (6.0 * b);
            let z2 = -1.0 / (12.0 * z1);
            for v in [z1 + 0.5, z2 + 0.5] {
                if v > 0.0 && v < 1.0 {
                    m = m.max((v * (1.0 - v) * (big_a + b * v)).abs());
                }
            }
        }
        t * t * m / bt2.sqrt()
    };

    // Start from the single-term solutions of the paper's transverse bound
    //     (|s2| / 4) t^2 + (|s3| / 2) t^3 = budget,
    // clamped to 1. Both are >= the root of the combined bound, and the exact
    // d_max differs from that bound only by the longitudinal correction, so
    // the start is within a small factor of the target on either side.
    let mut t = (4.0 * budget / s2.abs())
        .sqrt()
        .min((2.0 * budget / s3.abs()).cbrt())
        .min(1.0);
    let mut d = d_max(t);
    for _ in 0..4 {
        if t >= 0.995 && d <= budget {
            return 1.0;
        }
        let factor = if d > 0.0 && d.is_finite() {
            (budget / d).cbrt()
        } else {
            // Degenerate chord: back off and let the caller make progress.
            0.5
        };
        if (factor - 1.0).abs() < 0.02 {
            break;
        }
        t *= factor;
        d = d_max(t);
    }
    // The fixed point assumes d_max grows at most like t^3 between
    // evaluations, which a chord whose far end curls back toward the start
    // can violate. Never return a step whose certified error exceeds the
    // budget: shrink until it does (d_max vanishes like t^2, so this
    // terminates; the floor of 0.5 avoids over-shrinking when d is barely
    // over budget).
    let mut guard = 0;
    while d > budget && guard < 16 {
        let factor = (budget / d).cbrt();
        t *= if factor < 0.5 { 0.5 } else { factor };
        d = d_max(t);
        guard += 1;
    }

    // TODO: We start having floating point precision issues if this constant
    // is closer to 1.0 with a small enough tolerance threshold.
    if t >= 0.995 || t <= 0.0 {
        return 1.0;
    }

    return t;
}

// Find the inflection points of a cubic bezier curve.
pub(crate) fn find_cubic_bezier_inflection_points(
    bezier: &CubicBezierSegment,
    t1: &mut f32,
    t2: &mut f32,
) -> u32 {
    // Find inflection points.
    // See www.faculty.idc.ac.il/arik/quality/appendixa.html for an explanation
    // of this approach.
    let pa = bezier.ctrl1 - bezier.from;
    let pb = bezier.ctrl2.to_vector() - (bezier.ctrl1.to_vector() * 2.0) + bezier.from.to_vector();
    let pc = bezier.to.to_vector() - (bezier.ctrl2.to_vector() * 3.0)
        + (bezier.ctrl1.to_vector() * 3.0)
        - bezier.from.to_vector();

    let a = pb.cross(pc);
    let b = pa.cross(pc);
    let c = pa.cross(pb);

    if f32::abs(a) < EPSILON {
        // Not a quadratic equation.
        if f32::abs(b) < EPSILON {
            // Instead of a linear acceleration change we have a constant
            // acceleration change. This means the equation has no solution
            // and there are no inflection points, unless the constant is 0.
            // In that case the curve is a straight line, essentially that means
            // the easiest way to deal with is is by saying there's an inflection
            // point at t == 0. The inflection point approximation range found will
            // automatically extend into infinity.
            if f32::abs(c) < EPSILON {
                *t1 = 0.0;
                return 1;
            }
        } else {
            let t = -c / b;
            if in_range(t) {
                *t1 = t;
                return 1;
            }
        }

        return 0;
    }

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < 0.0 {
        return 0;
    }

    if discriminant < EPSILON * b * b {
        // The two roots are merged into a double root at the midpoint.
        // The threshold is relative to b (the root separation is sqrt(discriminant)/|a|,
        // the double root sits at -b/(2a), so discriminant < EPSILON * b * b means the
        // separation is less than sqrt(EPSILON) ~ 1% of the root location). An absolute
        // threshold misclassifies distinct roots of small-extent curves, whose whole
        // quadratic equation lives at a tiny scale, and produces a bogus inflection
        // point with a bogus approximation range.
        // Extrapolating the approximation range from a root outside of the
        // curve would produce points off the curve, so only use it when it
        // is in range and treat the curve as inflection-free otherwise.
        let t = -b / (2.0 * a);

        if in_range(t) {
            *t1 = t;
            return 1;
        }
        return 0;
    }

    // This code is derived from https://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf page 184.
    // Computing the roots this way avoids precision issues when a, c or both are small.
    let discriminant_sqrt = f32::sqrt(discriminant);
    let sign_b = if b >= 0.0 { 1.0 } else { -1.0 };
    let q = -0.5 * (b + sign_b * discriminant_sqrt);
    let mut first_inflection = q / a;
    let mut second_inflection = c / q;

    if first_inflection > second_inflection {
        std::mem::swap(&mut first_inflection, &mut second_inflection);
    }

    let mut next = t1;
    let mut count = 0;

    if in_range(first_inflection) {
        *next = first_inflection;
        next = t2;
        count += 1;
    }

    if in_range(second_inflection) {
        *next = second_inflection;
        count += 1;
    }

    return count;
}

// Find the range around the start of the curve where the curve can locally be approximated
// with a line segment, given a tolerance threshold.
fn inflection_approximation_range(
    bezier: &CubicBezierSegment,
    t: f32,
    tolerance: f32,
    t_min: &mut f32,
    t_max: &mut f32,
) {
    let next_curve = bezier.after_split(t);

    // Transform the curve such that it starts at the origin.
    let mut p1 = next_curve.ctrl1 - next_curve.from;
    let p2 = next_curve.ctrl2 - next_curve.from;
    let p3 = next_curve.to - next_curve.from;

    if p1.x == 0.0 && p1.y == 0.0 {
        p1 = p2
    }

    if p1.x == 0.0 && p1.y == 0.0 {
        // The first three control points coincide: the remaining curve is a
        // straight line parameterized by t^3, so its baseline chord is an
        // exact approximation. The range never extends to the left of the
        // inflection point, which nothing certifies.
        *t_min = t;
        *t_max = 1.0;
        return;
    }

    // Thus, curve(t) = t^3 * (3*p1 - 3*p2 + p3) + t^2 * (-6*p1 + 3*p2) + t * (3*p1).
    // Since curve(0) is an inflection point, cross(p1, p2) = 0, i.e. p1 and p2 are parallel.

    let cross_p2_p3 = p2.cross(p3);
    // The de Casteljau sub-division rounds every control point to a few ulps
    // of the curve's coordinate magnitude, and the cross product itself rounds,
    // so a |cross_p2_p3| below this bound is indistinguishable from the exactly
    // straight remaining curve it is meant to detect. Feeding such a value into
    // the s3 formula below would extrapolate pure rounding noise into an
    // approximation range covering the entire curve.
    let max_coord = bezier
        .from
        .x
        .abs()
        .max(bezier.from.y.abs())
        .max(bezier.ctrl1.x.abs())
        .max(bezier.ctrl1.y.abs())
        .max(bezier.ctrl2.x.abs())
        .max(bezier.ctrl2.y.abs())
        .max(bezier.to.x.abs())
        .max(bezier.to.y.abs());
    let cross_noise = f32::EPSILON
        * (max_coord * (p2.length() + p3.length()) * 8.0 + p2.length() * p3.length() * 4.0);
    if cross_p2_p3.abs() <= cross_noise {
        // s3 == 0 (up to rounding): the paper's model gives s(t) = s3 * t^3 = 0
        // for all t, i.e. the remaining curve stays on the tangent line and its
        // baseline chord is exact. The model assumes a genuine simple inflection
        // point (p1 and p2 parallel), which does not hold for degenerate curves
        // such as a cusp at an endpoint, so verify the chord instead: both
        // control points must stay within tolerance of the tangent line, which
        // bounds the tail's transverse displacement
        //     s(u) = u^3 * y3 + 3 * u^2 * (1 - u) * y2
        // by (3u^2 - 2u^3) * tolerance <= tolerance, with y2/y3 the transverse
        // offsets of p2/p3. The range never extends to the left of the
        // inflection point, which nothing certifies.
        if p1.cross(p2).abs() <= tolerance * p1.length()
            && p1.cross(p3).abs() <= tolerance * p1.length()
        {
            *t_min = t;
            *t_max = 1.0;
        } else {
            // No usable linear range; the caller flattens both sides.
            *t_min = t;
            *t_max = t;
        }
        return;
    }

    let s3 = cross_p2_p3 / p2.length();
    let r_next = (tolerance / s3).abs().powf(1.0 / 3.0);
    let r = r_next * (1.0 - t);
    *t_min = t - r;
    *t_max = t + r;
}
