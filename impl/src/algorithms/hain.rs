/// The algorithm implemented here is based on:
///     "Fast, precise flattening of cubic Bézier path and offset curves"
///     http://sibgrapi.sid.inpe.br/col/sid.inpe.br/banon/2004/08.13.18.12/doc/BezierOffsetRendering.pdf
///
/// The basic premise is that for a small t the third order term in the
/// equation of a cubic bezier curve is insignificantly small. This can
/// then be approximated by a quadratic equation for which the maximum
/// difference from a linear approximation can be much more easily determined.
///
/// The step size below follows the paper's quadratic bound; the cases the
/// premise does not cover are handled conservatively, and the quality caveats
/// that follow from the premise are documented in notes/hain.md.
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

// Distance from the curve's inner control points to the chord, which bounds
// the distance from the curve to the chord: the curve lies within the convex
// hull of its control points, the chord's endpoints are two of the hull's
// points, and the distance to a convex set is a convex function, so the curve
// stays within tolerance of the chord whenever its inner control points do.
fn chord_within_tolerance(bezier: &CubicBezierSegment, tolerance: f32) -> bool {
    let v = bezier.to - bezier.from;
    let vv = v.dot(v);
    let mut max_dist = 0.0f32;
    for p in [bezier.ctrl1, bezier.ctrl2] {
        let u = if vv > 0.0 { (p - bezier.from).dot(v) / vv } else { 0.0 };
        let proj = bezier.from + v * u.clamp(0.0, 1.0);
        max_dist = max_dist.max((p - proj).length());
    }
    max_dist <= tolerance
}

fn no_inflection_flattening_step(bezier: &CubicBezierSegment, tolerance: f32) -> f32 {
    let v1 = bezier.ctrl1 - bezier.from;
    let v2 = bezier.ctrl2 - bezier.from;

    // This function assumes that the bézier segment is not starting at an inflection point,
    // otherwise the following cross product may result in very small numbers which will hit
    // floating point precision issues.

    // The paper transforms the curve into the frame of its start tangent and
    // keeps only the quadratic transverse term, arguing that the cubic term is
    // insignificantly small for a small enough step:
    //     s(u) = s2 * u^2,     s2 = 3 * cross(v2, v1) / |v1|.
    // The largest deviation of a quadratic from the chord over [0, t] is
    // |s2| * t^2 / 4, reached at u = t/2, so the paper flattens the [0, t]
    // portion of the curve as a straight line with t = 2 * sqrt(tolerance / |s2|).
    //
    // Everything the premise drops is a source of over-tolerant steps at the
    // looser tolerances used here: the cubic transverse term, and the
    // longitudinal motion that makes the chord deviate from the tangent line
    // the transverse deviation is measured against. The paper targets much
    // tighter tolerances (around 0.0005) where these are negligible; see
    // notes/hain.md for the measured impact at 0.01 to 1.0.
    let v2_cross_v1 = v2.cross(v1);
    if v2_cross_v1 == 0.0 {
        // The transverse deviation of ctrl2 rounds to zero: either the curve
        // stays on its start tangent and the chord is an exact approximation,
        // or the cross product of two nearly parallel (possibly long) vectors
        // was rounded to zero for a curve that does deviate, which the
        // quadratic model cannot see.
        if chord_within_tolerance(bezier, tolerance) {
            return 1.0;
        }
        // The model is blind to this curve: take a conservative step to make
        // progress. The sub-curves' own frames typically leave this degenerate
        // configuration after a split or two.
        return 0.5;
    }

    // To remove divisions and check for divide-by-zero, this is optimized from:
    // s2 = 3. * (v2.x * v1.y - v2.y * v1.x) / hypot(v1.x, v1.y);
    // t = 2. * sqrt(tolerance / abs(s2));
    let s2inv = (v1.x * v1.x + v1.y * v1.y).sqrt() / (3.0 * v2_cross_v1);
    let t = 2.0 * f32::sqrt(tolerance * f32::abs(s2inv));

    // TODO: We start having floating point precision issues if this constant
    // is closer to 1.0 with a small enough tolerance threshold.
    if t >= 0.995 || t <= 0.0 {
        // The closed-form step claims the whole remaining curve is within
        // tolerance. That claim inherits every blind spot of the quadratic
        // model (the ignored cubic term, the longitudinal motion, and cross
        // products of nearly parallel vectors rounding to zero), so verify it
        // against the convex hull before emitting the whole curve as a line.
        if chord_within_tolerance(bezier, tolerance) {
            return 1.0;
        }
        return 0.5;
    }

    t
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
