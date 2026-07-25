use crate::{CubicBezierSegment, LineSegment, QuadraticBezierSegment, point};

use crate::{fast_ceil, fast_recip};

// Avoid two square roots using a lookup table that contains
// i^4 for  i in 1..25.
const N: usize = 24;
const LUT: [f32; N] = [
    1.0, 16.0, 81.0, 256.0, 625.0, 1296.0, 2401.0, 4096.0, 6561.0,
    10000.0, 14641.0, 20736.0, 28561.0, 38416.0, 50625.0, 65536.0,
    83521.0, 104976.0, 130321.0, 160000.0, 194481.0, 234256.0,
    279841.0, 331776.0
];

/// Computes the number of line segments required to build a flattened approximation
/// of the curve with segments placed at regular `t` intervals.
pub fn num_segments_cubic(curve: &CubicBezierSegment, tolerance: f32) -> f32 {
    let from = curve.from.to_vector();
    let ctrl1 = curve.ctrl1.to_vector();
    let ctrl2 = curve.ctrl2.to_vector();
    let to = curve.to.to_vector();
    let v1 = (from - ctrl1 * 2.0 + ctrl2) * 6.0;
    let v2 = (ctrl1 - ctrl2 * 2.0 + to) * 6.0;
    let l = v1.dot(v1).max(v2.dot(v2));
    let d = fast_recip(8.0 * tolerance);
    let err4 = l * d * d;

    // If the value we are looking for is within the LUT, take the fast path
    if err4 <= 331776.0 {
        #[allow(clippy::needless_range_loop)]
        for i in 0..N {
            if err4 <= LUT[i] {
                return (i + 1) as f32;
            }
        }
    }

    // Otherwise fall back to computing via two square roots.
    fast_ceil(err4.sqrt().sqrt()).max(1.0)
}

/// Computes the number of line segments required to build a flattened approximation
/// of the curve with segments placed at regular `t` intervals.
pub fn num_segments_quadratic(curve: &QuadraticBezierSegment, tolerance: f32) -> f32 {
    let from = curve.from.to_vector();
    let ctrl = curve.ctrl.to_vector();
    let to = curve.to.to_vector();
    let l = (from - ctrl * 2.0 + to) * 2.0;
    let l2 = l.dot(l);
    let d = fast_recip(8.0 * tolerance);
    let err4 = l2 * d * d;

    // If the value we are looking for is within the LUT, take the fast path
    if err4 <= 331776.0 {
        #[allow(clippy::needless_range_loop)]
        for i in 0..N {
            if err4 <= LUT[i] {
                return (i + 1) as f32;
            }
        }
    }

    // Otherwise fall back to computing via two square roots.
    fast_ceil(err4.sqrt().sqrt()).max(1.0)
}


/// Flatten the curve by precomputing a number of segments and splitting the curve
/// at regular `t` intervals.
pub fn flatten_cubic<F>(curve: &CubicBezierSegment, tolerance: f32, callback: &mut F)
    where
    F:  FnMut(&LineSegment)
{
    let poly = crate::polynomial_form_cubic(&curve);
    let n = num_segments_cubic(curve, tolerance);
    let step = fast_recip(n);
    let mut prev = 0.0;
    let mut from = curve.from;
    for _ in 0..(n as u32 - 1) {
        let t = prev + step;
        let to = poly.sample_fma(t);
        callback(&mut LineSegment { from, to });
        from = to;
        prev = t;
    }

    let to = curve.to;
    callback(&mut LineSegment { from, to });
}

/// Flatten the curve by precomputing a number of segments and splitting the curve
/// at regular `t` intervals.
pub fn flatten_quadratic<F>(curve: &QuadraticBezierSegment, tolerance: f32, callback: &mut F)
    where
    F:  FnMut(&LineSegment)
{
    let poly = crate::polynomial_form_quadratic(curve);
    let n = num_segments_quadratic(curve, tolerance);
    let step = fast_recip(n);
    let mut prev = 0.0;
    let mut from = curve.from;
    for _ in 0..(n as u32 - 1) {
        let t = prev + step;
        let to = poly.sample(t);
        callback(&mut LineSegment { from, to });
        from = to;
        prev = t;
    }

    let to = curve.to;
    callback(&mut LineSegment { from, to });
}


#[cfg_attr(target_arch = "x86_64", target_feature(enable = "avx"))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "fma"))]
pub unsafe fn flatten_cubic_simd4<F>(curve: &CubicBezierSegment, tolerance: f32, callback: &mut F)
    where
    F:  FnMut(&LineSegment)
{
    use crate::simd4::{vec4, splat, add, mul};

    let poly = crate::polynomial_form_cubic(&curve);
    let n = num_segments_cubic(curve, tolerance);
    let mut from = curve.from;

    let a0x = splat(poly.a0.x);
    let a0y = splat(poly.a0.y);
    let a1x = splat(poly.a1.x);
    let a1y = splat(poly.a1.y);
    let a2x = splat(poly.a2.x);
    let a2y = splat(poly.a2.y);
    let a3x = splat(poly.a3.x);
    let a3y = splat(poly.a3.y);
    let step = fast_recip(n);
    let step = splat(step);
    let mut t = mul(step, vec4(1.0, 2.0, 3.0, 4.0));
    let step4 = mul(step, splat(4.0));

    // minus one because we'll add the last point explicitly
    let mut n = n as i32 - 1;
    // Only the last iteration can be partial, so keep the hot loop free of the
    // `take` adapter and the per-element counter it needs.
    while n >= 4 {
        let (x, y) = crate::simd4::sample_cubic_horner_simd4(a0x, a0y, a1x, a1y, a2x, a2y, a3x, a3y, t);

        let x: [f32; 4] = std::mem::transmute(x);
        let y: [f32; 4] = std::mem::transmute(y);

        for (x, y) in x.iter().zip(y.iter()) {
            let p = point(*x, *y);
            callback(&LineSegment { from, to: p });
            from = p;
        }

        t = add(t, step4);
        n -= 4;
    }

    if n > 0 {
        let (x, y) = crate::simd4::sample_cubic_horner_simd4(a0x, a0y, a1x, a1y, a2x, a2y, a3x, a3y, t);

        let x: [f32; 4] = std::mem::transmute(x);
        let y: [f32; 4] = std::mem::transmute(y);

        for (x, y) in x.iter().zip(y.iter()).take(n as usize) {
            let p = point(*x, *y);
            callback(&LineSegment { from, to: p });
            from = p;
        }
    }

    let to = curve.to;
    callback(&mut LineSegment { from, to });
}

#[cfg_attr(target_arch = "x86_64", target_feature(enable = "avx"))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "fma"))]
pub unsafe fn flatten_quadratic_simd4<F>(curve: &QuadraticBezierSegment, tolerance: f32, callback: &mut F)
    where
    F:  FnMut(&LineSegment)
{
    use crate::simd4::{vec4, splat, add, mul};

    let poly = crate::polynomial_form_quadratic(&curve);
    let n = num_segments_quadratic(curve, tolerance);
    let mut from = curve.from;

    let a0x = splat(poly.a0.x);
    let a0y = splat(poly.a0.y);
    let a1x = splat(poly.a1.x);
    let a1y = splat(poly.a1.y);
    let a2x = splat(poly.a2.x);
    let a2y = splat(poly.a2.y);
    let step = fast_recip(n);
    let step = splat(step);
    let mut t = mul(step, vec4(1.0, 2.0, 3.0, 4.0));
    let step4 = mul(step, splat(4.0));

    // minus one because we'll add the last point explicitly
    let mut n = n as i32 - 1;
    // Only the last iteration can be partial, so keep the hot loop free of the
    // `take` adapter and the per-element counter it needs.
    while n >= 4 {
        let x = crate::simd4::sample_quadratic_horner_simd4(a0x, a1x, a2x, t);
        let y = crate::simd4::sample_quadratic_horner_simd4(a0y, a1y, a2y, t);

        let x: [f32; 4] = std::mem::transmute(x);
        let y: [f32; 4] = std::mem::transmute(y);

        for (x, y) in x.iter().zip(y.iter()) {
            let p = point(*x, *y);
            callback(&LineSegment { from, to: p });
            from = p;
        }

        t = add(t, step4);
        n -= 4;
    }

    if n > 0 {
        let x = crate::simd4::sample_quadratic_horner_simd4(a0x, a1x, a2x, t);
        let y = crate::simd4::sample_quadratic_horner_simd4(a0y, a1y, a2y, t);

        let x: [f32; 4] = std::mem::transmute(x);
        let y: [f32; 4] = std::mem::transmute(y);

        for (x, y) in x.iter().zip(y.iter()).take(n as usize) {
            let p = point(*x, *y);
            callback(&LineSegment { from, to: p });
            from = p;
        }
    }

    let to = curve.to;
    callback(&mut LineSegment { from, to });
}

#[test]
fn wang_simd_scalar() {
    let curve = CubicBezierSegment {
        from: point(100.0, 100.0),
        ctrl1: point(0.0, 100.0),
        ctrl2: point(100.0, 0.0),
        to: point(10.0, 0.0),
    };

    let mut scalar = Vec::new();
    flatten_cubic(&curve, 0.001, &mut |seg| {
        scalar.push(seg.to)
    });

    let mut simd = Vec::new();
    unsafe {
        flatten_cubic_simd4(&curve, 0.001, &mut |seg| {
            simd.push(seg.to)
        });
    }

    assert_eq!(scalar.len(), simd.len());
    for (i, (scalar, simd)) in scalar.iter().zip(simd.iter()).enumerate() {
        assert!(scalar.distance_to(*simd) < 0.01, "{i:?}: scalar {scalar:?} simd {simd:?}");
    }
}
