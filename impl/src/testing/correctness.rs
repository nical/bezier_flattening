#[cfg(test)]
use std::path::PathBuf;

use crate::{CubicBezierSegment, Point};
#[cfg(test)]
use crate::testing::table::{print_first_row_md, print_row_md};
#[cfg(test)]
use crate::testing::generate_bezier_curves;
use crate::Flatten;

fn compute_error(curve: &CubicBezierSegment, approximation: &[Point]) -> (f32, f32) {
    let mut quads = Vec::new();
    curve.for_each_quadratic_bezier(0.0001, &mut |quad| {
        quads.push(quad.to_f64());
    });

    let mut sum_error: f32 = 0.0;
    let mut max_error: f32 = 0.0;
    let mut count = 0.0;
    let mut prev: Option<Point> = None;
    for to in approximation {
        if let Some(from) = prev {
            let mid = from.to_f64().lerp(to.to_f64(), 0.5);
            let mut min_dist: f32 = 100000000.0;
            for quad in &quads {
                min_dist = min_dist.min(quad_distance_to_point(quad, mid) as f32);
            }
            max_error = max_error.max(min_dist);
            sum_error += min_dist;
            count += 1.0;
        }
        prev = Some(*to);
    }

    // The midpoint metric above is blind to a dropped tail: a polyline that
    // stops short of the curve's end point still has well-placed segment
    // midpoints. Charge the uncovered gap as error.
    let tail_error = match approximation.last() {
        Some(p) => (p.to_f64() - curve.to.to_f64()).length(),
        None => (curve.from.to_f64() - curve.to.to_f64()).length(),
    };
    max_error = max_error.max(tail_error as f32);

    let avg = sum_error / count;

    (max_error, avg)
}

// Distance from a point to a quadratic bézier curve.
//
// QuadraticBezierSegment::distance_to_point finds the closest point by solving
// the perpendicularity cubic with Cardano's formula, which loses precision when
// the quad is nearly straight (the quadratic coefficient nearly cancels): the
// root can come out far enough from the true one to inflate a near-zero
// distance by orders of magnitude. Refine the root with Newton iterations on
// (Q(t) - pos) · Q'(t) and keep whichever of the two roots is closer.
fn quad_distance_to_point(
    quad: &lyon_path::geom::QuadraticBezierSegment<f64>,
    pos: lyon_path::geom::euclid::default::Point2D<f64>,
) -> f64 {
    let t0 = quad.closest_point(pos);
    let d0 = (quad.sample(t0) - pos).length();

    let v = quad.ctrl - quad.from;
    let c = quad.from + quad.to.to_vector() - quad.ctrl * 2.0;
    let mut t = t0;
    for _ in 0..8 {
        let p = quad.sample(t);
        let dp = p - pos;
        // Q'(t) = 2 * (v + t * c), Q''(t) = 2 * c
        let qp = v + c * t;
        let g = dp.dot(qp);
        let gp = qp.dot(qp) * 2.0 + dp.dot(c) * 2.0;
        if gp == 0.0 {
            break;
        }
        let next = (t - g / gp).max(0.0).min(1.0);
        if next == t {
            break;
        }
        t = next;
    }

    d0.min((quad.sample(t) - pos).length())
}

fn compute_cubic_error<F: Flatten>(curves: &[CubicBezierSegment], tolerance: f32) -> f32 {
    let mut poly = Vec::new();
    let mut max_error: f32 = 0.0;
    for curve in curves {
        poly.push(curve.from);
        F::cubic(&curve, tolerance, &mut |seg| {
            poly.push(seg.to);
        });

        max_error = max_error.max(compute_error(curve, &poly).0);
        poly.clear();
    }

    max_error
}

#[test]
fn cubic_error() {
    let curves = generate_bezier_curves();

    let mut linear: Vec<f32> = Vec::new();
    let mut levien: Vec<f32> = Vec::new();
    let mut levien_simd: Vec<f32> = Vec::new();
    let mut levien_linear: Vec<f32> = Vec::new();
    let mut wang: Vec<f32> = Vec::new();
    let mut wang_simd: Vec<f32> = Vec::new();
    let mut hain: Vec<f32> = Vec::new();

    for tolerance in crate::TOLERANCES {
        linear.push(compute_cubic_error::<crate::Linear>(&curves, tolerance));
        levien.push(compute_cubic_error::<crate::Levien>(&curves, tolerance));
        levien_simd.push(compute_cubic_error::<crate::LevienSimd>(&curves, tolerance));
        levien_linear.push(compute_cubic_error::<crate::LevienLinear>(&curves, tolerance));
        wang.push(compute_cubic_error::<crate::Wang>(&curves, tolerance));
        wang_simd.push(compute_cubic_error::<crate::WangSimd4>(&curves, tolerance));
        hain.push(compute_cubic_error::<crate::Hain>(&curves, tolerance));
    }

    let out_name = crate::testing::table::get_flatten_output();

    let mut _std_out = None;
    let mut _out_file = None;
    let output: &mut dyn std::io::Write = match out_name {
        None => {
            _std_out = Some(std::io::stdout().lock());
            _std_out.as_mut().unwrap()
        }
        Some(file_name) => {
            let path: PathBuf = file_name.into();
            _out_file = Some(std::fs::File::create(path).unwrap());
            _out_file.as_mut().unwrap()
        }
    };

    print_first_row_md(output);
    print_row_md(output, "linear     ", &linear);
    print_row_md(output, "levien     ", &levien);
    print_row_md(output, "levien-simd", &levien_simd);
    print_row_md(output, "levien-linear", &levien_linear);
    print_row_md(output, "wang", &wang);
    print_row_md(output, "wang-simd", &wang_simd);
    print_row_md(output, "hain", &hain);
}

