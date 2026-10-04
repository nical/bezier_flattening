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

// Collect the error of each curve individually (the curve's maximum deviation
// from its own polyline), so that the test can aggregate max, mean and median
// over the corpus.
fn compute_cubic_curve_errors<F: Flatten>(
    curves: &[CubicBezierSegment],
    tolerance: f32,
) -> Vec<f32> {
    let mut poly = Vec::new();
    let mut errors = Vec::with_capacity(curves.len());
    for curve in curves {
        poly.push(curve.from);
        F::cubic(&curve, tolerance, &mut |seg| {
            poly.push(seg.to);
        });

        errors.push(compute_error(curve, &poly).0);
        poly.clear();
    }

    errors
}

// (max, mean, median) of a per-curve error sample.
fn error_stats(errors: &[f32]) -> (f32, f32, f32) {
    if errors.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut max = 0.0f32;
    let mut sum = 0.0f64;
    for e in errors {
        max = max.max(*e);
        sum += *e as f64;
    }
    let mut sorted = errors.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = sorted.len() / 2;
    let median = if sorted.len() % 2 == 1 {
        sorted[mid]
    } else {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    };

    (max, (sum / sorted.len() as f64) as f32, median)
}

#[test]
fn cubic_error() {
    let curves = generate_bezier_curves();

    // Per algorithm: (name, [max, mean, median] rows, one value per tolerance).
    // The aggregation population is the set of curves: each curve contributes
    // its maximum deviation from its own polyline.
    let mut results: Vec<(&'static str, [Vec<f32>; 3])> = Vec::new();
    results.push(("linear     ", [Vec::new(), Vec::new(), Vec::new()]));
    results.push(("levien     ", [Vec::new(), Vec::new(), Vec::new()]));
    results.push(("levien-simd", [Vec::new(), Vec::new(), Vec::new()]));
    results.push(("levien-linear", [Vec::new(), Vec::new(), Vec::new()]));
    results.push(("wang", [Vec::new(), Vec::new(), Vec::new()]));
    results.push(("wang-simd", [Vec::new(), Vec::new(), Vec::new()]));
    results.push(("hain", [Vec::new(), Vec::new(), Vec::new()]));

    fn measure<F: Flatten>(
        curves: &[CubicBezierSegment],
        tolerance: f32,
        rows: &mut [Vec<f32>; 3],
    ) {
        let errors = compute_cubic_curve_errors::<F>(curves, tolerance);
        let (max, mean, median) = error_stats(&errors);
        rows[0].push(max);
        rows[1].push(mean);
        rows[2].push(median);
    }

    for tolerance in crate::TOLERANCES {
        measure::<crate::Linear>(&curves, tolerance, &mut results[0].1);
        measure::<crate::Levien>(&curves, tolerance, &mut results[1].1);
        measure::<crate::LevienSimd>(&curves, tolerance, &mut results[2].1);
        measure::<crate::LevienLinear>(&curves, tolerance, &mut results[3].1);
        measure::<crate::Wang>(&curves, tolerance, &mut results[4].1);
        measure::<crate::WangSimd4>(&curves, tolerance, &mut results[5].1);
        measure::<crate::Hain>(&curves, tolerance, &mut results[6].1);
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

    let metric_names = ["max error", "mean error", "median error"];
    for (metric, metric_name) in metric_names.iter().enumerate() {
        let _ = writeln!(output, "{}", metric_name);
        print_first_row_md(output);
        for (name, rows) in &results {
            print_row_md(output, *name, &rows[metric]);
        }
        let _ = writeln!(output);
    }
}

