# Raph Levien's flattening algorithm

[Implementation in this repository](../impl/src/algorithms/levien.rs) ([simd version](../impl/src/algorithms/levien_simd.rs))


Raph explains the maths behind the flattening of quadratic bézier curves in this blog post: https://raphlinus.github.io/graphics/curves/2019/12/23/flatten-quadbez.html

Flattening cubic bézier curves is done by first approximating them with a series of quadratic curves. However there is a twist: Rather than considering each quadratic bézier sub-curve individually, the algorithm integrates the number of edges that must be produced over multiple sub-curves using a "fractional subdivision" scheme. This avoids the need to insert split points between each sub-curve and can even skip over entire sub-curves.


# Structure of the algorithm

A very simplified version, in pseudo-code:

```rust
fn flatten_cubic(curve, tolerance, callback) {
    let quads_tolerance = tolerance * 0.1;
    let flatten_tolerance = tolerance * 0.9;

    let num_quads = num_quadratics(curve, tolerance);

    // This part is fairly arithmetic-heavy. It maps very well to SIMD.
    for i in 0..num_quads {
        sub_curves.push(flattening_params(curve, flatten_tolerance, i));
    }

    for quad in sub_curves {
        while let Some(u) = fractional_subdivision(quad) {
            // This part is also fairly arithmetic-heavy. It is a little
            // more difficult to optimize with SIMD because for some of
            // the datasets, there are few samples per sub-curves. (See
            // the stats below).
            let t = map_curve_parameter_for(u);
            let p = quad.sample(t)
            callback(p)
        }
    }

    callback(curve.to)
}
```

The tolerance budget is split between the two stages: the quadratic
approximation introduces error of its own, so the flattening stage uses a
stricter tolerance for the total error to stay within the requested tolerance.

# Implementations in this repository

- `Levien`: the scalar implementation. Cubic and quadratic.
- `LevienQuads`: the algorithm without the fractional subdivision scheme; each
  quadratic sub-curve is flattened independently. Cubic and quadratic.
- `LevienSimd`: computes the sub-curve parameters four at a time with SIMD.
  Cubic and quadratic.
- `LevienSimd2` / `LevienSimd3` / `LevienSimdBuf`: so-far unsuccessful
  experiment variants, see
  [levien_experiments.rs](../impl/src/experiments/levien_experiments.rs).
  Cubic only.
- `LevienLinear`: dispatches by curve size, sending small curves to the
  [linear](linear.md) algorithm. Cubic and quadratic.
- `Kurbo`: the [kurbo](https://crates.io/crates/kurbo) crate's implementation
  of the same algorithm (in f64), used as an external reference. Cubic and
  quadratic.

Bench/edge-count labels: `levien`, `levien-quads`, `levien-simd`,
`levien-linear`.

# Performance characteristics

This algorithm spread its computation to every stage: per curve to decide the
number of quadratic sub-curves, per sub-curve to compute its flattening
parameters, and per edge to map the target edge position back to a curve parameter.

The upside is the closest to optimal output: generally the fewest line segments of all
algorithms for the same tolerance. The up-front cost amortizes well for curves
that produce many segments, but not for lots of very small curves; see for
example the [benchmark results for the font dataset](results/bench-cubic-font-threadripper.svg). `LevienLinear` implements the size-threshold idea: small
curves go to the cheaper algorithm.

# Stats

Some stats that were collected to guide optimization:

- [Number of quadratic sub-curves per cubic curve](levien_quads_per_cubics.md)
- [Number of lines per quadratic sub-curve](levien_lines_per_quads.md)
