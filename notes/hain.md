# Fast, Precise Flattening of Cubic Bézier Segment Offset Curves

[Paper](http://sibgrapi.sid.inpe.br/col/sid.inpe.br/banon/2004/08.13.18.12/doc/BezierOffsetRendering.pdf)
[Implementation in this repository](../impl/src/algorithms/hain.rs)


The idea behind this algorithm is that for a "small enough step", the third order term in the equation of a cubic bezier curve is insignificantly small. This can then be approximated by a quadratic equation for which the maximum difference from a linear approximation can be much more easily determined.

# Summary

In a nutshell, the algorithm is based on a circular approximation of the curve that is applied iteratively. This approximation fails near inflection points, so parts near the inflection points are handled separately and the sub-curve segments in-between are flattened iteratively: at each step, the longest prefix that the approximation deems flattenable with a line segment is computed, pushed into the result and removed from the curve.

# The step size

The implementation follows the paper's formula. In the frame of the start tangent, the transverse displacement of the curve is approximated by its quadratic term only,

    s(u) = s2 * u^2,    s2 = 3 * cross(v2, v1) / |v1|,

whose largest deviation from the chord over [0, t] is |s2| * t^2 / 4, reached at u = t/2. Setting that deviation to the tolerance gives the split at t = 2 * sqrt(tolerance / |s2|), which the implementation evaluates in closed form. The caller emits the prefix as a line segment, splits the original curve at the step boundary (a fresh split per chord, rather than chaining sub-curves, keeps the f32 rounding error of each chord endpoint down to the split itself), and repeats on the remainder.

Two cases fall outside the premise and are handled conservatively. When the cross product rounds to zero the model cannot see the curve at all, and when the closed-form step claims the entire remaining curve is within tolerance the claim inherits every blind spot of the model. In both cases the chord is verified against the convex hull of the control points: the curve lies within that hull and the distance to a convex set is a convex function, so control points within tolerance of the chord certify it. Otherwise the implementation splits the curve in half and continues.

An earlier revision scaled the paper's step by a tolerance-dependent fudge factor between 1 and 2 to compensate for the quality problems described below. That factor is not part of the paper and the current implementation does not use it.

# Performance

The up-front work is an inflection point analysis of the curve (a numerically stable quadratic root solve), after which each step costs a couple of subtractions, a cross product and a square root. On the full corpus ([table](results/edge-count-cubic-all.md)) hain emits the fewest segments of all algorithms at every tolerance, about 6% fewer than the levien family at 0.01. This is not a sign of efficiency: hain emits the fewest segments because its chords are the least certified.

In the [cubic benchmarks](bench-cubic-threadripper.md) hain takes about 1.2 to 3.4 times the levien family's time on the chord-heavy corpora (2.6 times at tolerance 0.01 on the full corpus), or 1.25 to 3.6 times per emitted segment. On the short-curve datasets (fonts, inkscape, tiger) the cheap step estimate and the up-front inflection analysis make it faster than the levien family overall, at 0.35 to 0.8 times the total time.

# Issues

The premise drops the cubic transverse term. With y2 and y3 the transverse offsets of the second control point and of the endpoint in the start tangent's frame, the exact transverse displacement of the curve from the chord is

    u (t - u) * (s2 + s3 * (u + t)),    s3 = y3 - 3 * y2,

so on top of the tolerance-budgeted |s2| * t^2 / 4 the ignored term contributes up to roughly |s3| * t^3 / 2. Nothing in the formula keeps the step small enough for that term to be negligible: it scales with the curve's transverse magnitude, and the looser the tolerance, the larger the steps and the weaker the premise. The paper's own experiments run at a tolerance of 0.0005 where the premise holds; at the tolerances swept here it frequently does not.

Measured with the `cubic_error` correctness test on the full corpus, the maximum error stays between 1.7 and 5.1 units across the entire tolerance sweep, 5 to 170 times the tolerance, and barely tracks the tolerance at all: it tracks the coordinate scale of the input instead. On the small-coordinate fonts dataset the error stays under 0.65 units (within tolerance at 0.5 and 1.0), while on the large-coordinate nehab datasets it reaches 1.7 to 5.2 units at every threshold.

The over-tolerance is not a small calibration offset but a heavy tail over the corpus: the median curve lands at or below the tolerance (1.02 times at tolerance 0.01, 0.52 at 0.1, 0.11 at 1.0) while 57% of curves exceed the tolerance at 0.01, 32% at 0.1 and 20% at 1.0, with a 99th percentile at 1.4 to 2.4 times the tolerance and a long tail up to the maximum. The mean is far less tail-sensitive than the maximum and still ranks hain last at every threshold, 1.2 to 1.4 times the levien family's (`cubic_error` prints max, mean and median tables).

Scaling the tolerance does not fix this. Dividing the tolerance by 8 leaves the maximum error at 1.3 to 4.1 units (at tolerance 0.01 it even rises from 1.7 to 2.6 units, since a smaller tolerance shifts which sub-curves are capped and hull-verified and which take one long unverified chord), and dividing it by 64 still leaves 0.26 to 3.4 units at every threshold while spending about 8 times more segments. The error is set by the geometry of each curve (the s3 term and the small-s2 tail near inflection ranges), not by the tolerance.

Two further approximations add to the gap: the bound only measures displacement transverse to the start tangent, assuming the chord runs along it, which under-estimates the true chord distance near cusps by up to a factor of two; and each f32 split rounds the chord endpoints, adding a small tolerance-independent slack proportional to the coordinate magnitude, which becomes visible at tolerance 0.01.

The degenerate-case handling is not cosmetic. The previous revision of this implementation, which flattened blindly in those cases and chained sub-curves (accumulating one rounding error per chord), measured 60 to 186 units of maximum error under the same harness.

For reference, an experiment in this repository's history (commit 5acf5e9) replaced the paper's estimate with an exact certification of the maximum chord distance at every step. That brought the error within tolerance at every threshold, at roughly 8.6 times the time per batch of curves on the chord-heavy datasets (20 to 30 times the levien family's time per emitted segment at tolerance 0.2). It demonstrates the price of the paper's shortcut on this corpus, and was reverted to keep this implementation faithful to the paper.