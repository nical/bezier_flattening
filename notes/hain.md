# Fast, Precise Flattening of Cubic Bézier Segment Offset Curves

[Paper](http://sibgrapi.sid.inpe.br/col/sid.inpe.br/banon/2004/08.13.18.12/doc/BezierOffsetRendering.pdf)
[Implementation in this repository](../impl/src/algorithms/hain.rs)


The idea behind this algorithm is that for a "small enough step", the third order term in the equation of a cubic bezier curve is insignificantly small. This can then be approximated by a quadratic equation for which the maximum difference from a linear approximation can be much more easily determined.

# Summary

In a nutshell, the algorithm is based on a circular approximation of the curve that is applied iteratively. This approximation fails near inflection points, so parts near the inflection points are handled separately and the sub-curve segments in-between are flattened using an iterative process in which at each step, the longest part of the curve that can be approximated with a line segment is computed, pushed into the result and removed from the curve.

# Performance

The up-front work is an inflection point analysis of the curve (solving for the
inflection parameters with a numerically stable quadratic root solver), after
which the sub-curves between inflection ranges are flattened iteratively by
extracting the longest flattenable prefix at each step.

In the [edge counts](edge_count.md) (label `hain`) the segment counts sit within
about 10% of the levien family at small tolerances and are the lowest of all
algorithms at 0.5 and 1.0 on the full corpus. Unlike the levien variants, whose
counts come with a certified error that exceeds the tolerance at several
thresholds, hain's error stays within tolerance across the sweep. Cubic curves
only.

Certification is not free: computing the exact maximum chord distance takes a
handful of square and cube roots per step where the paper's estimate needed a
single polynomial evaluation, and hain is the slowest algorithm in every cubic
[benchmark](bench-cubic-threadripper.md) table. At tolerance 0.2 it spends
roughly 20 to 30 times the levien family's time per emitted segment on the
chord-heavy datasets, where the paper's estimate cost 2.5 to 3.5 times; on
datasets of short curves fixed per-curve costs dominate and the gap is much
smaller.

# Issues

The paper derives its step size from the transverse displacement alone,
assuming the longitudinal position is `x(u) ~ 3u|v1|`, that the cubic term is
negligible for a "small enough" step, and that the next split lands at `2 * t0`.
That under-estimates the true chord distance wherever the start tangent does
not point along the direction of travel (near cusps the true distance is up to
twice the approximation), and the error grows with the tolerance as the
"small enough" assumption weakens. With the paper's step size the flattening
exceeded the tolerance at every threshold.

The implementation therefore no longer uses the paper's step size. At each step
it computes the exact maximum distance from the curve to the candidate chord
(the closed-form cubic `u (t - u) |A + B u| / |B(t)|`, maximized by solving a
quadratic in the substitution `v = u/t`, `z = v - 1/2` for numerical stability)
and takes the longest prefix whose certified distance meets the tolerance, with
part of the budget reserved for the f32 rounding of the split itself.
