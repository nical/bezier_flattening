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

In the [edge counts](edge_count.md) (label `hain`) the algorithm produces the
fewest segments of all algorithms at small tolerances (roughly 0.01 to 0.075 on
the tiger dataset), but degrades at larger tolerances where it produces the
most; its speed follows the same pattern in the [benchmarks](bench-cubic-threadripper.md).
Cubic curves only.

# Issues

Unfortunately the paper does not expand on how small is "small enough" for the approximation to work. I spent quite a bit of time trying to get this algorithm to work, but could not get it to produce correct results at tolerance thresholds in the order of 0.25.

The paper notes that testing was done with a flatness threshold of `0.0005`. Eyballing the results for small thresholds (for example `0.01`), it looks like the approximation indeed holds up so my current theory is that this algorithm works at low flatness thresholds (and is quite good at produce as few edges as possible under those parameters) but the approximation does not work for higher flatness thresholds.

TODO: It should not be too difficult to work out the math and show incorrect error approximations with examples.

# Workaround

In order to be fair, the implementation in this repository applies a scaling factor to the flattening step which halves it at `tolerance = 0.25` and higher and gradually scales back to the paper's flattening step for tolerances `0.01` and below. This isn't a great solution.
