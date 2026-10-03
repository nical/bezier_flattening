# Flatness criteria

The recursive and linear algorithms are generic over a flatness criterion: a
predicate deciding whether (the remainder of) a curve can be replaced with a
single line segment, given a tolerance. Implementations live in
[`impl/src/algorithms/flatness.rs`](../impl/src/algorithms/flatness.rs):

- `DefaultFlatness`: the default criterion. It measures the control points'
  distance to the baseline, with conservative scaling depending on whether
  the control points are on the same side of the baseline or not, and extra
  checks for control points that do not project onto the baseline.
- `HfdFlatness`: the criterion of the [hybrid forward differencing](hfd.md)
  flattener: the curve is flat when the magnitudes of its second differences
  are below the tolerance.
- `AggFlatness`: the criterion from [Antigrain
  Geometry](https://agg.sourceforge.net/antigrain.com/research/adaptive_bezier/index.html).

TODO: the exact origin of the default criterion is unclear (it is similar to
lyon's and uses fat-line-style bounds). This TODO used to live in
algorithms.md.