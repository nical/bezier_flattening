# Hybrid forward differencing

[Implementation in this repository](../impl/src/algorithms/hybrid_fwd_diff.rs)

An adaptive version of [forward differencing](fwd_diff.md): the step size can
change while flattening, but only by powers of two (halving or doubling). That
restriction is what makes the scheme practical: the curve is expressed in a
basis for which stepping, halving and doubling the step all have cheap exact
formulas, and the second-difference terms of that basis directly measure the
error, doubling as the flatness test.

At each step the flattener halves the step if the error exceeds the tolerance
(a single halving is provably sufficient), or doubles it repeatedly while the
error stays well within the tolerance.

## References

The algorithm is described in US patent 5367617, "System and method of hybrid
forward differencing to render Bezier splines" (Microsoft, 1995), with earlier
academic references:

- Lien, Shantz and Pratt, "Adaptive Forward Differencing for Rendering Curves
  and Surfaces", Computer Graphics, July 1987.
- Chang and Shantz, "Rendering Trimmed NURBS with Adaptive Forward
  Differencing", Computer Graphics, August 1988.

The code is an adaptation of Jeff Muizelaar's Rust port of WPF's C++
implementation (WPF's rasterizer predates Direct2D).

## Implementations in this repository

- `HybridFwdDiff`: cubic only. Appears in the [edge counts](edge_count.md)
  under the label `hfd` but is not benchmarked.