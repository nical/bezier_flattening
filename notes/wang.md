# Wang's formula (fixed step count)

[Implementation in this repository](../impl/src/algorithms/wang.rs)

Rather than adaptively looking for split points, this algorithm computes
up-front the number of line segments needed to approximate the curve within
the tolerance, then samples the curve at regular parameter intervals.

Wang's formula derives that count from a bound on the deviation of a uniformly
sampled polyline, based on the magnitude of the curve's second differences.
The same count is reused by the [forward differencing](fwd_diff.md) flattener
and, per sub-curve, by the [Yzerman algorithm](yzerman.md).

## Implementations in this repository

- `Wang`: cubic and quadratic.
- `WangSimd4`: cubic-only SIMD variant sampling four parameter values at a
  time.

Bench/edge-count labels: `wang`, `wang-simd`.

## Performance characteristics

- The cost is concentrated up-front, once per curve (a fourth-root
  computation); after that, each subdivision costs only a curve evaluation.
  This is the opposite of [Levien's algorithm](levien.md), which does a fair
  amount of work at every stage (per curve, per sub-curve, and per edge).
- The step count is a conservative bound rather than an estimate, so it
  produces the most segments of the "real" algorithms at equal tolerance (see
  the [edge counts](edge_count.md)), yet it is still the fastest family in the
  benchmarks.
- The sampling loop is regular and independent per sample, which maps well to
  SIMD: `wang-simd` is consistently faster than `wang`.