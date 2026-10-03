# Forward differencing

[Implementation in this repository](../impl/src/algorithms/fwd_diff.rs)

Forward differencing evaluates the curve incrementally: instead of evaluating
the polynomial at each parameter value, each successive point is obtained by
evolving the value and difference accumulators with a few additions.

The step count is fixed up-front using the same segment count as the
[Wang](wang.md) flattener, so the step size does not adapt to the curve's
shape.

## Caveats

- Positions accumulate additions over the whole curve, so floating point
  rounding errors compound over many steps.

## Implementations in this repository

- `FwdDiff`: cubic and quadratic.

Bench/edge-count label: `fwd-diff`.

## Performance characteristics

- Produces exactly the same segment counts as `Wang` (same up-front count),
  see the [edge counts](edge_count.md).
- The per-subdivision cost is a few additions, lower than evaluating the
  curve; in the benchmarks it is on par with `wang`, and both are the fastest
  family by a wide margin.
- The adaptive counterpart is [hybrid forward differencing](hfd.md), which
  trades a bit of per-step cost for fewer output segments.