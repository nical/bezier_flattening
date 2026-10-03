# Linear (sequential) flattening

[Implementation in this repository](../impl/src/algorithms/linear.rs)

A non-recursive variant of [recursive subdivision](recursive.md): instead of
splitting the curve in halves down a recursion tree, the algorithm repeatedly
finds the longest flat prefix of the remaining curve, emits it, and advances
toward the end of the curve. The search for each split point starts optimistic
and refines by doubling after a success and halving after a failure, so
flatness tests amortize over runs of similarly sized segments.

## Implementations in this repository

- `Linear`: cubic and quadratic.
- `LinearHfd`, `LinearAgg`: cubic-only variants using the other
  [flatness criteria](flatness.md). Not benched.

Bench/edge-count label: `linear`.

## Performance characteristics

- Like the recursive version, has no up-front cost: work is proportional to
  the number of produced segments.
- Avoids the recursion, and the split search amortizes flatness tests over
  several segments.
- Produces slightly fewer segments than the recursive algorithm: split points
  are not restricted to dyadic positions of the original curve.
- In the benchmarks the two trade places depending on dataset and tolerance;
  see the [benchmark results](bench-cubic-threadripper.md).