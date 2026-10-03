# Recursive subdivision

[Implementation in this repository](../impl/src/algorithms/recursive.rs)

The most straightforward adaptive flattening algorithm: recursively split the
curve in half until each remainder passes a [flatness test](flatness.md),
then emit it as a single line segment.

## Implementations in this repository

- `Recursive`: cubic and quadratic.
- `RecursiveHfd`, `RecursiveAgg`: cubic-only variants using the other
  [flatness criteria](flatness.md). Not benched.

Bench/edge-count label: `recursive`.

## Performance characteristics

- No up-front cost: the work is proportional to the number of produced
  segments, which makes it a good fit for small curves.
- Each subdivision costs a de Casteljau split plus a flatness test.
- Produces more segments than Levien's algorithm for the same tolerance (see
  the [edge counts](edge_count.md)): the dyadic splits cannot adapt to where
  the curve actually needs them.
- Uses the stack for recursion.