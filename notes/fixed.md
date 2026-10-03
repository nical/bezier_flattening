# Fixed splitting (baseline)

Implementation in this repository:
[`impl/src/traits.rs`](../impl/src/traits.rs)

Control "algorithms" that ignore the tolerance parameter, used as baselines:

- `Fixed1`: approximates every curve with its baseline (a single segment).
- `Fixed16`: splits every curve into 16 segments at regular parameter
  intervals.

Useful as a sanity floor in the [edge counts](edge_count.md) and for checking
that an adaptive algorithm's overhead is justified on trivial curves. Neither
is currently benched.