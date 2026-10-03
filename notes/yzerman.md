# Yzerman

[Paper](https://blend2d.com/research/simplify_and_offset_bezier_curves.pdf):
"Fast approaches to simplify and offset Bézier curves within specified error
limits", Fabian Yzerman.
[Implementation in this repository](../impl/src/algorithms/yzerman.rs)

Two-stage flattening, similar in spirit to [Levien's algorithm](levien.md):
first approximate the cubic curve with a few quadratic sub-curves (using a
fraction of the error budget), then flatten each quadratic with the
[Wang](wang.md) segment count.

A shortcut in the scalar version emits sub-curves that already pass the
flatness test as single segments; the SIMD version drops it.

## Implementations in this repository

- `Yzerman`: cubic only.
- `YzermanSimd4`: cubic-only SIMD variant computing the sub-curve parameters
  four at a time.

Bench/edge-count labels: `yzerman`, `yzerman-simd`.

## Performance characteristics

- Cheaper setup than Levien's algorithm: no integral approximations and no
  fractional subdivision logic, while keeping the two-stage structure.
- Several times faster than `levien` and somewhat slower than `wang` in the
  benchmarks.
- Its edge counts land between `wang` and `levien`: the cubic→quadratic
  approximation costs edges compared to Levien's fractional subdivision, see
  the [edge counts](edge_count.md).