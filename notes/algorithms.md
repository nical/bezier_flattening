# Algorithms

The notes are organized per algorithm family. An algorithm may have several
implementations in this repository (scalar, SIMD variants, experiments, or an
external reference implementation); each note lists them with the labels used
in the [benchmarks](readme.md) and [edge counts](edge_count.md).

- [Recursive subdivision](recursive.md): the classic adaptive bisection
  flattening, generic over a [flatness criterion](flatness.md).
  Implementations: `Recursive`, plus the `RecursiveHfd` / `RecursiveAgg`
  criterion variants.
- [Linear scanning](linear.md): a non-recursive variant that repeatedly emits
  the longest flat prefix of the remaining curve, also generic over the
  [flatness criterion](flatness.md).
  Implementations: `Linear`, plus the `LinearHfd` / `LinearAgg` variants.
- [Levien's algorithm](levien.md): approximates cubics with quadratics, and
  flattens them with a fractional subdivision scheme.
  Implementations: `Levien`, `LevienQuads`, `LevienSimd`, the
  `LevienSimd2/3/Buf` experiments, the size-dispatching hybrid `LevienLinear`,
  and `Kurbo` (the kurbo crate's implementation, used as an external
  reference).
- [Wang's formula](wang.md): computes a fixed segment count up-front and
  samples at regular parameter intervals.
  Implementations: `Wang`, `WangSimd4`.
- [Forward differencing](fwd_diff.md): evaluates the curve incrementally with
  a fixed step count (the same count as Wang).
  Implementations: `FwdDiff`.
- [Hybrid forward differencing](hfd.md): the adaptive variant of forward
  differencing, with step halving/doubling; the basis of WPF's rasterizer.
  Implementations: `HybridFwdDiff`.
- [Yzerman](yzerman.md): approximates cubics with quadratics, each flattened
  with Wang's segment count.
  Implementations: `Yzerman`, `YzermanSimd4`.
- [Hain](hain.md): flattens between inflection points with the paper's
  quadratic transverse step bound.
  Implementations: `Hain`.
- [Fixed splitting](fixed.md): tolerance-ignoring baselines.
  Implementations: `Fixed1`, `Fixed16`.

# Caveats

- A fair amount of effort has already gone into speeding up `levien` using 128-bit wide SIMD intrinsics. More effort is underway, most of the other algorithms did not get much performance work yet, so numbers will probably change.