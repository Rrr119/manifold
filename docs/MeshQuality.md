# Mesh quality cleanup

Manifold's Boolean operations guarantee a valid, oriented 2-manifold and the
correct solid volume. That contract does not imply a well-shaped triangulation:
near-coincident input faces, tiny intersection edges, and narrow clipped regions
can leave needles and slivers. `Manifold::ImproveQuality()` is an explicit
post-processing pass for those cases. It does not change the result of the
Boolean itself.

## API

```cpp
Manifold improved = result.ImproveQuality(tolerance);
Manifold improved = result.Cleanup(tolerance);  // alias
Manifold improved = result.RemoveDegenerates(tolerance);  // tolerance-aware overload
```

The same operation is exposed as:

- C: `manifold_improve_quality()` and `manifold_cleanup()`
- Python: `manifold.improve_quality()` and `manifold.cleanup()`
- TypeScript/JavaScript: `manifold.improveQuality()` and `manifold.cleanup()`

`tolerance` has units of length and defaults to the manifold's current epsilon.
Increasing it collapses longer short edges. It is not a smoothing radius and it
does not simplify well-shaped regions.

## What the pass does

1. **Removes degenerate triangles.** It reuses the same edge-collapse and
   edge-swap machinery as Boolean cleanup.
2. **Collapses short edges.** Edges no longer than `tolerance` are candidates
   only when the quadric collapse error is also within `tolerance`. Collapses
   are sorted by error and rejected if they would invert or overly distort a
   neighboring triangle.
3. **Improves triangle shape.** Every coplanar pair of adjacent triangles is
   tested with the alternative diagonal. The diagonal is flipped only when it
   raises the worse of the two triangle qualities. Quality is normalized so
   that 1 is equilateral and 0 is degenerate:

   ```text
   q = 2 * sqrt(3) * |(b - a) x (c - a)| / (|b - a|^2 + |c - b|^2 + |a - c|^2)
   ```

4. **Restores internal invariants.** Face/vertex normals are recomputed and the
   result is re-sorted into the usual compact Manifold layout.

## Guarantees

- The output is a valid oriented 2-manifold. Status is propagated from the
  input, and a round trip through `MeshGL`/`MeshGL64` remains valid.
- Coplanar edge flips keep exactly the same quadrilateral surface, so they do
  not change volume or area. Only flips whose local volume change is within
  numerical epsilon are accepted.
- A short-edge collapse moves the surface by at most `tolerance` and is only
  accepted when its quadric error is within `tolerance` as well.
- Non-coplanar quads are not flipped. A sliver whose neighbors cannot form a
  valid convex quad may therefore remain: quality improvement is best effort,
  while manifoldness and volume preservation are hard requirements.
- Vertex properties are carried through the existing edge-collapse and
  edge-swap interpolation paths. Property seams are not deliberately moved.

## Tests and benchmark

`test/quality_test.cpp` covers:

- short-edge and sliver removal on an extruded notch,
- volume, area, genus, and manifold round-trip preservation across a Boolean,
- quality improvement from `Simplify()`, and
- a 30k+ triangle Boolean workload used as the benchmark signal.

`Quality.ImproveQualityBenchmark` is included in the weekly benchmark suite by
`scripts/run_weekly_benchmarks.py`, so its runtime is tracked alongside the
existing Boolean and tessellation benchmarks.

## Example

```cpp
#include "manifold/manifold.h"

using namespace manifold;

Manifold MakeBracket() {
  Manifold body = Manifold::Cube({40, 20, 8}, true);
  Manifold hole = Manifold::Cylinder(12, 4, 4, 64);
  Manifold result = body - hole.Rotate(90, 0, 0);

  // Clean up the clipped rim without changing the solid.
  return result.ImproveQuality(result.GetEpsilon() * 100);
}
```
