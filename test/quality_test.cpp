// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#include "gtest/gtest.h"
#include "manifold/manifold.h"
#include "test.h"

using namespace manifold;

namespace {

// 1 for an equilateral triangle, approaching 0 for a degenerate one.
double TriQuality(const vec3& a, const vec3& b, const vec3& c) {
  const double perimeter2 = length2(b - a) + length2(c - b) + length2(a - c);
  if (!(perimeter2 > 0)) return 0;
  return 2 * std::sqrt(3.0) * length(cross(b - a, c - a)) / perimeter2;
}

struct QualityStats {
  double min = std::numeric_limits<double>::infinity();
  double mean = 0;
  int below01 = 0;
};

QualityStats MeasureQuality(const Manifold& manifold) {
  const MeshGL64 mesh = manifold.GetMeshGL64();
  QualityStats stats;
  if (mesh.NumTri() == 0) {
    stats.min = 0;
    return stats;
  }
  for (size_t tri = 0; tri < mesh.NumTri(); ++tri) {
    const vec3 a = mesh.GetVertPos(mesh.triVerts[3 * tri]);
    const vec3 b = mesh.GetVertPos(mesh.triVerts[3 * tri + 1]);
    const vec3 c = mesh.GetVertPos(mesh.triVerts[3 * tri + 2]);
    const double quality = TriQuality(a, b, c);
    stats.min = std::min(stats.min, quality);
    stats.mean += quality;
    if (quality < 0.1) ++stats.below01;
  }
  stats.mean /= mesh.NumTri();
  return stats;
}

void ExpectValidRoundTrip(const Manifold& manifold) {
  ASSERT_EQ(manifold.Status(), Manifold::Error::NoError);
  ASSERT_EQ(Manifold(manifold.GetMeshGL64()).Status(),
            Manifold::Error::NoError);
  EXPECT_EQ(manifold.NumDegenerateTris(), 0);
}

}  // namespace

TEST(Quality, CollapsesShortEdgesAndImprovesShape) {
  // A unit square with a needle-thin triangular notch. Extrusion gives the caps
  // two very short edges and a pair of sliver triangles.
  const double width = 1e-5;
  const Polygons section = {{
      {0, 0},
      {0.5 - width, 0},
      {0.5, width},
      {0.5 + width, 0},
      {1, 0},
      {1, 1},
      {0, 1},
  }};
  const Manifold notched = Manifold::Extrude(section, 1);
  ASSERT_EQ(notched.Status(), Manifold::Error::NoError);
  const QualityStats before = MeasureQuality(notched);
  EXPECT_LT(before.min, 1e-3);
  EXPECT_GT(before.below01, 0);

  const Manifold improved = notched.ImproveQuality(1e-3);
  ExpectValidRoundTrip(improved);
  const QualityStats after = MeasureQuality(improved);

  EXPECT_LT(improved.NumTri(), notched.NumTri());
  EXPECT_LT(improved.NumVert(), notched.NumVert());
  EXPECT_GT(after.min, 0.5);
  EXPECT_GT(after.mean, before.mean);
  EXPECT_EQ(after.below01, 0);
  EXPECT_EQ(improved.Genus(), notched.Genus());
  EXPECT_NEAR(improved.Volume(), 1.0, 1e-5);
  EXPECT_NEAR(improved.SurfaceArea(), 6.0, 1e-5);

  // Cleanup is the documented alias.
  const Manifold cleaned = notched.Cleanup(1e-3);
  ExpectValidRoundTrip(cleaned);
  EXPECT_EQ(cleaned.NumTri(), improved.NumTri());
  EXPECT_NEAR(cleaned.Volume(), improved.Volume(), 1e-12);

  const Manifold removed = notched.RemoveDegenerates(1e-3);
  ExpectValidRoundTrip(removed);
  EXPECT_EQ(removed.NumTri(), improved.NumTri());
}

TEST(Quality, BooleanResultRetainsVolumeAndManifoldness) {
  const Manifold a = Manifold::Cube({1, 1, 1}, true);
  const Manifold b =
      Manifold::Cube({1, 1, 1}, true).Rotate(0, 0, 30).Translate(
          {0.2, 0.1, 0.05});
  const Manifold booleanResult = a + b;
  ASSERT_EQ(booleanResult.Status(), Manifold::Error::NoError);
  const QualityStats before = MeasureQuality(booleanResult);

  // All edges here are much longer than the tolerance, so this exercises the
  // coplanar flip pass rather than short-edge collapse.
  const Manifold improved = booleanResult.ImproveQuality(1e-3);
  ExpectValidRoundTrip(improved);
  const QualityStats after = MeasureQuality(improved);

  EXPECT_EQ(improved.NumTri(), booleanResult.NumTri());
  EXPECT_EQ(improved.NumVert(), booleanResult.NumVert());
  EXPECT_EQ(improved.Genus(), booleanResult.Genus());
  EXPECT_GT(after.mean, before.mean);
  EXPECT_GE(after.min, before.min - 1e-12);
  EXPECT_NEAR(improved.Volume(), booleanResult.Volume(), 1e-9);
  EXPECT_NEAR(improved.SurfaceArea(), booleanResult.SurfaceArea(), 1e-9);
}

TEST(Quality, SimplifyLeavesBetterTriangles) {
  const Manifold cube = Manifold::Cube({1, 1, 1}, true);
  const Manifold other =
      Manifold::Cube({1, 1, 1}, true).Rotate(0, 0, 30).Translate(
          {0.2, 0.1, 0.05});
  const Manifold result = cube + other;
  const QualityStats before = MeasureQuality(result);
  const Manifold simplified = result.Simplify(0.01);

  ExpectValidRoundTrip(simplified);
  EXPECT_GT(MeasureQuality(simplified).mean, before.mean);
  EXPECT_NEAR(simplified.Volume(), result.Volume(), 1e-6);
}

TEST(Quality, PropagatesErrors) {
  MeshGL input = TetGL();
  input.vertProperties[2 * 3 + 1] = NAN;
  const Manifold errored(input);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.ImproveQuality(0.01).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.Cleanup(0.01).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.RemoveDegenerates(0.01).Status(),
            Manifold::Error::NonFiniteVertex);
}

// The weekly benchmark suite includes this test. Its gtest timing is the
// quality-pass benchmark signal; the assertions keep the workload honest.
TEST(Quality, ImproveQualityBenchmark) {
  const Manifold sphere = Manifold::Sphere(1, 256);
  const Manifold cube = Manifold::Cube({1, 1, 1}, true)
                            .Rotate(17, 23, 31)
                            .Translate({0.2, 0.1, 0.05});
  const Manifold input = sphere + cube;
  ASSERT_EQ(input.Status(), Manifold::Error::NoError);
  ASSERT_GT(input.NumTri(), 30000);

  const auto start = std::chrono::steady_clock::now();
  const Manifold improved = input.ImproveQuality(0.01);
  const auto end = std::chrono::steady_clock::now();
  const double elapsedMs =
      std::chrono::duration<double, std::milli>(end - start).count();

  ExpectValidRoundTrip(improved);
  EXPECT_LE(improved.NumTri(), input.NumTri());
  EXPECT_GT(MeasureQuality(improved).mean, MeasureQuality(input).mean);
  EXPECT_NEAR(improved.Volume(), input.Volume(), 1e-5);
  RecordProperty("improve_quality_ms", static_cast<int>(elapsedMs));
}
