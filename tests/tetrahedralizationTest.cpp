// This file is part of the GAMer software.
// Copyright (C) 2016-2021
// by Christopher T. Lee and contributors
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, see <http://www.gnu.org/licenses/>
// or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
// Boston, MA 02111-1307 USA

#include "gamer/SurfaceMesh.h"
#include "gamer/TetMesh.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

/// Namespace for all things gamer
namespace gamer {

class TetrahedralizationTest : public testing::Test {
protected:
  TetrahedralizationTest() {}
  ~TetrahedralizationTest() {}
  virtual void SetUp() {
    outermesh = sphere(0);
    scale(*outermesh, outerScaleFactor);

    auto &omGlobal = outermesh->get_simplex_up().data();
    omGlobal.ishole = false;
    omGlobal.marker = volumeMarker;
    for (auto &fdata : outermesh->get_level<3>()) {
      fdata.marker = outerSurfaceMarker;
    }

    for (auto &vdata : outermesh->get_level<1>()) {
      EXPECT_TRUE((vdata.position | vdata.position) -
                      (outerScaleFactor * outerScaleFactor) <
                  tolerance);
    }

    innermesh = sphere(0);
    scale(*innermesh, innerScaleFactor);
    auto &imGlobal = innermesh->get_simplex_up().data();
    imGlobal.ishole = true;

    for (auto &fdata : innermesh->get_level<3>()) {
      fdata.marker = innerSurfaceMarker;
    }

    for (auto &vdata : innermesh->get_level<1>()) {
      EXPECT_TRUE((vdata.position | vdata.position) -
                      (innerScaleFactor * innerScaleFactor) <
                  tolerance);
    }
  }
  virtual void TearDown() {}

  std::unique_ptr<SurfaceMesh> outermesh;
  std::unique_ptr<SurfaceMesh> innermesh;
  int volumeMarker = 25;
  int outerSurfaceMarker = 2;
  int innerSurfaceMarker = 7;
  double innerScaleFactor = 2;
  double outerScaleFactor = 7;
  double tolerance = 1e-3;
};

TEST_F(TetrahedralizationTest, tetrahedralize) {
  std::vector<SurfaceMesh const *> meshes{outermesh.get(), innermesh.get()};
  auto tetmesh = makeTetMesh(meshes, "q1.3/10a1O8/7AYCQ");
  EXPECT_TRUE(tetmesh->size<4>() > 0);
  EXPECT_TRUE(tetmesh->size<3>() > outermesh->size<3>() + innermesh->size<3>());
  EXPECT_TRUE(tetmesh->size<2>() > outermesh->size<2>() + innermesh->size<2>());
  EXPECT_TRUE(tetmesh->size<1>() > outermesh->size<1>() + innermesh->size<1>());

  for (auto &tetData : tetmesh->get_level<4>()) {
    EXPECT_EQ(tetData.marker, volumeMarker);
  }

  for (auto fID : tetmesh->get_level_id<3>()) {
    if (tetmesh->onBoundary(fID)) {
      EXPECT_TRUE(fID.data().marker == outerSurfaceMarker ||
                  fID.data().marker == innerSurfaceMarker);
    }
  }

  for (auto &vData : tetmesh->get_level<1>()) {
    double dist = vData.position | vData.position;
    EXPECT_TRUE(dist - (innerScaleFactor * innerScaleFactor) > -tolerance);
    EXPECT_TRUE(dist - (outerScaleFactor * outerScaleFactor) < tolerance);
  }
}

} // end namespace gamer
