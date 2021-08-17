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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "gamer/SurfaceMesh.h"
#include "gamer/TetMesh.h"

/// Namespace for all things gamer
namespace gamer {
void writeComsol(const std::string &filename,
                 const std::vector<SurfaceMesh const *> &meshes) {
  std::ofstream fout(filename);
  if (!fout.is_open()) {
    std::cerr << "File '" << filename << "' could not be written to."
              << std::endl;
    // exit(1);
    return;
  }

  fout << "# Generated using GAMer\n\n";
  fout << "# Major & minor version\n";
  fout << "0 1\n";

  fout << "1 # number of tags\n";
  fout << "# Tags\n";
  fout << "5 mesh1\n";
  // for (std::size_t idx = 0; idx < meshes.size(); ++idx) {
  //   std::stringstream ss;
  //   ss << "mesh" << idx;
  //   fout << ss.str().length() << " " << ss.str() << "\n";
  // }

  fout << "1 # number of types\n";
  fout << "# Types\n";
  // for (std::size_t idx = 0; idx < meshes.size(); ++idx)
  fout << "3 obj\n";

  fout << "\n# --------- Object 0 ----------\n\n";
  fout << "0 0 1\n";
  fout << "4 Mesh # class\n";
  fout << "4 # version\n";
  fout << "3 # sdim\n";

  std::size_t nvertices = 0;
  std::size_t ntri = 0;

  for (std::size_t idx = 0; idx < meshes.size(); ++idx) {
    const SurfaceMesh &mesh = *meshes[idx];
    nvertices += mesh.size<1>();
    ntri += mesh.size<3>();
  }

  fout << nvertices << " # number of mesh vertices\n";
  fout << "0 # lowest mesh vertex index\n\n";

  std::vector<
      std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType>>
      sigma;
  typename SurfaceMesh::KeyType cnt = 0;

  fout << "# Mesh point coordinates\n";
  fout.precision(10);
  for (std::size_t idx = 0; idx < meshes.size(); ++idx) {
    const SurfaceMesh &mesh = *meshes[idx];
    sigma.emplace_back();

    // Get the vertex data directly
    for (const auto vertexID : mesh.get_level_id<1>()) {
      sigma[idx][mesh.get_name(vertexID)[0]] = cnt++;
      auto vertex = *vertexID;

      fout << vertex[0] << " " << vertex[1] << " " << vertex[2] << " "
           << "\n";
    }
  }
  fout << "\n";

  fout << "1 # number of element types\n\n";
  fout << "# Type #0\n\n";
  fout << "3 tri # type name\n\n";
  fout << "3 # number of vertices per element\n";
  fout << ntri << " # number of elements\n";
  fout << "# Elements\n";

  for (std::size_t idx = 0; idx < meshes.size(); ++idx) {
    const SurfaceMesh &mesh = *meshes[idx];
    bool orientationError = false;

    // Get the face nodes
    for (auto faceID : mesh.get_level_id<3>()) {
      auto w = mesh.get_name(faceID);

      auto orientation = (*faceID).orientation;
      if (orientation == 1) {
        fout << sigma[idx][w[0]] << " " << sigma[idx][w[1]] << " "
             << sigma[idx][w[2]] << "\n";
      } else if (orientation == -1) {
        fout << sigma[idx][w[2]] << " " << sigma[idx][w[1]] << " "
             << sigma[idx][w[0]] << "\n";

      } else {
        orientationError = true;
        fout << sigma[idx][w[0]] << " " << sigma[idx][w[1]] << " "
             << sigma[idx][w[2]] << "\n";
      }
    }
    if (orientationError) {
      std::cerr << "WARNING(writeComsol): The orientation of one or more faces "
                << "is not defined. Did you run compute_orientation()?"
                << std::endl;
    }
  }
  fout << "\n" << ntri << " # number of geometric entity indices\n";
  fout << "# Geometric entity indices\n";
  for (std::size_t idx = 0; idx < meshes.size(); ++idx) {
    const SurfaceMesh &mesh = *meshes[idx];
    for (auto faceID : mesh.get_level_id<3>()) {
      fout << (*faceID).marker << "\n";
    }
  }
  fout.close();
} // namespace gamer

void writeComsol(const std::string &filename, const SurfaceMesh &mesh) {
  std::vector<SurfaceMesh const *> v{&mesh};
  writeComsol(filename, v);
}

void writeComsol(const std::string &filename, const TetMesh &mesh) {

  if ((*mesh.get_simplex_up()).higher_order == true) {
    gamer_runtime_error(
        "Comsol output does not support higher order meshes at this point.");
  }

  std::ofstream fout(filename);
  if (!fout.is_open()) {
    std::stringstream ss;
    ss << "File '" << filename << "' could not be written to.";
    gamer_runtime_error(ss.str());
  }

  fout << "# Generated using GAMer\n\n";
  fout << "# Major & minor version\n";
  fout << "0 1\n";

  fout << "1 # number of tags\n";
  fout << "# Tags\n";
  fout << "5 mesh1\n";
  // for (std::size_t idx = 0; idx < meshes.size(); ++idx) {
  //   std::stringstream ss;
  //   ss << "mesh" << idx;
  //   fout << ss.str().length() << " " << ss.str() << "\n";
  // }

  fout << "1 # number of types\n";
  fout << "# Types\n";
  // for (std::size_t idx = 0; idx < meshes.size(); ++idx)
  fout << "3 obj\n";

  fout << "\n# --------- Object 0 ----------\n\n";
  fout << "0 0 1\n";
  fout << "4 Mesh # class\n";
  fout << "4 # version\n";
  fout << "3 # sdim\n";

  std::map<typename TetMesh::KeyType, typename TetMesh::KeyType> sigma;
  size_t cnt = 0;

  fout << mesh.size<1>() << " # number of mesh vertices\n";
  fout << "0 # lowest mesh vertex index\n\n";
  fout << "# Mesh point coordinates\n";

  fout.precision(10);
  for (const auto vertexID : mesh.get_level_id<1>()) {
    size_t idx = cnt;
    sigma[mesh.get_name(vertexID)[0]] = cnt++;
    auto vertex = *vertexID;
    fout << vertex[0] << " " << vertex[1] << " " << vertex[2] << " "
         << "\n";
  }
  fout << "\n";

  fout << "2 # number of element types\n\n";
  fout << "# Type #0\n\n";
  fout << "3 tet # type name\n\n";
  fout << "4 # number of vertices per element\n";
  fout << mesh.size<4>() << " # number of elements\n";
  fout << "# Elements\n";
  // Print out Tetrahedra
  cnt = 0;
  std::vector<std::tuple<std::size_t, std::size_t, int>> faceMarkerList;
  std::vector<std::tuple<std::size_t, int>> cellMarkerList;
  bool orientationError = false;

  for (const auto tetID : mesh.get_level_id<4>()) {
    std::size_t idx = cnt++;
    auto w = mesh.get_name(tetID);
    auto orientation = (*tetID).orientation;

    if (orientation == 1) {
      fout << std::setw(4) << sigma[w[0]] << " " << std::setw(4) << sigma[w[1]]
           << " " << std::setw(4) << sigma[w[2]] << " " << std::setw(4)
           << sigma[w[3]] << "\n";
    } else if (orientation == -1) {
      fout << std::setw(4) << sigma[w[3]] << " " << std::setw(4) << sigma[w[1]]
           << " " << std::setw(4) << sigma[w[2]] << " " << std::setw(4)
           << sigma[w[0]] << "\n";
    } else {
      orientationError = true;
      fout << std::setw(4) << sigma[w[0]] << " " << std::setw(4) << sigma[w[1]]
           << " " << std::setw(4) << sigma[w[2]] << " " << std::setw(4)
           << sigma[w[3]] << "\n";
    }
  }
  if (orientationError) {
    std::cerr << "WARNING(writeComsol): The orientation of one or more cells "
              << "is not defined. Did you run compute_orientation()?"
              << std::endl;
  }

  fout << "\n" << mesh.size<4>() << " # number of geometric entity indices\n";
  fout << "# Geometric entity indices\n";
  for (const auto tetID : mesh.get_level_id<4>()) {
    fout << (*tetID).marker << "\n";
  }

  fout << "\n# Type #1\n\n";
  fout << "3 tri # type name\n\n";
  fout << "3 # number of vertices per element\n";
  fout << mesh.size<3>() << " # number of elements\n";
  fout << "# Elements\n";
  orientationError = false;
  // Get the face nodes
  for (auto faceID : mesh.get_level_id<3>()) {
    auto w = mesh.get_name(faceID);
    fout << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";
  }
  if (orientationError) {
    std::cerr << "WARNING(writeComsol): The orientation of one or more faces "
              << "is not defined. Did you run compute_orientation()?"
              << std::endl;
  }

  fout << "\n" << mesh.size<3>() << " # number of geometric entity indices\n";
  fout << "# Geometric entity indices\n";
  for (const auto faceID : mesh.get_level_id<3>()) {
    fout << (*faceID).marker << "\n";
  }
  fout.close();
}

} // namespace gamer
