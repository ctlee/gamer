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

/**
 * @file gamer.h
 * @brief Contains various global definitions and type definitions used in
 *        the overall project.
 */

#pragma once

#include "gamer/tensor.h"
#include <Eigen/Dense>
#include <sstream>

namespace gamer {
namespace detail {
template <typename... T>
void throw_runtime_error(const char *function, const char *file, const int line,
                         T &&...ts) {
  std::stringstream ss;
  ss << "Error: ";
  // https://stackoverflow.com/questions/25680461/variadic-template-pack-expansion/25683817#25683817
  int dummy[] = {0, ((ss << std::forward<T>(ts) << " "), 0)...};
  static_cast<void>(dummy); // Avoid warning for unused variable
  ss << " in function " << function << " at " << file << ":" << line;
  throw std::runtime_error(ss.str());
}
} // namespace detail
} // namespace gamer

#ifdef _MSC_VER 
  #define __PRETTY_FUNCTION__ __FUNCSIG__ 
#endif

#define gamer_runtime_error(...)                                               \
  gamer::detail::throw_runtime_error(__PRETTY_FUNCTION__, __FILE__, __LINE__,  \
                                     __VA_ARGS__);

/// Namespace for all things gamer
namespace gamer {
/// Blurring blobbyness to use in conversion from PDB/PQR to 3D volumes
#define BLOBBYNESS -0.2f

/// Discretization rate of 3D volumes
#define DIM_SCALE 1.99

/// The minimal volume (in voxels) of islands to be automatically removed
#define MIN_VOLUME 333333

#ifdef SINGLE
/// Defines REAL to be float
#define REAL float
#else
/// Defines REAL to be double
#define REAL double
#endif

/// Floating point vector with precision defined by user at compile time
using Vector = tensor<REAL, 3, 1>;
/// 3 Vector of doubles
using Vector3d = tensor<double, 3, 1>;
/// 3 Vector of floats
using Vector3f = tensor<float, 3, 1>;
/// 3 Vector of integers
using Vector3i = tensor<int, 3, 1>;
/// 3 Vector of std::size_t
using Vector3szt = tensor<std::size_t, 3, 1>;

/// Eigen 3x3 Matrix of type REAL
using EigenMatrix = Eigen::Matrix<REAL, 3, 3>;
/// Eigen 3 Vector of type REAL
using EigenVector = Eigen::Matrix<REAL, 3, 1>;
/// Dynamically sized Eigen Matrix
using EigenMatrixN = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>;
/// Dynamically sized Eigen Vector
using EigenVectorN = Eigen::Matrix<REAL, Eigen::Dynamic, 1>;

/**
 * @brief      Convert 3D array indices to the index of a flat
 *             array.
 *
 * As a convention use the following loop structure to optimize
 * cache efficiency:
 * ~~~~~~~~~~~~~~~{.cpp}
 *  for (int k;...; k++){
 *    for (int j;...; j++){
 *      for (int i;...; i++){
 *  }}}
 * ~~~~~~~~~~~~~~~
 * In short, index i should occupy the inner most scope followed
 * by j then k.
 *
 * @param[in]  i     Value of the first index
 * @param[in]  j     Value of the second index
 * @param[in]  k     Value of the third index
 * @param[in]  dim   Dimensions of the 3D array
 *
 * @return     Index of flat array corresponding to indices in
 *             3D array.
 */
inline std::size_t Vect2Index(const std::size_t i, const std::size_t j,
                              const std::size_t k, const Vector3szt &dim) {
  return k * dim[0] * dim[1] + j * dim[0] + i;
}
} // end namespace gamer
