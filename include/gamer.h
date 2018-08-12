/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2017
 * by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
 *    and Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * ***************************************************************************
 */

#pragma once

//#include <libraries/triangle/triangle.h>
//#define TRILIBRARY
#include "tensor.h"

// The number of rings to use to compute local structure tensor
#define RINGS 1

/** @brief Blurring blobyness used in conversion from PDB/PQR to 3D volumes */
#define BLOBBYNESS        -0.2f

/** @brief Discretization rate of 3D volumes */
#define DIM_SCALE         1.99

/** @brief The minimal volumes (in voxels) of islands to be removed */
#define MIN_VOLUME        333333


// // Other definitions and data structures
// /** @brief Other definition */
// #define _LITTLE_ENDIAN   1


using Vector = tensor<double,3,1>;
using f3Vector = tensor<float,3,1>;
using i3Vector = tensor<int,3,1>;

template <class T>
std::size_t Vect2Index(const T i, const T j, const T k, const i3Vector& dim){
    return k*dim[0]*dim[1] + j*dim[0] + i;
}