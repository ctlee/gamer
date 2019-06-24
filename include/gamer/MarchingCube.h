/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2018
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

/**
 * @file MarchingCube.h
 * @brief Implementation of marching cubes algorithm
 */

#pragma once

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <queue>
#include "gamer/gamer.h"
#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{

/// Marching cubes edge table
static const int edgeTable[256] = {
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

/// Marching cubes triangle table
static const int triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};


/**
 * @brief      Marching cubes algorithm
 *
 * @param      dataset    Voxel array to mesh
 * @param[in]  maxval     Maximum value in the dataset
 * @param[in]  dim        Dimension of the dataset
 * @param[in]  span       Real space size of a voxel
 * @param[in]  isovalue   Isovalue to contour at
 * @param[in]  holelist   Inserter to append holes
 *
 * @tparam     NumType    Numerical typename
 * @tparam     <unnamed>  Check to ensure NumType is numerical
 * @tparam     Inserter   Typename of the inserter
 *
 * @return     Surface mesh
 */
template <typename NumType, typename = std::enable_if_t<std::is_arithmetic<NumType>::value>,
          class Inserter>
std::unique_ptr<SurfaceMesh> marchingCubes(
    NumType       * dataset,
    NumType         maxval,
    const i3Vector &dim,
    const f3Vector &span,
    NumType         isovalue,
    Inserter        holelist
    )
{
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);
    bool* mask = new bool[dim[0]*dim[1]*dim[2]];
    for (int i = 0; i < dim[0]*dim[1]*dim[2]; ++i){
        mask[i] = false;
    }

    std::cout << "Isolating isosurface" << std::endl;

    // TODO: (4) is it necessary to go through this three step masking process?

    // Mask all of the vertices connected to {0,0,0} outside of isosurface
    std::deque<i3Vector> visit = {i3Vector({0, 0, 0})};
    while (!visit.empty())
    {
        const i3Vector tmp = visit.front();
        visit.pop_front();

        // Look at the neighbor vertices
        for (int k = std::max(tmp[2]-1, 0); k <= std::min(tmp[2]+1, dim[2]-1); ++k)
        {
            for (int j = std::max(tmp[1]-1, 0); j <= std::min(tmp[1]+1, dim[1]-1); ++j)
            {
                for (int i = std::max(tmp[0]-1, 0); i <= std::min(tmp[0]+1, dim[0]-1); ++i)
                {
                    if ((dataset[Vect2Index(i, j, k, dim)] < isovalue)
                        && !mask[Vect2Index(i, j, k, dim)])
                    {
                        mask[Vect2Index(i, j, k, dim)] = true;
                        visit.push_back(i3Vector({i, j, k}));
                    }
                }
            }
        }
    }

    // Find internal holes
    for (int n = 0; n < dim[2]; n++)
    {
        for (int m = 0; m < dim[1]; m++)
        {
            for (int l = 0; l < dim[0]; l++)
            {
                if ((dataset[Vect2Index(l, m, n, dim)] < isovalue)
                    && !mask[Vect2Index(l, m, n, dim)])
                {
                    int holesize = 1;
                    visit.push_back(i3Vector({l, m, n}));

                    std::vector<int> holevoxels;

                    while (!visit.empty())
                    {
                        const i3Vector tmp = visit.front();
                        visit.pop_front();

                        // Look at the neighbor vertices
                        for (int i = std::max(tmp[0]-1, 0); i <= std::min(tmp[0]+1, dim[0]-1); ++i)
                        {
                            for (int j = std::max(tmp[1]-1, 0); j <= std::min(tmp[1]+1, dim[1]-1); ++j)
                            {
                                for (int k = std::max(tmp[2]-1, 0); k <= std::min(tmp[2]+1, dim[2]-1); ++k)
                                {
                                    if ((dataset[Vect2Index(i, j, k, dim)] < isovalue)
                                        && !mask[Vect2Index(i, j, k, dim)])
                                    {
                                        int idx = Vect2Index(i, j, k, dim);
                                        holevoxels.push_back(idx);

                                        mask[idx] = true;
                                        visit.push_back(i3Vector({i, j, k}));
                                        holesize++;
                                    }
                                }
                            }
                        }
                    }
                    std::cout << "Hole size: " << holesize << std::endl;

                    if (holesize < MIN_VOLUME)
                    {
                        for (auto idx : holevoxels)
                        {
                            dataset[idx] = maxval;
                            mask[idx] = false;
                        }
                    }
                    else
                    {
                        Vector v = Vector({static_cast<double>(l),
                                           static_cast<double>(m),
                                           static_cast<double>(n)}).ElementwiseProduct(span);
                        *holelist++ = v;
                        std::cout << "Hole real size: " << v << std::endl;
                    }
                }
            }
        }
    }
    std::cout << "Done isolating isosurface" << std::endl;

    size_t              vertexNum = 0;
    size_t              triNum = 0;
    std::array<int, 3>* triangles = new std::array<int, 3>[dim[0]*dim[1]*dim[2]];
    i3Vector          * edges = new i3Vector[dim[0]*dim[1]*dim[2]];
    Vector            * vertices = new Vector[dim[0]*dim[1]*dim[2]];

    // This section in particular is weird...
    for (int k = 0; k < dim[2]-1; k++)
    {
        for (int j = 0; j < dim[1]-1; j++)
        {
            for (int i = 0; i < dim[0]-1; i++)
            {
                int idx = Vect2Index(i, j, k, dim);
                // If isovalue is within tolerance make it bigger
                if ((dataset[idx] > isovalue - 0.0001) && (dataset[idx] < isovalue + 0.0001))
                    dataset[idx] = isovalue + 0.0001;
                edges[idx] = i3Vector({-1, -1, -1});

                if (dataset[idx] >= isovalue)
                {
                    mask[idx] = false;
                }
                else
                {
                    mask[idx] = true;
                }
            }
        }
    }

    std::cout << "Marching..." << std::endl;
    // Marching cubes vertex indices and edges convention
    //		   v4_________e4_________v5
    //			/|                  /|
    //		e7 / |                 / |
    //		  /  |             e5 /  |
    //		 /   | e8            /	 | e9
    //	  v7/____|_____e6_______/v6	 |
    //		|	 |              |	 |
    //	    |  v0|______e0______|____|v1
    //	e11 |	/               |   /
    //		|  /			e10	|  /
    //		| /	e3				| / e1
    //		|/					|/
    //	  v3/_________e2________/v2
    //
    for (int i = 0; i < dim[0]-1; i++)
    {
        for (int j = 0; j < dim[1]-1; j++)
        {
            for (int k = 0; k < dim[2]-1; k++)
            {

                // List of vertices for the current cell
                int cellVertices[12];
                std::fill_n(cellVertices, 12, -1);

                // Table of integer indices to overall dataset array
                int indexTable[8];
                indexTable[0] = Vect2Index(i, j, k, dim);
                indexTable[1] = Vect2Index(i, j+1, k, dim);
                indexTable[2] = Vect2Index(i+1, j+1, k, dim);
                indexTable[3] = Vect2Index(i+1, j, k, dim);
                indexTable[4] = Vect2Index(i, j, k+1, dim);
                indexTable[5] = Vect2Index(i, j+1, k+1, dim);
                indexTable[6] = Vect2Index(i+1, j+1, k+1, dim);
                indexTable[7] = Vect2Index(i+1, j, k+1, dim);

                int cellIndex = 0;  // Bitmask for intersections
                for (int idx = 0; idx < 8; ++idx)
                {
                    if (mask[indexTable[idx]])
                    {
                        cellIndex |= (1 << idx);
                    }
                }
                // std::cout << std::bitset<8>(cellIndex) << std::endl;

                NumType den1, den2;

                if (edgeTable[cellIndex] & (1 << 0))
                {
                    auto &edgeIdx = edges[indexTable[0]][1];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[0]];
                        den2 = dataset[indexTable[1]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i),
                                                      static_cast<double>(j) + ratio,
                                                      static_cast<double>(k)
                                                     }).ElementwiseProduct(span);
                        cellVertices[0] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[0] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 1))
                {
                    auto &edgeIdx = edges[indexTable[1]][0];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[1]];
                        den2 = dataset[indexTable[2]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + ratio,
                                                      static_cast<double>(j) + 1,
                                                      static_cast<double>(k)
                                                     }).ElementwiseProduct(span);
                        cellVertices[1] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[1] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 2))
                {
                    auto &edgeIdx = edges[indexTable[3]][1];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[3]];
                        den2 = dataset[indexTable[2]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + 1,
                                                      static_cast<double>(j) + ratio,
                                                      static_cast<double>(k)
                                                     }).ElementwiseProduct(span);
                        cellVertices[2] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[2] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 3))
                {
                    auto &edgeIdx = edges[indexTable[0]][0];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[0]];
                        den2 = dataset[indexTable[3]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + ratio,
                                                      static_cast<double>(j),
                                                      static_cast<double>(k)
                                                     }).ElementwiseProduct(span);
                        cellVertices[3] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[3] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 4))
                {
                    auto &edgeIdx = edges[indexTable[4]][1];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[4]];
                        den2 = dataset[indexTable[5]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i),
                                                      static_cast<double>(j) + ratio,
                                                      static_cast<double>(k) + 1
                                                     }).ElementwiseProduct(span);
                        cellVertices[4] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[4] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 5))
                {
                    auto &edgeIdx = edges[indexTable[5]][0];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[5]];
                        den2 = dataset[indexTable[6]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + ratio,
                                                      static_cast<double>(j) + 1,
                                                      static_cast<double>(k) + 1
                                                     }).ElementwiseProduct(span);
                        cellVertices[5] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[5] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 6))
                {
                    auto &edgeIdx = edges[indexTable[7]][1];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[7]];
                        den2 = dataset[indexTable[6]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + 1,
                                                      static_cast<double>(j) + ratio,
                                                      static_cast<double>(k) + 1
                                                     }).ElementwiseProduct(span);
                        cellVertices[6] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[6] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 7))
                {
                    auto &edgeIdx = edges[indexTable[4]][0];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[4]];
                        den2 = dataset[indexTable[7]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + ratio,
                                                      static_cast<double>(j),
                                                      static_cast<double>(k) + 1
                                                     }).ElementwiseProduct(span);
                        cellVertices[7] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[7] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 8))
                {
                    auto &edgeIdx = edges[indexTable[0]][2];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[0]];
                        den2 = dataset[indexTable[4]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i),
                                                      static_cast<double>(j),
                                                      static_cast<double>(k) + ratio
                                                     }).ElementwiseProduct(span);
                        cellVertices[8] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[8] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 9))
                {
                    auto &edgeIdx = edges[indexTable[1]][2];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[1]];
                        den2 = dataset[indexTable[5]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i),
                                                      static_cast<double>(j) + 1,
                                                      static_cast<double>(k) + ratio
                                                     }).ElementwiseProduct(span);
                        cellVertices[9] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[9] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 10))
                {
                    auto &edgeIdx = edges[indexTable[2]][2];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[2]];
                        den2 = dataset[indexTable[6]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + 1,
                                                      static_cast<double>(j) + 1,
                                                      static_cast<double>(k) + ratio
                                                     }).ElementwiseProduct(span);
                        cellVertices[10] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[10] = edgeIdx;
                    }
                }

                if (edgeTable[cellIndex] & (1 << 11))
                {
                    auto &edgeIdx = edges[indexTable[3]][2];
                    if (edgeIdx == -1)
                    {
                        den1 = dataset[indexTable[3]];
                        den2 = dataset[indexTable[7]];
                        NumType ratio = (den1 != den2) ? (isovalue-den1)/(den2-den1) : 0;
                        // std::cout << den1 << " " << den2 << " " << ratio <<
                        // std::endl;
                        vertices[vertexNum] = Vector({static_cast<double>(i) + 1,
                                                      static_cast<double>(j),
                                                      static_cast<double>(k) + ratio
                                                     }).ElementwiseProduct(span);
                        cellVertices[11] = vertexNum;
                        edgeIdx = vertexNum;
                        vertexNum++;
                    }
                    else
                    {
                        cellVertices[11] = edgeIdx;
                    }
                }

                int ii = 0;
                while (triTable[cellIndex][ii] != -1)
                {
                    triangles[triNum] = {cellVertices[triTable[cellIndex][ii++]],
                                         cellVertices[triTable[cellIndex][ii++]],
                                         cellVertices[triTable[cellIndex][ii++]]};
                    triNum++;
                }
            }
        }
    }

    for (int i = 0; i < vertexNum; ++i)
    {
        mesh->insert<1>({i}, SMVertex(vertices[i]));
    }
    for (int i = 0; i < triNum; ++i)
    {
        mesh->insert<3>(triangles[i]);
    }

    delete[] vertices;
    delete[] triangles;
    delete[] mask;
    delete[] edges;

    compute_orientation(*mesh);
    return mesh;
}
} // end namespace gamer
