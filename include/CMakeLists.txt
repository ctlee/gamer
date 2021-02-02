# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2021
# by Christopher T. Lee and contributors

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see <http://www.gnu.org/licenses/>
# or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# ***************************************************************************


set(GAMER_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/gamer"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/SurfaceMesh.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/TetMesh.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/EigenDiagonalization.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/MarchingCube.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/OsculatingJets.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/PDBReader.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/Vertex.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/gamer.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/stringutil.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/tensor.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/version.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/gamer/visibility_utils.h"
PARENT_SCOPE)
