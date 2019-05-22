/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2019
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SurfaceMesh.h"

namespace py = pybind11;

// Forward function declarations
void init_SurfaceMesh(py::module &);
void init_SMSimplexID(py::module &);
void init_SMVertex(py::module &);
void init_SMFace(py::module &);
void init_SMFunctions(py::module &);


// Initialize the main `pygamer` module
PYBIND11_MODULE(pygamer, pygamer) {
    pygamer.doc() = "Python wrapper around the GAMer C++ library.";

    // pygamer.surfacemesh submodule defs
    py::module SurfMeshMod = pygamer.def_submodule("surfacemesh",
        "Submodule containing SurfaceMesh along with related objects and methods.");

    init_SurfaceMesh(SurfMeshMod);  // SurfaceMesh class
    init_SMVertex(SurfMeshMod);     // Vertex class
    init_SMFace(SurfMeshMod);       // Face class
    init_SMSimplexID(SurfMeshMod);  // SurfaceMesh::SimplexID class
    // Functions applied to pygamer.SurfaceMesh.SurfaceMesh
    init_SMFunctions(SurfMeshMod);
}