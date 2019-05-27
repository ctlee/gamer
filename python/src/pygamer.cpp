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
void init_SMGlobal(py::module &);
void init_SMVertex(py::module &);
void init_SMEdge(py::module &);
void init_SMFace(py::module &);
void init_SMSimplexID(py::module &);
void init_SMFunctions(py::module &);
void init_SurfaceMesh(py::module &);

// Initialize the main `pygamer` module
PYBIND11_MODULE(pygamer, pygamer) {
    pygamer.doc() = "Python wrapper around the GAMer C++ library.";

    // pygamer.surfacemesh submodule defs
    py::module SurfMeshMod = pygamer.def_submodule("surfacemesh",
        "Submodule containing SurfaceMesh along with related objects and methods.");

    // WARNING: the order of initialization matters for dependent calls.
    init_SMGlobal(SurfMeshMod);     // Global class
    init_SMVertex(SurfMeshMod);     // Vertex class
    init_SMEdge(SurfMeshMod);       // Edge class
    init_SMFace(SurfMeshMod);       // Face class
    init_SMSimplexID(SurfMeshMod);  // SurfaceMesh::SimplexID class
    init_SMFunctions(SurfMeshMod);  // Functions pyg.sm.SM.*
    init_SurfaceMesh(SurfMeshMod);  // SurfaceMesh class


    pygamer.def("readOFF", &readOFF,
        py::arg("filename"),
        R"delim(
            Read OFF file to mesh

            Args:
                filename (str): Filename to read from
            Returns:
                :py:class:`SurfaceMesh`: Mesh of interest
        )delim"
    );


    pygamer.def("writeOFF", &writeOFF,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in OFF format

            Args:
                filename (str): Filename to write to
                mesh (:py:class:`SurfaceMesh`): Mesh of interest
        )delim"
    );


    pygamer.def("readOBJ", &readOBJ,
        py::arg("filename"),
        R"delim(
            Read OBJ file to mesh

            Args:
                filename (str): Filename to read from
            Returns:
                :py:class:`SurfaceMesh`: Mesh
        )delim"
    );


    pygamer.def("writeOBJ", &writeOBJ,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in OBJ format

            Args:
                filename (str): Filename to write to
                mesh (:py:class:`SurfaceMesh`): Mesh of interest
        )delim"
    );

}