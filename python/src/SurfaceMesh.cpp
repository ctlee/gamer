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

using VertexID_Iterator = SurfaceMesh::SimplexIDIterator<1>;
using VertexData_Iterator = SurfaceMesh::DataIterator<1>;

using FaceID_Iterator = SurfaceMesh::SimplexIDIterator<3>;
using FaceData_Iterator = SurfaceMesh::DataIterator<3>;

void init_SurfaceMesh(py::module& mod){
    // Bindings for SurfaceMesh
    py::class_<SurfaceMesh> SurfMesh(mod, "SurfaceMesh",
        R"delim(
            Python wrapper around :cpp:type:`SurfaceMesh`.
        )delim"
    );
    SurfMesh.def(py::init<>(), "Default constructor");

    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::insert<1>),
        R"delim(
            Inserts a vertex based on key.

            Args:
                key (array): Array [1] of vertex key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 1>&, const Vertex&>(&SurfaceMesh::insert<1>),
        R"delim(
            Inserts a vertex and data based on key.

            Args:
                key (array): Array [1] of vertex key.
                data (:py:class:`Vertex`): Vertex data.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::insert<2>),
        R"delim(
            Insert an edge into the mesh

            Args:
                key (array): Array [2] of edge key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&, const Edge&>(&SurfaceMesh::insert<2>),
        R"delim(
            Insert an edge into the mesh

            Args:
                key (array): Array [2] of edge key.
                data (Edge): Edge data.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::insert<3>),
        R"delim(
            Insert a face into the mesh

            Args:
                key (array): Array [3] of face key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 3>&, const Face&>(&SurfaceMesh::insert<3>),
        R"delim(
            Insert a face into the mesh

            Args:
                arg1 (array): Array [3] of edge key.
                arg2 (Face): Face data.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::get_simplex_up<1>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [1] of vertex key.

            Returns:
                SimplexID (:py:class:`SMVertexID`): The object.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::get_simplex_up<2>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [2] of edge key.

            Returns:
                SimplexID (:py:class:`SMEdgeID`): The object.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::get_simplex_up<3>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [3] of vertex key.

            Returns:
                SimplexID (:py:class:`SMFaceID`): The object.
        )delim"
    );

    SurfMesh.def("print", &print, "Print a surface mesh");
    SurfMesh.def("__repr__", &print);


}