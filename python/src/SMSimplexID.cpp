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

using VertexID = SurfaceMesh::SimplexID<1>;
using EdgeID = SurfaceMesh::SimplexID<2>;
using FaceID = SurfaceMesh::SimplexID<3>;

void init_SMSimplexID(py::module& mod){
    // Bindings for VertexID
    py::class_<VertexID> vid(mod, "VertexID",
        R"delim(
            Wrapper around :cpp:type:`SurfaceMesh`::SimplexID<1> object. :py:

            This is a token to represent a 1-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    vid.def("data",
            py::overload_cast<>(&VertexID::data),
            py::return_value_policy::reference_internal,
            R"delim(
                Access the data stored on the 1-simplex

                Returns:
                    Vertex (:py:class:`Vertex`): Vertex data
            )delim"
        );

    // Bindings for EdgeID
    py::class_<EdgeID> eid(mod, "EdgeID",
        R"delim(
            Wrapper around :cpp:type:`SurfaceMesh`::SimplexID<2> object

            This is a token to represent a 2-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    eid.def("data", py::overload_cast<>(&EdgeID::data),
        py::return_value_policy::reference_internal,
        R"delim(
            Access the data stored on the edge.
        )delim"
    );

    // Bindings for FaceID
    py::class_<FaceID> fid(mod, "FaceID",
        R"delim(
            Wrapper around :cpp:type:`SurfaceMesh`::SimplexID<3> object

            This is a token to represent a 3-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    fid.def("data",
        py::overload_cast<>(&FaceID::data),
        py::return_value_policy::reference_internal,
        R"delim(
            Access the data stored on the 3-simplex

            Returns:
                Face (:py:class:`Face`): Face data
        )delim"
    );
}