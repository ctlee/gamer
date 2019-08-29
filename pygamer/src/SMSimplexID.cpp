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

#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

using SMVertexID = SurfaceMesh::SimplexID<1>;
using SMEdgeID = SurfaceMesh::SimplexID<2>;
using SMFaceID = SurfaceMesh::SimplexID<3>;

void init_SMSimplexID(py::module& mod){
    // Bindings for SMVertexID
    py::class_<SMVertexID> vid(mod, "VertexID",
        R"delim(
            Wrapper around :cpp:type:`gamer::SurfaceMesh`::SimplexID<1> object.

            This is a token to represent a 1-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    vid.def(py::init<>());


    vid.def("data",
            py::overload_cast<>(&SMVertexID::data),
            py::return_value_policy::reference_internal,
            R"delim(
                Access the data stored on the 1-simplex.

                Returns:
                    :py:class:`surfacemesh.Vertex`: Vertex data
            )delim"
        );


    vid.def("isValid",
        [](const SMVertexID& lhs) {return lhs != nullptr;},
        R"delim(
            Checks if VertexID refers to a valid simplex.

            This is useful for validating if :py:func:`get_simplex_up` has
            returned a valid SimplexID.

            Returns:
                bool: True if valid.
        )delim"
    );


    vid.def("indices",
        &SMVertexID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );


    vid.def("cover",
        &SMVertexID::cover,
        R"delim(
            Get the coboundary relations of this simplex.

            Returns:
                list: Indices of coboundary relations
        )delim"
    );
    vid.def("__repr__",
        [](const SMVertexID vid){
            std::ostringstream out;
            out << vid;
            return out.str();
        }
    );


    // Bindings for SMEdgeID
    py::class_<SMEdgeID> eid(mod, "EdgeID",
        R"delim(
            Wrapper around :cpp:type:`gamer::SurfaceMesh`::SimplexID<2> object

            This is a token to represent a 2-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    eid.def_property_readonly("dataro",
        py::overload_cast<>(&SMEdgeID::data),
        R"delim(
            Access the data stored on the edge.

            Returns:
                :py:class:`Edge`: Edge data
        )delim"
    );
    eid.def("data", py::overload_cast<>(&SMEdgeID::data),
        py::return_value_policy::reference_internal,
            R"delim(
                Access the data stored on the edge.

                Returns:
                    :py:class:`Edge`: Edge data
            )delim"
    );
    eid.def("isValid",
        [](const SMEdgeID& lhs) {return lhs != nullptr;},
        R"delim(
            Checks if EdgeID refers to a valid simplex.

            Returns:
                bool: True if valid.
        )delim"
    );
    eid.def("indices",
        &SMEdgeID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );
    eid.def("cover",
        &SMEdgeID::cover,
        R"delim(
            Get the coboundary relations of this simplex.

            Returns:
                list: Indices of coboundary relations
        )delim"
    );
    eid.def("__repr__",
        [](const SMEdgeID eid){
            std::ostringstream out;
            out << eid;
            return out.str();
        }
    );


    // Bindings for SMFaceID
    py::class_<SMFaceID> fid(mod, "FaceID",
        R"delim(
            Wrapper around :cpp:type:`gamer::SurfaceMesh`::SimplexID<3> object

            This is a token to represent a 3-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    fid.def_property_readonly("dataro",
        py::overload_cast<>(&SMFaceID::data),
        R"delim(
            Access the data stored on the 3-simplex

            Returns:
                :py:class:`Face`: Face data
        )delim"
    );
    fid.def("data",
        py::overload_cast<>(&SMFaceID::data),
        py::return_value_policy::reference_internal,
        R"delim(
            Access the data stored on the 3-simplex

            Returns:
                :py:class:`Face`: Face data
        )delim"
    );
    fid.def("isValid",
        [](const SMFaceID& lhs) {return lhs != nullptr;},
        R"delim(
            Checks if FaceID refers to a valid simplex.

            Returns:
                bool: True if valid.
        )delim"
    );
    fid.def("indices",
        &SMFaceID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );
    fid.def("__repr__",
        [](const SMFaceID fid){
            std::ostringstream out;
            out << fid;
            return out.str();
        }
    );
}

} // end namespace gamer