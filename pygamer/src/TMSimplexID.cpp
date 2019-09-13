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

#include "gamer/TetMesh.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

using TMVertexID = TetMesh::SimplexID<1>;
using TMEdgeID = TetMesh::SimplexID<2>;
using TMFaceID = TetMesh::SimplexID<3>;
using TMCellID = TetMesh::SimplexID<4>;

void init_TMSimplexID(py::module& mod){


    /************************************
     *  Binding TMVertexID
     ************************************/
    py::class_<TMVertexID> vid(mod, "VertexID",
        R"delim(
            Wrapper around :cpp:type:`gamer::TetMesh`::SimplexID<1> object.

            This is a token to represent a 1-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    vid.def(py::init<>());

    vid.def("data",
            py::overload_cast<>(&TMVertexID::data),
            py::return_value_policy::reference_internal,
            R"delim(
                Access the data stored on the 1-simplex.

                Returns:
                    :py:class:`Vertex`: Vertex data
            )delim"
        );
    vid.def("isValid",
        [](const TMVertexID lhs) {return lhs != nullptr;},
        R"delim(
            Checks if VertexID refers to a valid simplex.

            Returns:
                bool: True if valid.
        )delim"
    );
    vid.def("indices",
        &TMVertexID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );
    vid.def("cover",
        &TMVertexID::cover,
        R"delim(
            Get the coboundary relations of this simplex.

            Returns:
                list: Indices of coboundary relations
        )delim"
    );
    vid.def("__repr__",
        [](const TMVertexID vid){
            std::ostringstream out;
            out << vid;
            return out.str();
        }
    );

    /************************************
     *  Binding TMEdgeID
     ************************************/
    py::class_<TMEdgeID> eid(mod, "EdgeID",
        R"delim(
            Wrapper around :cpp:type:`gamer::TetMesh`::SimplexID<2> object

            This is a token to represent a 2-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    eid.def("data", py::overload_cast<>(&TMEdgeID::data),
        py::return_value_policy::reference_internal,
            R"delim(
                Access the data stored on the edge.

                Returns:
                    :py:class:`Edge`: Edge data
            )delim"
    );
    eid.def("isValid",
        [](const TMEdgeID lhs) {return lhs != nullptr;},
        R"delim(
            Checks if EdgeID refers to a valid simplex.

            Returns:
                bool: True if valid.
        )delim"
    );
    eid.def("indices",
        &TMEdgeID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );
    eid.def("cover",
        &TMEdgeID::cover,
        R"delim(
            Get the coboundary relations of this simplex.

            Returns:
                list: Indices of coboundary relations
        )delim"
    );
    eid.def("__repr__",
        [](const TMEdgeID eid){
            std::ostringstream out;
            out << eid;
            return out.str();
        }
    );


    /************************************
     *  Binding TMFaceID
     ************************************/
    py::class_<TMFaceID> fid(mod, "FaceID",
        R"delim(
            Wrapper around :cpp:type:`gamer::TetMesh`::SimplexID<3> object

            This is a token to represent a 3-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    fid.def("data",
        py::overload_cast<>(&TMFaceID::data),
        py::return_value_policy::reference_internal,
        R"delim(
            Access the data stored on the 3-simplex

            Returns:
                :py:class:`Face`: Face data
        )delim"
    );
    fid.def("isValid",
        [](const TMFaceID lhs) {return lhs != nullptr;},
        R"delim(
            Checks if FaceID refers to a valid simplex.

            Returns:
                bool: True if valid.
        )delim"
    );
    fid.def("indices",
        &TMFaceID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );
    fid.def("cover",
        &TMFaceID::cover,
        R"delim(
            Get the coboundary relations of this simplex.

            Returns:
                list: Indices of coboundary relations
        )delim"
    );
    fid.def("__repr__",
        [](const TMFaceID fid){
            std::ostringstream out;
            out << fid;
            return out.str();
        }
    );


    /************************************
     *  Binding TMCellID
     ************************************/
    py::class_<TMCellID> cid(mod, "CellID",
        R"delim(
            Wrapper around :cpp:type:`gamer::TetMesh`::SimplexID<4> object

            This is a token to represent a 4-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    cid.def("data",
        py::overload_cast<>(&TMCellID::data),
        py::return_value_policy::reference_internal,
        R"delim(
            Access the data stored on the 4-simplex

            Returns:
                :py:class:`Cell`: Cell data
        )delim"
    );
    cid.def("isValid",
        [](const TMCellID lhs) {return lhs != nullptr;},
        R"delim(
            Checks if CellID refers to a valid simplex.

            Returns:
                bool: True if valid.
        )delim"
    );
    cid.def("indices",
        &TMCellID::indices,
        R"delim(
            Get the indices of this simplex.

            Returns:
                list: Indices
        )delim"
    );
    cid.def("__repr__",
        [](const TMCellID cid){
            std::ostringstream out;
            out << cid;
            return out.str();
        }
    );
}

} // end namespace gamer