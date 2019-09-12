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

void init_TMVertex(py::module& mod){
    py::class_<TMVertex> vertex(mod, "Vertex",
        R"delim(
            Wrapper around a :cpp:class:`gamer::TMVertex`.
        )delim"
    );

    vertex.def(py::init<double, double, double, int, bool>(),
        py::arg("x") = 0,
        py::arg("y") = 0,
        py::arg("z") = 0,
        py::arg("marker") = -1,
        py::arg("selected") = false,
        "Constructor defining coordinates, marker, and selection status."
    );


    vertex.def("__getitem__",
        [](const TMVertex &v, std::size_t i) -> const double& {
            if (i >= 3) throw py::index_error("Vector only contains three coordinates");
            return v[i];
        },
        "Get value of coordinate position"
    );
    vertex.def("__setitem__",
        [](TMVertex &v, size_t i, double val) {
            if (i >= 3) throw py::index_error("Vector only contains three coordinates");
            v[i] = val;
        },
        "Set coordinate value"
    );

    vertex.def_readwrite("position", &TMVertex::position, "Position of the vertex");
    vertex.def_readwrite("marker", &TMVertex::marker, "Boundary marker value");
    vertex.def_readwrite("selected", &TMVertex::selected, "Selection status of vertex");
    vertex.def_readwrite("error", &TMVertex::error, "Error value");
    vertex.def("__repr__", &TMVertex::to_string);
}

} // end namespace gamer