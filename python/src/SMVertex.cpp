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

#include "Vertex.h"

namespace py = pybind11;

void init_SMVertex(py::module& mod){
    py::class_<Vertex> vertex(mod, "Vertex",
        R"delim(
            Wrapper around a :cpp:class:`Vertex`.
        )delim"
    );
    vertex.def(py::init<>(), "Default constructor.");
    vertex.def(py::init<double, double, double>(),
        py::arg("x"), py::arg("y"), py::arg("z"),
        "Constructor defining coordinates.");
    vertex.def(py::init<double, double, double, int, bool>(),
        py::arg("x"), py::arg("y"), py::arg("z"), py::arg("marker"), py::arg("selected"),
        "Constructor defining doordinates, marker, and selection status.");

    vertex.def("__getitem__",
        [](const Vertex &v, std::size_t i) -> const double& {
            if (i >= 3) throw py::index_error("Vector only contains three coordinates");
            return v[i];
        },
        "Get coordinate of position");
    vertex.def("__setitem__",
        [](Vertex &v, size_t i, double val) {
            if (i >= 3) throw py::index_error("Vector only contains three coordinates");
            v[i] = val;
        },
        "Set coordinate...");
    vertex.def_readwrite("marker", &Vertex::marker, "Boundary marker value");
    vertex.def_readwrite("selected", &Vertex::selected, "Selection status of vertex");
    vertex.def("__repr__", &Vertex::to_string);
}