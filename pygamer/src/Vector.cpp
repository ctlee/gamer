// This file is part of the GAMer software.
// Copyright (C) 2016-2021
// by Christopher T. Lee and contributors
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, see <http://www.gnu.org/licenses/>
// or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
// Boston, MA 02111-1307 USA

#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gamer/gamer.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

void init_Vector(py::module& mod){
    py::class_<Vector> vector(mod, "Vector", py::buffer_protocol(),
        R"delim(
            Wrapper around a :cpp:type:`gamer::Vector`
        )delim"
    );

    vector.def(py::init<>(), "Default constructor");
    vector.def("__init__",
        [](Vector& v, double x, double y, double z){
            v = Vector({x,y,z});
            return;
        },
        "Constructor with values"
    );

    // vector.def_buffer([](Vector &v) -> py::buffer_info {
    //     return py::buffer_info(
    //         v.data(),           // Pointer to buffer
    //         sizeof(double),     // Size of one element
    //         py::format_descriptor<double>::value,
    //         1,                  // Number of dimensions
    //         { 3 },              // Buffer dimensions
    //         { sizeof(double) }  // Strides (in bytes) for each index
    //         );
    //     }
    // );

    vector.def("__neg__",
        [](Vector& v) { return -1*v; }
    );

    vector.def("__getitem__",
        [](const Vector &v, std::size_t i) -> const double& {
            if (i >= 3) throw py::index_error("Vector only contains three coordinates");
            return v[i];
        },
        "Get coordinate"
    );

    vector.def("__setitem__",
        [](Vector &v, size_t i, double val) {
            if (i >= 3) throw py::index_error("Vector only contains three coordinates");
            v[i] = val;
        },
        "Set coordinate."
    );

    vector.def("__repr__",
        [](const Vector& v){
            std::ostringstream out;
            out << "Vector("
                << v[0] << ","
                << v[1] << ","
                << v[2] << ")";
            return out.str();
        }
    );
}

} // end namespace gamer
