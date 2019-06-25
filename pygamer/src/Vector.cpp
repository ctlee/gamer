
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
            Wrapper around a Vector`
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