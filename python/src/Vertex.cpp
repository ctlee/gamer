
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Vertex.h"

namespace py = pybind11;

void init_vertex(py::module& pygamer){
    py::class_<Vertex> vertex(pygamer, "Vertex",
    R"delim(
    Wrapper around a :cpp:class:`Vertex`.
    )delim"
    );
    vertex.def(py::init<>(), "Default constructor");
    vertex.def(py::init<double, double, double>(), "Other constructor");

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