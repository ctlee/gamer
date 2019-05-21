
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
    vertex.def_readwrite("position", (std::array<double, 3> (Vertex::*)) &Vertex::position, "Vertex position");
    vertex.def_readwrite("marker", &Vertex::marker, "Boundary marker value");
    vertex.def_readwrite("selected", &Vertex::selected, "Selection status of vertex");
    vertex.def("__repr__", &Vertex::to_string);
}