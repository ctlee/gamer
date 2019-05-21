
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SurfaceMesh.h"

namespace py = pybind11;

void init_face(py::module& pygamer){
    py::class_<Face> face(pygamer, "Face",
        R"delim(
            Wrapper around a :cpp:class:`Face`.
        )delim"
    );
    face.def(py::init<>(), "Default constructor");
    face.def(py::init<int, bool>(), "Construct with marker and selection");
    face.def(py::init<int, int, bool>(), "Construct with orientation, marker, and selection");
    face.def_readwrite("orientation", &Face::orientation, "The orientation of the face");
    face.def_readwrite("marker", &Face::marker, "Boundary marker value");
    face.def_readwrite("selected", &Face::selected, "Selection status of face");
    face.def("__repr__", &Face::to_string, "Pretty print");
}