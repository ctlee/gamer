
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SurfaceMesh.h"

namespace py = pybind11;


PYBIND11_MODULE(pygamer, m) {
    py::class_<SurfaceMesh>(m, "SurfaceMesh")
        .def(py::init<>());
        // .def("setName", &Pet::setName)
        // .def("getName", &Pet::getName);
}