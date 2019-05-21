
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SurfaceMesh.h"

namespace py = pybind11;

void init_vertex(py::module &);
void init_face(py::module &);
void init_surfacemesh(py::module &);


PYBIND11_MODULE(pygamer, pygamer) {
    pygamer.doc() = R"delim(
        Python wrapper around the GAMer C++ library.
    )delim";

    init_vertex(pygamer);
    init_face(pygamer);
    init_surfacemesh(pygamer);
}