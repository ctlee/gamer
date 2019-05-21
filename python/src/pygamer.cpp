
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

    py::module SurfMeshMod = pygamer.def_submodule("surfmesh", "Methods associated with SurfaceMesh");

    SurfMeshMod.def("readOFF", &readOFF);
    SurfMeshMod.def("writeOFF", &writeOFF);
    SurfMeshMod.def("readOBJ", &readOBJ);
    SurfMeshMod.def("writeOBJ", &writeOBJ);
    SurfMeshMod.def("getMeanCurvature", &getMeanCurvature);
    SurfMeshMod.def("getGaussianCurvature", &getGaussianCurvature);
    SurfMeshMod.def("smoothMesh", &smoothMesh);
    SurfMeshMod.def("coarse", &coarse);
    SurfMeshMod.def("normalSmooth", &normalSmooth);
    SurfMeshMod.def("cube", &cube);
    SurfMeshMod.def("sphere", &sphere);
}