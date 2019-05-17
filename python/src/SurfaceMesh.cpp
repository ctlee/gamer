

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SurfaceMesh.h"

namespace py = pybind11;

using VertexID_Iterator = SurfaceMesh::SimplexIDIterator<1>;
using VertexData_Iterator = SurfaceMesh::DataIterator<1>;

using FaceID_Iterator = SurfaceMesh::SimplexIDIterator<3>;
using FaceData_Iterator = SurfaceMesh::DataIterator<3>;

using VertexID = SurfaceMesh::SimplexID<1>;
using EdgeID = SurfaceMesh::SimplexID<2>;
using FaceID = SurfaceMesh::SimplexID<3>;

PYBIND11_MODULE(pygamer, pygamer) {

    py::class_<VertexID> vid(pygamer, "VertexID");
    vid.def(py::init<>());
    vid.def("data", py::overload_cast<>(&VertexID::data));

    py::class_<EdgeID> eid(pygamer, "EdgeID");
    eid.def(py::init<>());
    eid.def("data", py::overload_cast<>(&EdgeID::data));

    py::class_<FaceID> fid(pygamer, "FaceID");
    fid.def(py::init<>());
    fid.def("data", py::overload_cast<>(&FaceID::data));

    py::class_<SurfaceMesh> SurfMesh(pygamer, "SurfaceMesh");
    SurfMesh.def(py::init<>());

    SurfMesh.def("insert", py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::insert<1>));
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 1>&, const Vertex&>(&SurfaceMesh::insert<1>));
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::insert<2>));
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 2>&, const Edge&>(&SurfaceMesh::insert<2>));
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::insert<3>));
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 3>&, const Face&>(&SurfaceMesh::insert<3>));

    SurfMesh.def("get_simplex_up", py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::get_simplex_up<1>, py::const_));

    SurfMesh.def("print", &print, "Print a surface mesh");
}