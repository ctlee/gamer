

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
    pygamer.doc() = R"delim(
        Python wrapper around the GAMer C++ library.

        Classes starting with SM are associated with SurfaceMesh.
    )delim";

    // Bindings for VertexID
    py::class_<VertexID> vid(pygamer, "SMVertexID",
        R"delim(
            SurfaceMesh::SimplexID<1> object

            This is a token to represent a 1-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    vid.def(py::init<>(), "Default constructor");
    vid.def("data", py::overload_cast<>(&VertexID::data), "Access the data.");

    // Bindings for EdgeID
    py::class_<EdgeID> eid(pygamer, "SMEdgeID",
        R"delim(
            SurfaceMesh::SimplexID<2> object

            This is a token to represent a 2-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    eid.def(py::init<>(), "Default constructor");
    eid.def("data", py::overload_cast<>(&EdgeID::data), "Access the data.");

    // Bindings for FaceID
    py::class_<FaceID> fid(pygamer, "SMFaceID",
        R"delim(
            SurfaceMesh::SimplexID<3> object

            This is a token to represent a 3-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    fid.def(py::init<>(), "Default constructor");
    fid.def("data", py::overload_cast<>(&FaceID::data), "Access the data.");

    // Bindings for SurfaceMesh
    py::class_<SurfaceMesh> SurfMesh(pygamer, "SurfaceMesh",
        R"delim(
            Python wrapper around SurfaceMesh class.
        )delim"
    );
    SurfMesh.def(py::init<>(), "Default constructor");

    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::insert<1>),
        R"delim(
            Inserts a vertex based on key.

            Args:
                key (array): Array [1] of vertex key.

            Returns:
                void
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 1>&, const Vertex&>(&SurfaceMesh::insert<1>),
        R"delim(
            Inserts a vertex and data based on key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::insert<2>)
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&, const Edge&>(&SurfaceMesh::insert<2>)
    );
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::insert<3>));
    SurfMesh.def("insert", py::overload_cast<const std::array<int, 3>&, const Face&>(&SurfaceMesh::insert<3>));

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::get_simplex_up<1>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [1] of vertex key.

            Returns:
                SimplexID (SMVertexID): The object.
        )delim"
    );

    SurfMesh.def("print", &print, "Print a surface mesh");
}