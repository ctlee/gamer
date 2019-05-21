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


void init_surfacemesh(py::module& pygamer){
    // Bindings for VertexID
    py::class_<VertexID> vid(pygamer, "SMVertexID",
        R"delim(
            Wrapper around :cpp:type:`SurfaceMesh`::SimplexID<1> object. :py:

            This is a token to represent a 1-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    // vid.def(py::init<>(), "Default constructor");
    vid.def("data",
            py::overload_cast<>(&VertexID::data),
            py::return_value_policy::reference_internal,
            R"delim(
                Access the data stored on the 1-simplex

                Returns:
                    Vertex (:py:class:`Vertex`): Vertex data
            )delim"
        );

    // Bindings for EdgeID
    py::class_<EdgeID> eid(pygamer, "SMEdgeID",
        R"delim(
            Wrapper around :cpp:type:`SurfaceMesh`::SimplexID<2> object

            This is a token to represent a 2-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    // eid.def(py::init<>(), "Default constructor");
    eid.def("data", py::overload_cast<>(&EdgeID::data),
        py::return_value_policy::reference_internal,
        R"delim(
            Access the data stored on the edge.
        )delim"
    );

    // Bindings for FaceID
    py::class_<FaceID> fid(pygamer, "SMFaceID",
        R"delim(
            Wrapper around :cpp:type:`SurfaceMesh`::SimplexID<3> object

            This is a token to represent a 3-simplex object. It serves as a
            reference to the actual object.
        )delim"
    );
    // fid.def(py::init<>(), "Default constructor");
    fid.def("data",
        py::overload_cast<>(&FaceID::data),
        py::return_value_policy::reference_internal, "Access the data.");

    // Bindings for SurfaceMesh
    py::class_<SurfaceMesh> SurfMesh(pygamer, "SurfaceMesh",
        R"delim(
            Python wrapper around :cpp:type:`SurfaceMesh`.
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
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::insert<2>),
        R"delim(
            Inserts an edge based on key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&, const Edge&>(&SurfaceMesh::insert<2>),
        R"delim(
            Inserts an edge and data based on key.
        )delim"
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

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::get_simplex_up<2>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [1] of vertex key.

            Returns:
                SimplexID (SMVertexID): The object.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::get_simplex_up<3>, py::const_),
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