/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2019
 * by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
 *    and Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * ***************************************************************************
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gamer/TetMesh.h"
#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

void init_TetMesh(py::module& mod){
    // Bindings for TetMesh
    py::class_<TetMesh> TetMeshCls(mod, "TetMesh",
        R"delim(
            Python wrapper around a :cpp:type:`gamer::TetMesh`.
        )delim"
    );
    TetMeshCls.def(py::init<>(), "Default constructor.");

    /************************************
     *  INSERT/REMOVE
     ************************************/
    TetMeshCls.def("addVertex",
        py::overload_cast<const TMVertex&>(&TetMesh::add_vertex),
        py::arg("data"),
        R"delim(
            Add a vertex to the mesh without specifying the key.
        )delim"
    );


    TetMeshCls.def("insertVertex",
        py::overload_cast<const std::array<int, 1>&, const TMVertex&>(&TetMesh::insert<1>),
        py::arg("key"), py::arg("data") = TMVertex(),
        R"delim(
            Inserts a vertex and data based on key.

            Args:
                key (:py:class:`list`): Array [1] of vertex key.
                data (:py:class:`Vertex`): Vertex data.
        )delim"
    );


    TetMeshCls.def("insertEdge",
        py::overload_cast<const std::array<int, 2>&, const TMEdge&>(&TetMesh::insert<2>),
        py::arg("key"), py::arg("data") = TMEdge(),
        R"delim(
            Insert an edge into the mesh.

            Args:
                key (:py:class:`list`): Array [2] of edge key.
                data (tetmesh.Edge): Edge data.
        )delim"
    );


    TetMeshCls.def("insertFace",
        py::overload_cast<const std::array<int, 3>&, const TMFace&>(&TetMesh::insert<3>),
        py::arg("key"), py::arg("data") = TMFace(),
        R"delim(
            Insert a face into the mesh.

            Args:
                key (:py:class:`list`): Array [3] of face key.
                data (tetmesh.Face): Face data.
        )delim"
    );


    TetMeshCls.def("insertCell",
        py::overload_cast<const std::array<int, 4>&, const TMCell&>(&TetMesh::insert<4>),
        py::arg("key"), py::arg("data") = TMCell(),
        R"delim(
            Insert a face into the mesh.

            Args:
                key (:py:class:`list`): Array [4] of cell key.
                data (Cell): Cell data.
        )delim"
    );


    TetMeshCls.def("removeVertex",
        py::overload_cast<const std::array<int, 1>&>(&TetMesh::remove<1>),
        py::arg("key"),
        R"delim(
            Remove a vertex from the mesh.

            Args:
                key (:py:class:`list`): Array [1] of vertex key.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeVertex",
        py::overload_cast<TetMesh::SimplexID<1>>(&TetMesh::remove<1>),
        py::arg("key"),
        R"delim(
            Remove a vertex from the mesh.

            Args:
                key (:py:class:`VertexID`): SimplexID of vertex.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeEdge",
        py::overload_cast<const std::array<int, 2>&>(&TetMesh::remove<2>),
        py::arg("key"),
        R"delim(
            Remove an edge from the mesh.

            Args:
                key (:py:class:`list`): Array [2] of edge key.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeEdge",
        py::overload_cast<TetMesh::SimplexID<2>>(&TetMesh::remove<2>),
        py::arg("key"),
        R"delim(
            Remove an edge from the mesh.

            Args:
                key (:py:class:`EdgeID`): SimplexID of edge.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeFace",
        py::overload_cast<const std::array<int, 3>&>(&TetMesh::remove<3>),
        py::arg("key"),
        R"delim(
            Remove a face from the mesh.

            Args:
                key (:py:class:`list`): Array [3] of face key.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeFace",
        py::overload_cast<TetMesh::SimplexID<3>>(&TetMesh::remove<3>),
        py::arg("key"),
        R"delim(
            Remove a face from the mesh.

            Args:
                key (:py:class:`FaceID`): SimplexID of face.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeCell",
        py::overload_cast<const std::array<int, 4>&>(&TetMesh::remove<4>),
        py::arg("key"),
        R"delim(
            Remove a cell from the mesh.

            Args:
                key (:py:class:`list`): Array [3] of cell key.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    TetMeshCls.def("removeCell",
        py::overload_cast<TetMesh::SimplexID<4>>(&TetMesh::remove<4>),
        py::arg("key"),
        R"delim(
            Remove a cell from the mesh.

            Args:
                key (:py:class:`CellID`): SimplexID of cell.

            Returns:
                removed (int): Number of simplices removed.
        )delim"
    );


    /************************************
     * GET SIMPLEXID
     ************************************/
    // NOTE: these don't need py::return_value_policy::reference_internal
    // because they return pointers
    TetMeshCls.def("getVertex",
        py::overload_cast<const std::array<int, 1>&>(&TetMesh::get_simplex_up<1>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (:py:class:`list`): Array [1] of vertex key.

            Returns:
                :py:class:`VertexID`: The object.
        )delim"
    );

    TetMeshCls.def("getEdge",
        py::overload_cast<const std::array<int, 2>&>(&TetMesh::get_simplex_up<2>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (:py:class:`list`): Array [2] of edge key.

            Returns:
                :py:class:`EdgeID`: The object.
        )delim"
    );

    TetMeshCls.def("getFace",
        py::overload_cast<const std::array<int, 3>&>(&TetMesh::get_simplex_up<3>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (:py:class:`list`): Array [3] of vertex key.

            Returns:
                :py:class:`FaceID`: The object.
        )delim"
    );


    TetMeshCls.def("getCell",
        py::overload_cast<const std::array<int, 4>&>(&TetMesh::get_simplex_up<4>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (:py:class:`list`): Array [4] of vertex key.

            Returns:
                :py:class:`CellID`: The object.
        )delim"
    );


    TetMeshCls.def("getRoot",
        [](TetMesh &mesh) -> TMGlobal&{
            return mesh.get_simplex_up().data();
        },
        py::return_value_policy::reference_internal,
        R"delim(
            Get the global metadata.

            Returns:
                :py:class:`Global`: Global mesh metadata.
        )delim"
    );


    /************************************
     *  GET SIMPLEX NAMES
     ************************************/
    TetMeshCls.def("getName",
        py::overload_cast<TetMesh::SimplexID<1>>(&TetMesh::get_name<1>, py::const_),
        R"delim(
            Get the name of the vertex.

            Args:
                SimplexID  (:py:class:`VertexID`): VertexID to get the name of.

            Returns:
                :py:class:`list`: Name of the vertex.
        )delim"
    );


    TetMeshCls.def("getName",
        py::overload_cast<TetMesh::SimplexID<2>>(&TetMesh::get_name<2>, py::const_),
        R"delim(
            Get the name of the edge.

            Args:
                SimplexID  (:py:class:`EdgeID`): EdgeID to get the name of.

            Returns:
                :py:class:`list`: List of vertex indices which make up the edge.
        )delim"
    );


    TetMeshCls.def("getName",
        py::overload_cast<TetMesh::SimplexID<3>>(&TetMesh::get_name<3>, py::const_),
        R"delim(
            Get the name of the face.

            Args:
                SimplexID  (:py:class:`FaceID`): FaceID to get the name of.

            Returns:
                :py:class:`list`: List of vertex indices which make up the face.
        )delim"
    );


    TetMeshCls.def("getName",
        py::overload_cast<TetMesh::SimplexID<4>>(&TetMesh::get_name<4>, py::const_),
        R"delim(
            Get the name of the cell.

            Args:
                SimplexID  (:py:class:`CellID`): CellID to get the name of.

            Returns:
                :py:class:`list`: List of vertex indices which make up the cell.
        )delim"
    );

    /************************************
     *  UTILITY
     ************************************/
    TetMeshCls.def_property_readonly("nVertices",
        py::overload_cast<>(&TetMesh::size<1>, py::const_),
        R"delim(
            Get the number of vertices.

            Returns:
                :py:class:`int`: Number of vertices.
        )delim"
    );


    TetMeshCls.def_property_readonly("nEdges",
        py::overload_cast<>(&TetMesh::size<2>, py::const_),
        R"delim(
            Get the number of edges.

            Return:
                :py:class:`int`: Number of edges.
        )delim"
    );


    TetMeshCls.def_property_readonly("nFaces",
        py::overload_cast<>(&TetMesh::size<3>, py::const_),
        R"delim(
            Get the number of faces.

            Return:
                :py:class:`int`: Number of faces.
        )delim"
    );


    TetMeshCls.def_property_readonly("nCells",
        py::overload_cast<>(&TetMesh::size<4>, py::const_),
        R"delim(
            Get the number of cells.

            Return:
                :py:class:`int`: Number of cells.
        )delim"
    );

    TetMeshCls.def("onBoundary",
        py::overload_cast<const TetMesh::SimplexID<1>>(&TetMesh::onBoundary<1>, py::const_),
        R"delim(
            Check if a vertex is on a boundary.

            Args:
                SimplexID (:py:class:`VertexID`): VertexID to check.

            Returns:
                :py:class:`bool`: True if vertex is a member of a face on a boundary.
        )delim"
    );


    TetMeshCls.def("onBoundary",
        py::overload_cast<const TetMesh::SimplexID<2>>(&TetMesh::onBoundary<2>, py::const_),
        R"delim(
            Check if an edge is on a boundary.

            Args:
                SimplexID (:py:class:`EdgeID`): EdgeID to check.

            Returns:
                :py:class:`bool`: True if edge is a member of a face on a boundary.
        )delim"
    );

    TetMeshCls.def("onBoundary",
        py::overload_cast<const TetMesh::SimplexID<3>>(&TetMesh::onBoundary<3>, py::const_),
        R"delim(
            Check if a face is on a boundary.

            Args:
                SimplexID (:py:class:`FaceID`): FaceID to check.

            Returns:
                :py:class:`bool`: True if face is on a boundary.
        )delim"
    );

    TetMeshCls.def("onBoundary",
        py::overload_cast<const TetMesh::SimplexID<4>>(&TetMesh::onBoundary<4>, py::const_),
        R"delim(
            Check if a cell is on a boundary.

            Args:
                SimplexID (:py:class:`CellID`): CellID to check.

            Returns:
                :py:class:`bool`: True if cell has a face on a boundary.
        )delim"
    );

    TetMeshCls.def("nearBoundary",
        py::overload_cast<const TetMesh::SimplexID<1>>(&TetMesh::nearBoundary<1>, py::const_),
        R"delim(
            Check if a vertex is near a boundary.

            Args:
                SimplexID (:py:class:`VertexID`): VertexID to check.

            Returns:
                :py:class:`bool`: True if vertex touches a boundary.
        )delim"
    );

    TetMeshCls.def("nearBoundary",
        py::overload_cast<const TetMesh::SimplexID<2>>(&TetMesh::nearBoundary<2>, py::const_),
        R"delim(
            Check if an edge is on a boundary.

            Args:
                SimplexID (:py:class:`EdgeID`): EdgeID to check.

            Returns:
                :py:class:`bool`: True if edge touches a boundary.
        )delim"
    );

    TetMeshCls.def("nearBoundary",
        py::overload_cast<const TetMesh::SimplexID<3>>(&TetMesh::nearBoundary<3>, py::const_),
        R"delim(
            Check if a face is on a boundary.

            Args:
                SimplexID (:py:class:`FaceID`): FaceID to check.

            Returns:
                :py:class:`bool`: True if face touches a boundary.
        )delim"
    );

    TetMeshCls.def("nearBoundary",
        py::overload_cast<const TetMesh::SimplexID<4>>(&TetMesh::nearBoundary<4>, py::const_),
        R"delim(
            Check if a cell is on a boundary.

            Args:
                SimplexID (:py:class:`CellID`): CellID to check.

            Returns:
                :py:class:`bool`: True if cell touches a boundary.
        )delim"
    );

    TetMeshCls.def("extractSurface",
        &extractSurface,
        R"delim(
            Extract the surface of the TetMesh.

            Args:
                tetmesh (TetMesh): Tetrahedral mesh to extract from.

            Returns:
                :py:class:`SurfaceMesh`: Surface mesh of surface.
        )delim"
    );

    /************************************
     *  ITERATORS
     ************************************/
    TetMeshCls.def_property_readonly("vertexIDs",
        [](const TetMesh& m){
            auto it = m.get_level_id<1>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            A convenience iterator over :py:class:`VertexID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`VertexID`.
        )delim"
    );


    TetMeshCls.def("getVertexIDIterator",
        [](const TetMesh& m){
            auto it = m.get_level_id<1>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`VertexID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`VertexID`.
        )delim"
    );


    TetMeshCls.def_property_readonly("edgeIDs",
        [](const TetMesh& m){
            auto it = m.get_level_id<2>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            A convenience iterator over :py:class:`EdgeID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`EdgeID`.
        )delim"
    );


    TetMeshCls.def("getEdgeIDIterator",
        [](const TetMesh& m){
            auto it = m.get_level_id<2>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`EdgeID`.


            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`EdgeID`.
        )delim"
    );


    TetMeshCls.def_property_readonly("faceIDs",
        [](const TetMesh& m){
            auto it = m.get_level_id<3>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            An convenience iterator over :py:class:`FaceID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`FaceID`.
        )delim"
    );


    TetMeshCls.def("getFaceIDIterator",
        [](const TetMesh& m){
            auto it = m.get_level_id<3>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`FaceID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`FaceID`.
        )delim"
    );


    TetMeshCls.def_property_readonly("cellIDs",
        [](const TetMesh& m){
            auto it = m.get_level_id<4>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            An convenience iterator over :py:class:`CellID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`CellID`.
        )delim"
    );


    TetMeshCls.def("getCellIDIterator",
        [](const TetMesh& m){
            auto it = m.get_level_id<4>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`CellID`.

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`CellID`.
        )delim"
    );


    /************************************
     *  ORIENTATION
     ************************************/
    TetMeshCls.def("init_orientation",
        py::overload_cast<TetMesh&>(&casc::init_orientation<TetMesh>),
        R"delim(
            Initialize mesh orientations.
        )delim"
    );


    TetMeshCls.def("check_orientation",
        py::overload_cast<TetMesh&>(&casc::check_orientation<TetMesh>),
        R"delim(
            Check consistency and assign Face orientations.

            :py:func:`TetMesh.init_orientation` must be called before this
            function can operate.

            Returns:
                :py:class:`tuple` (:py:class:`int`, :py:class:`bool`, :py:class:`bool`): Tuple containing: number of connected components, orientability status, and manifoldness.
        )delim"
    );


    TetMeshCls.def("clear_orientation",
        py::overload_cast<TetMesh&>(&casc::clear_orientation<TetMesh>),
        R"delim(
            Clear orientation of all faces.
        )delim"
    );


    TetMeshCls.def("compute_orientation",
        py::overload_cast<TetMesh&>(&casc::compute_orientation<TetMesh>),
        R"delim(
            Compute orientations of the mesh.

            This is equivalent to the following series of calls:

            >>> init_orientation(mesh)
            >>> clear_orientation(mesh)
            >>> check_orientation(mesh)


            Returns:
                :py:class:`tuple` (:py:class:`int`, :py:class:`bool`, :py:class:`bool`): Tuple containing: number of connected components, orientability status, and manifoldness.
        )delim"
    );
}

} // end namespace gamer