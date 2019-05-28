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

#include "SurfaceMesh.h"
#include "Vertex.h"

namespace py = pybind11;

void init_SurfaceMesh(py::module& mod){
    // Bindings for SurfaceMesh
    py::class_<SurfaceMesh> SurfMeshCls(mod, "SurfaceMesh",
        R"delim(
            Python wrapper around a :cpp:type:`SurfaceMesh`.
        )delim"
    );
    SurfMeshCls.def(py::init<>(), "Default constructor.");


    /************************************
     *  INSERT/REMOVE
     ************************************/
    SurfMeshCls.def("addVertex",
        py::overload_cast<const Vertex&>(&SurfaceMesh::add_vertex),
        py::arg("data"),
        R"delim(
            Add a vertex to the mesh without specifying the key.
        )delim"
    );


    SurfMeshCls.def("insertVertex",
        py::overload_cast<const std::array<int, 1>&, const Vertex&>(&SurfaceMesh::insert<1>),
        py::arg("key"), py::arg("data") = Vertex(0,0,0,0,false),
        R"delim(
            Inserts a vertex and data based on key.

            Args:
                key (list): Array [1] of vertex key.
                data (:py:class:`Vertex`, optional): Vertex data.
        )delim"
    );


    SurfMeshCls.def("insertEdge",
        py::overload_cast<const std::array<int, 2>&, const Edge&>(&SurfaceMesh::insert<2>),
        py::arg("key"), py::arg("data") = Edge(false),
        R"delim(
            Insert an edge into the mesh.

            Args:
                key (list): Array [2] of edge key.
                data (Edge): Edge data.
        )delim"
    );


    SurfMeshCls.def("insertFace",
        py::overload_cast<const std::array<int, 3>&, const Face&>(&SurfaceMesh::insert<3>),
        py::arg("key"), py::arg("data") = Face(0,0,false),
        R"delim(
            Insert a face into the mesh.

            Args:
                key (list): Array [3] of edge key.
                data (Face): Face data.
        )delim"
    );


    SurfMeshCls.def("removeVertex",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::remove<1>),
        py::arg("key"),
        R"delim(
            Remove a vertex from the mesh.

            Args:
                key (list): Array [1] of vertex key.

            Returns:
                removed (int): Number of simplices removed
        )delim"
    );


    SurfMeshCls.def("removeVertex",
        py::overload_cast<SurfaceMesh::SimplexID<1>>(&SurfaceMesh::remove<1>),
        py::arg("key"),
        R"delim(
            Remove a vertex from the mesh.

            Args:
                key (:py:class:`VertexID`): SimplexID of vertex.

            Returns:
                removed (int): Number of simplices removed
        )delim"
    );


    SurfMeshCls.def("removeEdge",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::remove<2>),
        py::arg("key"),
        R"delim(
            Remove an edge from the mesh.

            Args:
                key (list): Array [2] of edge key.

            Returns:
                removed (int): Number of simplices removed
        )delim"
    );


    SurfMeshCls.def("removeEdge",
        py::overload_cast<SurfaceMesh::SimplexID<2>>(&SurfaceMesh::remove<2>),
        py::arg("key"),
        R"delim(
            Remove an edge from the mesh.

            Args:
                key (:py:class:`EdgeID`): SimplexID of edge.

            Returns:
                removed (int): Number of simplices removed
        )delim"
    );


    SurfMeshCls.def("removeFace",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::remove<3>),
        py::arg("key"),
        R"delim(
            Remove a face from the mesh.

            Args:
                key (list): Array [3] of face key.

            Returns:
                removed (int): Number of simplices removed
        )delim"
    );


    SurfMeshCls.def("removeFace",
        py::overload_cast<SurfaceMesh::SimplexID<3>>(&SurfaceMesh::remove<3>),
        py::arg("key"),
        R"delim(
            Remove a face from the mesh

            Args:
                key (:py:class:`FaceID`): SimplexID of face.

            Returns:
                removed (int): Number of simplices removed
        )delim"
    );


    /************************************
     * GET SIMPLEXID
     ************************************/
    // NOTE: these don't need py::return_value_policy::reference_internal
    // because they return pointers
    SurfMeshCls.def("getVertex",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::get_simplex_up<1>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (list): Array [1] of vertex key.

            Returns:
                SimplexID (:py:class:`VertexID`): The object.
        )delim"
    );

    SurfMeshCls.def("getEdge",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::get_simplex_up<2>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (list): Array [2] of edge key.

            Returns:
                SimplexID (:py:class:`EdgeID`): The object.
        )delim"
    );

    SurfMeshCls.def("getFace",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::get_simplex_up<3>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (list): Array [3] of vertex key.

            Returns:
                SimplexID (:py:class:`FaceID`): The object.
        )delim"
    );


    SurfMeshCls.def("getRoot",
        [](SurfaceMesh &mesh) -> Global&{
            return mesh.get_simplex_up().data();
        },
        py::return_value_policy::reference_internal,
        R"delim(
            Get the global metadata.

            Returns:
                data (:py:class:`Global`): Global mesh metadata.
        )delim"
    );


    /************************************
     *  GET SIMPLEX NAMES
     ************************************/
    SurfMeshCls.def("getName",
        py::overload_cast<SurfaceMesh::SimplexID<1>>(&SurfaceMesh::get_name<1>, py::const_),
        R"delim(
            Get the name of the vertex

            Args:
                SimplexID  (:py:class`VertexID`): VertexID to get the name of.

            Returns:
                key (list): Name of the vertex.
        )delim"
    );


    SurfMeshCls.def("getName",
        py::overload_cast<SurfaceMesh::SimplexID<2>>(&SurfaceMesh::get_name<2>, py::const_),
        R"delim(
            Get the name of the vertex

            Args:
                SimplexID  (:py:class`VertexID`): VertexID to get the name of.

            Returns:
                key (list): Name of the vertex.
        )delim"
    );


    SurfMeshCls.def("getName",
        py::overload_cast<SurfaceMesh::SimplexID<3>>(&SurfaceMesh::get_name<3>, py::const_),
        R"delim(
            Get the name of the vertex

            Args:
                SimplexID  (:py:class`FaceID`): FaceID to get the name of.

            Returns:
                key (list): List of vertex indices which make up the face.
        )delim"
    );


    /************************************
     *  UTILITY
     ************************************/
    SurfMeshCls.def("__repr__", &print, "Pretty print the mesh.");

    SurfMeshCls.def_property_readonly("nVertices",
        py::overload_cast<>(&SurfaceMesh::size<1>, py::const_),
        R"delim(
            Get the number of vertices.
        )delim"
    );


    SurfMeshCls.def_property_readonly("nEdges",
        py::overload_cast<>(&SurfaceMesh::size<2>, py::const_),
        R"delim(
            Get the number of edges.
        )delim"
    );


    SurfMeshCls.def_property_readonly("nFaces",
        py::overload_cast<>(&SurfaceMesh::size<3>, py::const_),
        R"delim(
            Get the number of faces.
        )delim"
    );


    /************************************
     *  ITERATORS
     ************************************/
    SurfMeshCls.def_property_readonly("vertexIDs",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<1>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            A convenience iterator over :py:class:`VertexID`.
        )delim"
    );


    SurfMeshCls.def("getVertexIDIterator",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<1>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`VertexID`.
        )delim"
    );


    SurfMeshCls.def_property_readonly("edgeIDs",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<2>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            A convenience iterator over :py:class:`EdgeID`.
        )delim"
    );


    SurfMeshCls.def("getEdgeIDIterator",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<2>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`EdgeID`.
        )delim"
    );


    SurfMeshCls.def_property_readonly("faceIDs",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<3>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            An convenience iterator over :py:class:`FaceID`.
        )delim"
    );


    SurfMeshCls.def("getFaceIDIterator",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<3>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`FaceID`.
        )delim"
    );


    /************************************
     *  FUNCTIONS
     ************************************/
    SurfMeshCls.def("getMeanCurvature", &getMeanCurvature,
        py::arg("VertexID"),
        R"delim(
            Compute the mean curvature around a vertex.

            Args:
                VertexID (:py:class:`VertexID`): SimplexID of the vertex

            Returns:
                float: Local mean curvature
        )delim"
        );


    SurfMeshCls.def("getGaussianCurvature", &getGaussianCurvature,
        py::arg("VertexID"),
        R"delim(
            Compute the Gaussian curvature around a vertex.

            Args:
                VertexID (:py:class:`VertexID`): SimplexID of the vertex

            Returns:
                float: Local Gaussian curvature
        )delim"
        );


    SurfMeshCls.def("smoothMesh", &smoothMesh,
        py::arg("maxIter"), py::arg("preserveRidges"), py::arg("verbose"),
        R"delim(
            Perform mesh smoothing.

            This operation performs an iterative weightedVertexSmooth
            and edge flipping approach.

            Args:
                maxIter (int): Maximum number of smoothing iterations
                preserveRidges (bool):  Prevent flipping of edges along ridges.
                verbose (bool): Print details to std::out
        )delim"
        );


    SurfMeshCls.def("coarse", &coarse,
        py::arg("coarseRate"),
        py::arg("flatRate"),
        py::arg("denseWeight"),
        R"delim(
            Coarsen a surface mesh.

            This operation selects vertices to decimate then removes them and
            retriangulates the resulting hole.

            Args:
                coarseRate (float): Threshold value for coarsening
                flatRate (float): Priority of decimating flat regions
                denseWeight (float): Priority of decimating dense regions
        )delim"
        );


    SurfMeshCls.def("coarse_flat",
        [](SurfaceMesh& mesh, double rate, int niter){
            for(int i = 0; i < niter; ++i) coarseIT(mesh, rate, 0.5, 0);
        },
        py::arg("rate"), py::arg("numiter"),
        R"delim(
            Coarsen flat regions of a surface mesh.

            Args:
                rate (float): Threshold value
                numiter (int): Number of iterations to run
        )delim"
    );


    SurfMeshCls.def("coarse_dense",
        [](SurfaceMesh& mesh, double rate, int niter){
            for(int i = 0; i < niter; ++i) coarseIT(mesh, rate, 0, 10);
        },
        py::arg("rate"), py::arg("numiter"),
        R"delim(
            Coarsen dense regions of a surface mesh.

            Args:
                rate (float): Threshold value
                numiter (int): Number of iterations to run
        )delim"
    );


    SurfMeshCls.def("normalSmooth", &normalSmooth,
        R"delim(
            Perform smoothing of mesh face normals.
        )delim"
    );


    SurfMeshCls.def("fillHoles", &fillHoles,
        R"delim(
            Fill holes in the mesh.
        )delim"
    );


    /************************************
     *  ORIENTATION
     ************************************/
    SurfMeshCls.def("init_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::init_orientation<SurfaceMesh>),
        R"delim(
            Initialize mesh orientations.
        )delim"
    );


    SurfMeshCls.def("check_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::check_orientation<SurfaceMesh>),
        R"delim(
            Check consistency and assign Face orientations.

            :py:func:`SurfaceMesh.init_orientation` must be called before this
            function can operate.
        )delim"
    );


    SurfMeshCls.def("clear_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::clear_orientation<SurfaceMesh>),
        R"delim(
            Clear orientation of all faces.
        )delim"
    );


    SurfMeshCls.def("compute_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::compute_orientation<SurfaceMesh>),
        R"delim(
            Compute orientations of the mesh.
        )delim"
    );


    SurfMeshCls.def("flipNormals", &flipNormals,
        R"delim(
            Flip all of the mesh normals.
        )delim"
    );


    SurfMeshCls.def("correctNormals",
        [](SurfaceMesh& mesh){
            if(getVolume(mesh) < 0){
                for(auto &face : mesh.get_level<3>())
                    face.orientation *= -1;
            }
        },
        R"delim(
            Set normals to be outward facing.
        )delim"
    );
}