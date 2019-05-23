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

namespace py = pybind11;

using VertexID_Iterator = SurfaceMesh::SimplexIDIterator<1>;
using VertexData_Iterator = SurfaceMesh::DataIterator<1>;

using FaceID_Iterator = SurfaceMesh::SimplexIDIterator<3>;
using FaceData_Iterator = SurfaceMesh::DataIterator<3>;

void init_SurfaceMesh(py::module& mod){
    // Bindings for SurfaceMesh
    py::class_<SurfaceMesh> SurfMesh(mod, "SurfaceMesh",
        R"delim(
            Python wrapper around a :cpp:type:`SurfaceMesh`.
        )delim"
    );
    SurfMesh.def(py::init<>(), "Default constructor.");

    SurfMesh.def("addVertex",
        py::overload_cast<const Vertex&>(&SurfaceMesh::add_vertex),
        R"delim(
            Add a vertex to the mesh without specifying the key.
        )delim"
    );

    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::insert<1>),
        R"delim(
            Inserts a vertex based on key.

            Args:
                key (array): Array [1] of vertex key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 1>&, const Vertex&>(&SurfaceMesh::insert<1>),
        R"delim(
            Inserts a vertex and data based on key.

            Args:
                key (array): Array [1] of vertex key.
                data (:py:class:`Vertex`): Vertex data.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::insert<2>),
        R"delim(
            Insert an edge into the mesh

            Args:
                key (array): Array [2] of edge key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 2>&, const Edge&>(&SurfaceMesh::insert<2>),
        R"delim(
            Insert an edge into the mesh

            Args:
                key (array): Array [2] of edge key.
                data (Edge): Edge data.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::insert<3>),
        R"delim(
            Insert a face into the mesh

            Args:
                key (array): Array [3] of face key.
        )delim"
    );
    SurfMesh.def("insert",
        py::overload_cast<const std::array<int, 3>&, const Face&>(&SurfaceMesh::insert<3>),
        R"delim(
            Insert a face into the mesh

            Args:
                arg1 (array): Array [3] of edge key.
                arg2 (Face): Face data.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::get_simplex_up<1>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [1] of vertex key.

            Returns:
                SimplexID (:py:class:`VertexID`): The object.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::get_simplex_up<2>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [2] of edge key.

            Returns:
                SimplexID (:py:class:`EdgeID`): The object.
        )delim"
    );

    SurfMesh.def("get_simplex_up",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::get_simplex_up<3>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (array): Array [3] of vertex key.

            Returns:
                SimplexID (:py:class:`FaceID`): The object.
        )delim"
    );

    // SurfMesh.def("print", &print, "Print a surface mesh");
    SurfMesh.def("__repr__", &print);


    /************************************
     *  ITERATORS
     ************************************/
    SurfMesh.def_property_readonly("vertexIDs",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<1>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            A convenience iterator over :py:class:`VertexID`.
        )delim"
    );

    // SurfMesh.def_property_readonly("vertices",
    //     [](const SurfaceMesh& m){
    //         auto it = m.get_level<1>();
    //         return py::make_iterator(it.begin(), it.end());
    //     },
    //     py::keep_alive<0, 1>(),
    //     R"delim(
    //         A convenience iterator over :py:class:`Vertex`s.
    //     )delim"
    // );

    SurfMesh.def("getVertexIDIterator",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<1>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`VertexID`.
        )delim"
    );

    // SurfMesh.def("getVertexDataIterator",
    //     [](const SurfaceMesh& m){
    //         auto it = m.get_level<1>();
    //         return py::make_iterator(it.begin(), it.end());
    //     },
    //     py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
    //     R"delim(
    //         Get an iterator over :py:class:`Vertex`s.
    //     )delim"
    // );

    SurfMesh.def_property_readonly("edgeIDs",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<2>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            A convenience iterator over :py:class:`EdgeID`.
        )delim"
    );

    SurfMesh.def("getEdgeIDIterator",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<2>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`EdgeID`.
        )delim"
    );

    SurfMesh.def_property_readonly("faceIDs",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<3>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(),
        R"delim(
            An convenience iterator over :py:class:`FaceID`.
        )delim"
    );

    // SurfMesh.def_property_readonly("faces",
    //     [](const SurfaceMesh& m){
    //         auto it = m.get_level<3>();
    //         return py::make_iterator(it.begin(), it.end());
    //     },
    //     py::keep_alive<0, 1>(),
    //     R"delim(
    //         A convenience iterator over :py:class:`Vertex`s.
    //     )delim"
    // );

    SurfMesh.def("getFaceIDIterator",
        [](const SurfaceMesh& m){
            auto it = m.get_level_id<3>();
            return py::make_iterator(it.begin(), it.end());
        },
        py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
        R"delim(
            Get an iterator over :py:class:`FaceID`.
        )delim"
    );

    // SurfMesh.def("getFaceDataIterator",
    //     [](const SurfaceMesh& m){
    //         auto it = m.get_level<3>();
    //         return py::make_iterator(it.begin(), it.end());
    //     },
    //     py::keep_alive<0, 1>(), /* Essential: keep object alive while iterator exists */
    //     R"delim(
    //         Get an iterator over :py:class:`Face`s.
    //     )delim"
    // );


    /************************************
     *  FUNCTIONS
     ************************************/
    SurfMesh.def("getMeanCurvature", &getMeanCurvature,
        py::arg("VertexID"),
        R"delim(
            Compute the mean curvature around a vertex.

            Args:
                VertexID (:py:class:`VertexID`): SimplexID of the vertex

            Returns:
                float: Local mean curvature
        )delim"
        );


    SurfMesh.def("getGaussianCurvature", &getGaussianCurvature,
        py::arg("VertexID"),
        R"delim(
            Compute the Gaussian curvature around a vertex.

            Args:
                VertexID (:py:class:`VertexID`): SimplexID of the vertex

            Returns:
                float: Local Gaussian curvature
        )delim"
        );


    SurfMesh.def("smoothMesh", &smoothMesh,
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


    SurfMesh.def("coarse", &coarse,
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


    SurfMesh.def("coarse_flat",
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


    SurfMesh.def("coarse_dense",
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


    SurfMesh.def("normalSmooth", &normalSmooth,
        R"delim(
            Perform smoothing of mesh face normals.
        )delim"
    );


    SurfMesh.def("fillHoles", &fillHoles,
        R"delim(
            Fill holes in the mesh
        )delim"
    );


    /************************************
     *  ORIENTATION
     ************************************/
    SurfMesh.def("init_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::init_orientation<SurfaceMesh>),
        R"delim(
            Initialize mesh orientations.
        )delim"
    );


    SurfMesh.def("check_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::check_orientation<SurfaceMesh>),
        R"delim(
            Check consistency and assign Face orientations.

            :py:func:`SurfaceMesh.init_orientation` must be called before this
            function can operate.
        )delim"
    );


    SurfMesh.def("clear_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::clear_orientation<SurfaceMesh>),
        R"delim(
            Clear orientation of all faces.
        )delim"
    );


    SurfMesh.def("compute_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::compute_orientation<SurfaceMesh>),
        R"delim(
            Compute orientations of the mesh.
        )delim"
    );


    SurfMesh.def("flipNormals", &flipNormals,
        R"delim(
            Flip all of the mesh normals.
        )delim"
    );


    SurfMesh.def("correctNormals",
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