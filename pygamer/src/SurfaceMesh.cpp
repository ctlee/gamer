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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

void init_SurfaceMesh(py::module& mod){
    // Bindings for SurfaceMesh
    py::class_<SurfaceMesh> SurfMeshCls(mod, "SurfaceMesh",
        R"delim(
            Python wrapper around a :cpp:type:`gamer::SurfaceMesh`.
        )delim"
    );
    SurfMeshCls.def(py::init<>(), "Default constructor.");


    /************************************
     *  INSERT/REMOVE
     ************************************/
    SurfMeshCls.def("addVertex",
        py::overload_cast<const SMVertex&>(&SurfaceMesh::add_vertex),
        py::arg("data"),
        R"delim(
            Add a vertex to the mesh without specifying the key.

            Returns:
                :py:class:`int`: Key value of the new vertex.
        )delim"
    );


    SurfMeshCls.def("addVertex",
        [](SurfaceMesh& mesh, double x, double y, double z, int marker, bool sel){
            mesh.add_vertex(SMVertex(x,y,z,marker,sel));
        },
        py::arg("x") = 0,
        py::arg("y") = 0,
        py::arg("z") = 0,
        py::arg("marker") = -1,
        py::arg("selected") = false,
        R"delim(
            Add a vertex to the mesh without specifying the key.

            Args:
                x (:py:class:`float`): X position.
                y (:py:class:`float`): Y position.
                z (:py:class:`float`): Z position.
                marker (:py:class:`int`): Marker value.
                selected (:py:class:`bool`): Selection status.

            Returns:
                :py:class:`int`: Name of the new vertex.
        )delim"
    );


    SurfMeshCls.def("insertVertex",
        py::overload_cast<const std::array<int, 1>&, const SMVertex&>(&SurfaceMesh::insert<1>),
        py::arg("key"), py::arg("data") = SMVertex(0,0,0,-1,false),
        R"delim(
            Inserts a vertex and data based on key.

            Args:
                key (:py:class:`list`): Array [1] of vertex key.
                data (:py:class:`Vertex`, optional): Vertex data.
        )delim"
    );


    SurfMeshCls.def("insertEdge",
        py::overload_cast<const std::array<int, 2>&, const SMEdge&>(&SurfaceMesh::insert<2>),
        py::arg("key"), py::arg("data") = SMEdge(false),
        R"delim(
            Insert an edge into the mesh.

            Args:
                key (:py:class:`list`): Array [2] of edge key.
                data (:py:class:`Edge`): Edge data.
        )delim"
    );


    SurfMeshCls.def("insertFace",
        py::overload_cast<const std::array<int, 3>&, const SMFace&>(&SurfaceMesh::insert<3>),
        py::arg("key"), py::arg("data") = SMFace(0,0,false),
        R"delim(
            Insert a face into the mesh.

            Args:
                key (:py:class:`list`): Array [3] of edge key.
                data (:py:class:`Face`): Face data.
        )delim"
    );


    SurfMeshCls.def("removeVertex",
        py::overload_cast<const std::array<int, 1>&>(&SurfaceMesh::remove<1>),
        py::arg("key"),
        R"delim(
            Remove a vertex from the mesh.

            Args:
                key (:py:class:`list`): Array [1] of vertex key.

            Returns:
                :py:class:`int`: Number of simplices removed
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
                :py:class:`int`: Number of simplices removed.
        )delim"
    );


    SurfMeshCls.def("removeEdge",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::remove<2>),
        py::arg("key"),
        R"delim(
            Remove an edge from the mesh.

            Args:
                key (:py:class:`list`): Array [2] of edge key.

            Returns:
                :py:class:`int`: Number of simplices removed.
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
                :py:class:`int`: Number of simplices removed.
        )delim"
    );


    SurfMeshCls.def("removeFace",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::remove<3>),
        py::arg("key"),
        R"delim(
            Remove a face from the mesh.

            Args:
                key (:py:class:`list`): Array [3] of face key.

            Returns:
                :py:class:`int`: Number of simplices removed.
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
                :py:class:`int`: Number of simplices removed.
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
                key (:py:class:`list`): Array [1] of vertex key.

            Returns:
                :py:class:`VertexID`: Simplex of interest.
        )delim"
    );

    SurfMeshCls.def("getEdge",
        py::overload_cast<const std::array<int, 2>&>(&SurfaceMesh::get_simplex_up<2>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (:py:class:`list`): Array [2] of edge key.

            Returns:
               :py:class:`EdgeID`: Simplex of interest.
        )delim"
    );

    SurfMeshCls.def("getFace",
        py::overload_cast<const std::array<int, 3>&>(&SurfaceMesh::get_simplex_up<3>, py::const_),
        R"delim(
            Get a simplex by key.

            Args:
                key (:py:class:`list`): Array [3] of vertex key.

            Returns:
                :py:class:`FaceID`: The object.
        )delim"
    );


    SurfMeshCls.def("getRoot",
        [](SurfaceMesh &mesh) -> SMGlobal&{
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
    SurfMeshCls.def("getName",
        py::overload_cast<SurfaceMesh::SimplexID<1>>(&SurfaceMesh::get_name<1>, py::const_),
        R"delim(
            Get the name of the vertex.

            Args:
                SimplexID (:py:class:`VertexID`): VertexID to get the name of.

            Returns:
                :py:class:`list`: Name of the vertex.
        )delim"
    );


    SurfMeshCls.def("getName",
        py::overload_cast<SurfaceMesh::SimplexID<2>>(&SurfaceMesh::get_name<2>, py::const_),
        R"delim(
            Get the name of the edge.

            Args:
                SimplexID (:py:class:`EdgeID`): VertexID to get the name of.

            Returns:
                :py:class:`list`: List of vertex keys which make up the edge.
        )delim"
    );


    SurfMeshCls.def("getName",
        py::overload_cast<SurfaceMesh::SimplexID<3>>(&SurfaceMesh::get_name<3>, py::const_),
        R"delim(
            Get the name of the face.

            Args:
                SimplexID (:py:class:`FaceID`): FaceID to get the name of.

            Returns:
                :py:class:`list`: List of vertex keys which make up the face.
        )delim"
    );


    SurfMeshCls.def("getCover",
        py::overload_cast<SurfaceMesh::SimplexID<1>>(&SurfaceMesh::get_cover<1>, py::const_),
        R"delim(
            Get the cover of a vertex.

            Args:
                SimplexID (:py:class:`VertexID`): VertexID to get the cover of.

            Returns:
                :py:class:`list`: Cover of the vertex.
        )delim"
    );


    SurfMeshCls.def("getCover",
        py::overload_cast<SurfaceMesh::SimplexID<2>>(&SurfaceMesh::get_cover<2>, py::const_),
        R"delim(
            Get the cover of an edge.

            Args:
                SimplexID (:py:class:`EdgeID`): EdgeID to get the cover of.

            Returns:
                :py:class:`list`: Cover of the edge.
        )delim"
    );

    /************************************
     *  UTILITY
     ************************************/
    SurfMeshCls.def_property_readonly("nVertices",
        py::overload_cast<>(&SurfaceMesh::size<1>, py::const_),
        R"delim(
            Get the number of vertices in the mesh.

            Returns:
                :py:class:`int`: Number of vertices.
        )delim"
    );


    SurfMeshCls.def_property_readonly("nEdges",
        py::overload_cast<>(&SurfaceMesh::size<2>, py::const_),
        R"delim(
            Get the number of edges in the mesh.

            Returns:
                :py:class:`int`: Number of edges.
        )delim"
    );


    SurfMeshCls.def_property_readonly("nFaces",
        py::overload_cast<>(&SurfaceMesh::size<3>, py::const_),
        R"delim(
            Get the number of faces in the mesh.

            Returns:
                :py:class:`int`: Number of faces.
        )delim"
    );

    SurfMeshCls.def("to_ndarray",
        [](const SurfaceMesh& mesh){
            std::map<typename SurfaceMesh::KeyType,typename SurfaceMesh::KeyType> sigma;

            double *vertices = new double[3*mesh.size<1>()];
            int    *edges = new int[2*mesh.size<2>()];
            int    *faces = new int[3*mesh.size<3>()];

            std::size_t i = 0;
            for(const auto vertexID : mesh.get_level_id<1>()){
                std::size_t o = 3*i;
                sigma[mesh.get_name(vertexID)[0]] = i++;
                auto vertex = vertexID.data();
                vertices[o]     = vertex[0];
                vertices[o+1]   = vertex[1];
                vertices[o+2]   = vertex[2];
            }

            i = 0;
            for(const auto edgeID : mesh.get_level_id<2>()){
                std::size_t o = 2*i;
                auto name = mesh.get_name(edgeID);

                edges[o]    = sigma[name[0]];
                edges[o+1]  = sigma[name[1]];
                ++i;
            }

            i = 0;
            for(const auto faceID : mesh.get_level_id<3>()){
                std::size_t o = 3*i;
                auto name = mesh.get_name(faceID);
                int orient = (*faceID).orientation;
                if(orient == 1){
                    faces[o]    = sigma[name[0]];
                    faces[o+1]  = sigma[name[1]];
                    faces[o+2]  = sigma[name[2]];
                }
                else{
                    faces[o]    = sigma[name[2]];
                    faces[o+1]  = sigma[name[1]];
                    faces[o+2]  = sigma[name[0]];
                }
                ++i;
            }

            auto free_vertices  = py::capsule(
                                    vertices,
                                    [](void *vertices) {
                                        delete[] reinterpret_cast<double*>(vertices);
                                 });
            auto free_edges     = py::capsule(
                                    edges,
                                    [](void *edges) {
                                        delete[] reinterpret_cast<int*>(edges);
                                  });
            auto free_faces     = py::capsule(
                                    faces,
                                    [](void *faces) {
                                        delete[] reinterpret_cast<int*>(faces);
                                  });

            return make_tuple(
                        py::array_t<double>(
                            std::array<std::size_t, 2>({mesh.size<1>(), 3}),
                            {3*sizeof(double), sizeof(double)},
                            vertices,
                            free_vertices
                        ),
                        py::array_t<int>(
                            std::array<std::size_t, 2>({mesh.size<2>(), 2}),
                            {2*sizeof(int), sizeof(int)},
                            edges,
                            free_edges
                        ),
                        py::array_t<int>(
                            std::array<std::size_t, 2>({mesh.size<3>(), 3}),
                            {3*sizeof(int), sizeof(int)},
                            faces,
                            free_faces
                        ));
        },
        R"delim(
            Converts the Surface Mesh into sets of numpy arrays.

            Returns:
                tuple(:py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`): Tuple of (nVertices, 3) array of vertex coordinates, (nEdges, 2) array of indices of vertices making up edges, (nFaces, 3) array of indices of vertices making up faces.
        )delim"
    );


    SurfMeshCls.def("onBoundary",
        py::overload_cast<const SurfaceMesh::SimplexID<1>>(&SurfaceMesh::onBoundary<1>, py::const_),
        R"delim(
            Check if a vertex is on a boundary.

            Args:
                SimplexID (:py:class:`VertexID`): VertexID to check.

            Returns:
                :py:class:`bool`: True if vertex is a member of an edge on a boundary.
        )delim"
    );

    SurfMeshCls.def("onBoundary",
        py::overload_cast<const SurfaceMesh::SimplexID<2>>(&SurfaceMesh::onBoundary<2>, py::const_),
        R"delim(
            Check if an edge is on a boundary.

            Args:
                SimplexID (:py:class:`EdgeID`): EdgeID to check.

            Returns:
                :py:class:`bool`: True if edge is a boundary.
        )delim"
    );

    SurfMeshCls.def("onBoundary",
        py::overload_cast<const SurfaceMesh::SimplexID<3>>(&SurfaceMesh::onBoundary<3>, py::const_),
        R"delim(
            Check if a face is on a boundary.

            Args:
                SimplexID (:py:class:`FaceID`): FaceID to check.

            Returns:
                :py:class:`bool`: True if face has an edge on a boundary.
        )delim"
    );

    SurfMeshCls.def("nearBoundary",
        py::overload_cast<const SurfaceMesh::SimplexID<1>>(&SurfaceMesh::nearBoundary<1>, py::const_),
        R"delim(
            Check if a vertex is near a boundary.

            Args:
                SimplexID (:py:class:`VertexID`): VertexID to check.

            Returns:
                :py:class:`bool`: True if vertex is a member of an edge on a boundary.
        )delim"
    );

    SurfMeshCls.def("nearBoundary",
        py::overload_cast<const SurfaceMesh::SimplexID<2>>(&SurfaceMesh::nearBoundary<2>, py::const_),
        R"delim(
            Check if an edge is near a boundary.

            Args:
                SimplexID (:py:class:`EdgeID`): EdgeID to check.

            Returns:
                :py:class:`bool`: True if edge is a boundary.
        )delim"
    );

    SurfMeshCls.def("nearBoundary",
        py::overload_cast<const SurfaceMesh::SimplexID<3>>(&SurfaceMesh::nearBoundary<3>, py::const_),
        R"delim(
            Check if a face is near a boundary.

            Args:
                SimplexID (:py:class:`FaceID`): FaceID to check.

            Returns:
                :py:class:`bool`: True if face has an edge on a boundary.
        )delim"
    );

    SurfMeshCls.def("splitSurfaces",
        &splitSurfaces,
        py::call_guard<py::scoped_ostream_redirect,
                py::scoped_estream_redirect>(),
        R"delim(
            Split disconnected surfaces into separate meshes.

            Args:
                SurfaceMesh (:py:class:`SurfaceMesh`): Surface mesh to split into multiple surfaces.

            Returns:
                :py:class:`list` (:py:class:`SurfaceMesh`): List of :py:class:`SurfaceMesh` containing separated surfaces.
        )delim"
    );


    SurfMeshCls.def("curvatureViaMDSB",
        [](const SurfaceMesh& mesh, std::size_t nIter){
            double* kh;
            double* kg;
            double* k1;
            double* k2;
            std::map<typename SurfaceMesh::KeyType,typename SurfaceMesh::KeyType> sigma;

            std::tie(kh,kg,k1,k2,sigma) = curvatureViaMDSB(mesh);

            auto free_kh  = py::capsule(
                                kh,
                                [](void *kh) {
                                    delete[] reinterpret_cast<double*>(kh);
                             });
            auto free_kg  = py::capsule(
                                kg,
                                [](void *kg) {
                                    delete[] reinterpret_cast<double*>(kg);
                             });
            auto free_k1  = py::capsule(
                                k1,
                                [](void *k1) {
                                    delete[] reinterpret_cast<double*>(k1);
                             });
            auto free_k2  = py::capsule(
                                k2,
                                [](void *k2) {
                                    delete[] reinterpret_cast<double*>(k2);
                             });

            if(nIter > 0){
                double* smoothedkh = new double[mesh.size<1>()];
                double* smoothedkg = new double[mesh.size<1>()];
                double* smoothedk1 = new double[mesh.size<1>()];
                double* smoothedk2 = new double[mesh.size<1>()];

                for(std::size_t round = 0; round < nIter; ++round){
                    // Smoothing
                    for(const auto vertexID : mesh.get_level_id<1>()) {
                        std::size_t i = sigma[vertexID.indices()[0]];
                        // double value = scale*kh[i];

                        double valuekh = kh[i];
                        double valuekg = kg[i];
                        double valuek1 = k1[i];
                        double valuek2 = k2[i];

                        auto neighbors = vertexID.cover();
                        for(std::size_t nbor : neighbors) {
                            valuekh += kh[sigma[nbor]];
                            valuekg += kg[sigma[nbor]];
                            valuek1 += k1[sigma[nbor]];
                            valuek2 += k2[sigma[nbor]];
                            // value += (1-scale)/static_cast<double>(neighbors.size())*kh[sigma[nbor]];
                        }
                        // smoothed[i] = value;
                        smoothedkh[i] = valuekh/static_cast<double>(neighbors.size()+1);
                        smoothedkg[i] = valuekg/static_cast<double>(neighbors.size()+1);
                        smoothedk1[i] = valuek1/static_cast<double>(neighbors.size()+1);
                        smoothedk2[i] = valuek2/static_cast<double>(neighbors.size()+1);
                    }
                    std::copy(smoothedkh, smoothedkh+mesh.size<1>(), kh);
                    std::copy(smoothedkg, smoothedkg+mesh.size<1>(), kg);
                    std::copy(smoothedk1, smoothedk1+mesh.size<1>(), k1);
                    std::copy(smoothedk2, smoothedk2+mesh.size<1>(), k2);
                }

                // // Smooth Gaussian curvature
                // // double scale = 0.7;
                // for(std::size_t round = 0; round < nIter; ++round){
                //     // Smoothing
                //     for(const auto vertexID : mesh.get_level_id<1>()) {
                //         std::size_t i = sigma[vertexID.indices()[0]];
                //         // double value = scale*kg[i];
                //         double value = kg[i];
                //         auto neighbors = vertexID.cover();
                //         for(std::size_t nbor : neighbors) {
                //             value += kg[sigma[nbor]];
                //             // value += (1-scale)/static_cast<double>(neighbors.size())*kg[sigma[nbor]];
                //         }
                //         // smoothed[i] = value;
                //         smoothed[i] = value/static_cast<double>(neighbors.size()+1);
                //     }
                //     for(std::size_t i = 0; i < mesh.size<1>(); ++i){
                //         kg[i] = smoothed[i];
                //     }
                // }

                delete[] smoothedkh;
                delete[] smoothedkg;
                delete[] smoothedk1;
                delete[] smoothedk2;
            }

            return  std::make_tuple(
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            kh,
                            free_kh),
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            kg,
                            free_kg),
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            k1,
                            free_k1),
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            k2,
                            free_k2)
                    );
        },
        R"delim(
            Compute the mean, Gaussian, and principal curvatures of the mesh.

            Args:
                nIter (:py:class:`int`): Number of smoothing iterations to run

            Returns:
                tuple(:py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`): Tuple of arrays containing Mean, Gaussian, First Prinicipal, and Second Principal curvatures.
        )delim"
    );

    SurfMeshCls.def("curvatureViaJets",
        [](const SurfaceMesh& mesh, std::size_t nIter){
            double* kh;
            double* kg;
            double* k1;
            double* k2;
            std::map<typename SurfaceMesh::KeyType,typename SurfaceMesh::KeyType> sigma;

            std::tie(kh,kg,k1,k2,sigma) = curvatureViaJets(mesh, 2, 2);

            auto free_kh  = py::capsule(
                                kh,
                                [](void *kh) {
                                    delete[] reinterpret_cast<double*>(kh);
                             });
            auto free_kg  = py::capsule(
                                kg,
                                [](void *kg) {
                                    delete[] reinterpret_cast<double*>(kg);
                             });
            auto free_k1  = py::capsule(
                                k1,
                                [](void *k1) {
                                    delete[] reinterpret_cast<double*>(k1);
                             });
            auto free_k2  = py::capsule(
                                k2,
                                [](void *k2) {
                                    delete[] reinterpret_cast<double*>(k2);
                             });

           if(nIter > 0){
                double* smoothedkh = new double[mesh.size<1>()];
                double* smoothedkg = new double[mesh.size<1>()];
                double* smoothedk1 = new double[mesh.size<1>()];
                double* smoothedk2 = new double[mesh.size<1>()];

                for(std::size_t round = 0; round < nIter; ++round){
                    // Smoothing
                    for(const auto vertexID : mesh.get_level_id<1>()) {
                        std::size_t i = sigma[vertexID.indices()[0]];
                        // double value = scale*kh[i];

                        double valuekh = kh[i];
                        double valuekg = kg[i];
                        double valuek1 = k1[i];
                        double valuek2 = k2[i];

                        auto neighbors = vertexID.cover();
                        for(std::size_t nbor : neighbors) {
                            valuekh += kh[sigma[nbor]];
                            valuekg += kg[sigma[nbor]];
                            valuek1 += k1[sigma[nbor]];
                            valuek2 += k2[sigma[nbor]];
                            // value += (1-scale)/static_cast<double>(neighbors.size())*kh[sigma[nbor]];
                        }
                        // smoothed[i] = value;
                        smoothedkh[i] = valuekh/static_cast<double>(neighbors.size()+1);
                        smoothedkg[i] = valuekg/static_cast<double>(neighbors.size()+1);
                        smoothedk1[i] = valuek1/static_cast<double>(neighbors.size()+1);
                        smoothedk2[i] = valuek2/static_cast<double>(neighbors.size()+1);
                    }
                    std::copy(smoothedkh, smoothedkh+mesh.size<1>(), kh);
                    std::copy(smoothedkg, smoothedkg+mesh.size<1>(), kg);
                    std::copy(smoothedk1, smoothedk1+mesh.size<1>(), k1);
                    std::copy(smoothedk2, smoothedk2+mesh.size<1>(), k2);
                }

                // // Smooth Gaussian curvature
                // // double scale = 0.7;
                // for(std::size_t round = 0; round < nIter; ++round){
                //     // Smoothing
                //     for(const auto vertexID : mesh.get_level_id<1>()) {
                //         std::size_t i = sigma[vertexID.indices()[0]];
                //         // double value = scale*kg[i];
                //         double value = kg[i];
                //         auto neighbors = vertexID.cover();
                //         for(std::size_t nbor : neighbors) {
                //             value += kg[sigma[nbor]];
                //             // value += (1-scale)/static_cast<double>(neighbors.size())*kg[sigma[nbor]];
                //         }
                //         // smoothed[i] = value;
                //         smoothed[i] = value/static_cast<double>(neighbors.size()+1);
                //     }
                //     for(std::size_t i = 0; i < mesh.size<1>(); ++i){
                //         kg[i] = smoothed[i];
                //     }
                // }

                delete[] smoothedkh;
                delete[] smoothedkg;
                delete[] smoothedk1;
                delete[] smoothedk2;
            }

            return  std::make_tuple(
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            kh,
                            free_kh),
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            kg,
                            free_kg),
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            k1,
                            free_k1),
                        py::array_t<double>(
                            std::array<std::size_t, 1>({mesh.size<1>()}),
                            {sizeof(double)},
                            k2,
                            free_k2)
                    );
        },
        R"delim(
            Compute the mean, Gaussian, and principal curvatures of the mesh.

            Args:
                nIter (:py:class:`int`): Number of smoothing iterations to run

            Returns:
                tuple(:py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`): Tuple of arrays containing Mean, Gaussian, First Prinicipal, and Second Principal curvatures.
        )delim"
    );


    SurfMeshCls.def("getNormal",
        py::overload_cast<const SurfaceMesh&, SurfaceMesh::SimplexID<3>>(&getNormal),
        R"delim(
            Return the normal of the face.

            Args:
                mesh  (:py:class:`SurfaceMesh`): Mesh of interest.

                FaceID (:py:class:`FaceID`): FaceID of interest.

            Returns:
                :py:class:`Vector`: Vector normal to the face.
        )delim"
    );


    SurfMeshCls.def("getNormal",
        py::overload_cast<const SurfaceMesh&, SurfaceMesh::SimplexID<1>>(&getNormal),
        R"delim(
            Return the normal of a vertex.

            This is computed as the average normal of all incident faces.

            Args:
                mesh  (:py:class:`SurfaceMesh`): Mesh of interest.

                VertexID (:py:class:`VertexID`): VertexID of interest.

            Returns:
                :py:class:`Vector`: Vector normal to the vertex.
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

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`VertexID`.
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

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`VertexID`.
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

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`EdgeID`.
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

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`EdgeID`.
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

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`FaceID`.
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

            Returns:
                :ref:`Iterator <python:iterator-objects>`: An iterator over all :py:class:`FaceID`.
        )delim"
    );


    /************************************
     *  FUNCTIONS
     ************************************/
    SurfMeshCls.def("smooth", &smoothMesh,
        py::arg("max_iter")=6, py::arg("preserve_ridges")=false, py::arg("rings")=2, py::arg("verbose")=false,
        py::call_guard<py::scoped_ostream_redirect,
                py::scoped_estream_redirect>(),
        R"delim(
            Perform mesh smoothing.

            This operation performs an iterative weightedVertexSmooth
            and edge flipping approach.

            Args:
                maxIter (int): Maximum number of smoothing iterations.
                preserveRidges (bool):  Prevent flipping of edges along ridges.
                rings (int): Number of LST rings to consider.
                verbose (bool): Print details.
        )delim"
    );


    SurfMeshCls.def("coarse", &coarse,
        py::arg("rate"), py::arg("flatRate"), py::arg("denseWeight"), py::arg("rings")=2, py::arg("verbose")=false,
        R"delim(
            Coarsen a surface mesh.

            This operation selects vertices to decimate then removes them and
            retriangulates the resulting hole.

            Args:
                rate (float): Threshold value for coarsening.
                flatRate (float): Priority of decimating flat regions.
                denseWeight (float): Priority of decimating dense regions.
                rings (int): Number of LST rings to consider.
                verbose (bool): Print details.
        )delim"
    );


    SurfMeshCls.def("coarse_flat",
        [](SurfaceMesh& mesh, double rate, int niter, std::size_t rings, bool verbose){
            for(int i = 0; i < niter; ++i) coarse_flat(mesh, rate, 0.5, rings, verbose);
        },
        py::arg("rate")=0.016, py::arg("numiter")=1, py::arg("rings")=2, py::arg("verbose")=false,
        R"delim(
            Coarsen flat regions of a surface mesh.

            Args:
                rate (float): Threshold value.
                numiter (int): Number of iterations to run.
                rings (int): Number of LST rings to consider.
                verbose (bool): Print details.
        )delim"
    );


    SurfMeshCls.def("coarse_dense",
        [](SurfaceMesh& mesh, double rate, int niter, std::size_t rings, bool verbose){
            for(int i = 0; i < niter; ++i) coarse_dense(mesh, rate, 10, rings, verbose);
        },
        py::arg("rate")=1.6, py::arg("numiter")=1, py::arg("rings")=2, py::arg("verbose")=false,
        R"delim(
            Coarsen dense regions of a surface mesh.

            Args:
                rate (float): Threshold value.
                numiter (int): Number of iterations to run.
                rings (int): Number of LST rings to consider.
                verbose (bool): Print details.
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


    SurfMeshCls.def("getVolume", &getVolume,
        R"delim(
            Get the volume of the surface mesh.

            Returns:
                :py:class:`float`: Volume of the mesh.
        )delim"
    );

    SurfMeshCls.def("getCenterRadius", &getCenterRadius,
        R"delim(
            Get the center and radius of a surface mesh.

            Returns:
                :py:class:`tuple` (:py:class:`Vector`, :py:class:`float`): Center of the mesh and radius.
        )delim"
    );

    SurfMeshCls.def("translate",
        py::overload_cast<SurfaceMesh&, Vector>(&translate),
        R"delim(
            Translate the mesh.

           Args:
                vec (:py:class:`Vector`): Vector translation.
        )delim"
    );

    SurfMeshCls.def("scale", py::overload_cast<SurfaceMesh&, double>(&scale),
        R"delim(
            Scale mesh by a scale factor.

            Args:
                scale (:py:class:`float`): scale factor.
        )delim"
    );

    /************************************
     *  ORIENTATION
     ************************************/
    SurfMeshCls.def("init_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::init_orientation<SurfaceMesh>),
        R"delim(
            Initialize the orientations of all mesh simplices.

            This function must be called prior to definition of face orientations.
        )delim"
    );


    SurfMeshCls.def("check_orientation",
        py::overload_cast<SurfaceMesh&>(&casc::check_orientation<SurfaceMesh>),
        R"delim(
            Check consistency and assign Face orientations.

            :py:func:`SurfaceMesh.init_orientation` must be called before this
            function can operate.


            Returns:
                :py:class:`tuple` (:py:class:`int`, :py:class:`bool`, :py:class:`bool`): Tuple containing: number of connected components, orientability status, and manifoldness.
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
            Compute the orientations of the mesh.

            This is equivalent to the following series of calls:

            >>> init_orientation(mesh)
            >>> clear_orientation(mesh)
            >>> check_orientation(mesh)

            Returns:
                :py:class:`tuple` (:py:class:`int`, :py:class:`bool`, :py:class:`bool`): Tuple containing: number of connected components, orientability status, and manifoldness.
        )delim"
    );


    SurfMeshCls.def("flipNormals", &flipNormals,
        R"delim(
            Flip the orientation of face normals.
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
            Set normals such that mesh volume is positive.
        )delim"
    );

    SurfMeshCls.def("getBettiNumbers", &getBettiNumbers,
        R"delim(
            Compute the Betti numbers of the mesh.
        )delim"
    );
}

} // end namespace gamer