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

void init_SMFunctions(py::module& mod){
    mod.def("readOFF", &readOFF);
    mod.def("writeOFF", &writeOFF);
    mod.def("readOBJ", &readOBJ);
    mod.def("writeOBJ", &writeOBJ);
    mod.def("getMeanCurvature", &getMeanCurvature);
    mod.def("getGaussianCurvature", &getGaussianCurvature,
        py::arg("mesh"), py::arg("VertexID"),
        R"delim(
            Compute the Gaussian curvature around a vertex

            Args:
                mesh (:py:class:`SurfaceMesh`): Mesh of interest
                VertexID (:py:class:`VertexID`): SimplexID of the vertex

            Returns:
                double: Local Gaussian curvature
        )delim"
        );
    mod.def("smoothMesh", &smoothMesh,
        py::arg("mesh"), py::arg("maxIter"), py::arg("preserveRidges"), py::arg("verbose"),
        R"delim(
            Smooth a :py:class:`SurfaceMesh`

            This operation performs an iterative weightedVertexSmooth
            and edge flipping approach.

            Args:
                mesh (:py:class:`SurfaceMesh`): Surface mesh to smooth
                maxIter (int): Maximum number of smoothing iterations
                preserveRidges (bool):  Prevent flipping of edges along ridges.
                verbose (bool): Print details to std::out
        )delim"
        );
    mod.def("coarse", &coarse,
        py::arg("mesh"), py::arg("coarseRate"), py::arg("flatRate"),
        py::arg("denseWeight"),
        R"delim(
            Coarsen a surface mesh

            This operation selects vertices to decimate then removes them and
            retriangulates the resulting hole.

            Args:
                mesh (:py:class:`SurfaceMesh`): Surface mesh to coarsen
                coarseRate (double): Threshold value for coarsening
                flatRate (double): Priority of decimating flat regions
                denseWeight (double): Priority of decimating dense regions
        )delim"
        );
    mod.def("normalSmooth", &normalSmooth,
        py::arg("mesh"),
        R"delim(
            Perform smoothing of mesh face normals

            Args:
                mesh (:py:class:`SurfaceMesh`): SurfaceMesh to smooth
        )delim"
        );
    mod.def("cube", &cube,
        py::arg("order"),
        R"delim(
            Construct a triangulated cube mesh

            Args:
                order (int): Number of subdivisions

            Returns:
                mesh: Cubical :py:class:`SurfaceMesh`
        )delim"
        );
    mod.def("sphere", &sphere,
        py::arg("order"),
        R"delim(
            Construct a triangulated icosahedral sphere mesh

            Args:
                order (int): Number of subdivisions

            Returns:
                mesh: Spherical :py:class:`SurfaceMesh`
        )delim"
        );
    mod.def("print", &print,
        py::arg("mesh"),
        R"delim(
            Print a surface mesh.

            Args:
                mesh (:py:class:`SurfaceMesh`): SurfaceMesh to print
        )delim"
        );
}