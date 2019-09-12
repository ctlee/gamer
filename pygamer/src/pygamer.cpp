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
#include <pybind11/iostream.h>

#include "gamer/SurfaceMesh.h"
#include "gamer/TetMesh.h"
#include "gamer/PDBReader.h"
#include "gamer/version.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

// Forward function declarations
void init_Vector(py::module &);
void init_SMGlobal(py::module &);
void init_SMVertex(py::module &);
void init_SMEdge(py::module &);
void init_SMFace(py::module &);
void init_SMSimplexID(py::module &);
void init_SMFunctions(py::module &);
void init_SurfaceMesh(py::module &);

void init_TMGlobal(py::module &);
void init_TMVertex(py::module &);
void init_TMEdge(py::module &);
void init_TMFace(py::module &);
void init_TMCell(py::module &);
void init_TMSimplexID(py::module &);
void init_TetMesh(py::module &);

// Initialize the main `pygamer` module
PYBIND11_MODULE(pygamer, pygamer) {
    pygamer.doc() = "Python wrapper around the GAMer C++ library.";

    init_Vector(pygamer);           // Vector class

    /************************************
     *  SURFACEMESH SUBMODULE
     ************************************/
    py::module SurfMeshMod = pygamer.def_submodule("surfacemesh",
        "Submodule containing SurfaceMesh along with related objects and methods.");

    // WARNING: the order of initialization matters for dependent calls.
    init_SMGlobal(SurfMeshMod);     // Global class
    init_SMVertex(SurfMeshMod);     // Vertex class
    init_SMEdge(SurfMeshMod);       // Edge class
    init_SMFace(SurfMeshMod);       // Face class
    init_SMSimplexID(SurfMeshMod);  // SurfaceMesh::SimplexID class
    init_SMFunctions(SurfMeshMod);  // Functions pyg.sm.*
    init_SurfaceMesh(SurfMeshMod);  // SurfaceMesh class


    /************************************
     *  TETMESH SUBMODULE
     ************************************/
    py::module TetMeshMod = pygamer.def_submodule("tetmesh",
        "Submodule containing TetMesh along with related objects and methods.");

    // WARNING: the order of initialization matters for dependent calls.
    init_TMVertex(TetMeshMod);     // Vertex class
    init_TMEdge(TetMeshMod);       // Edge class
    init_TMFace(TetMeshMod);       // Face class
    init_TMCell(TetMeshMod);       // Cell class
    init_TMSimplexID(TetMeshMod);  // TetMesh::SimplexID class
    init_TetMesh(TetMeshMod);      // TetMesh class

    /************************************
     *  PYGAMER FUNC/OBJECT DEFS
     ************************************/
    pygamer.def("readOFF", &readOFF,
        py::arg("filename"),
        R"delim(
            Read OFF file to mesh

            Args:
                filename (:py:class:`str`): Filename to read from.

            Returns:
                :py:class:`surfacemesh.SurfaceMesh`: Mesh of interest/
        )delim"
    );


    pygamer.def("readPDB_molsurf", &readPDB_molsurf,
        py::arg("filename"),
        R"delim(
            Read a PDB file into a mesh

            Args:
                filename (:py:class:`str`): PDB file to read.

            Returns:
                :py:class:`surfacemesh.SurfaceMesh`: Meshed object.
        )delim"
    );

    pygamer.def("readPDB_gauss", &readPDB_gauss,
        py::arg("filename"),
        py::arg("blobbyness") = -0.2,
        py::arg("isovalue") = 2.5,
        R"delim(
            Read a PDB file into a mesh

            Args:
                filename (:py:class:`str`): PDB file to read.
                blobbyness (:py:class:`float`): Blobbiness of the Gaussian.
                isovalue (:py:class:`float`): The isocontour value to mesh.

            Returns:
                :py:class:`surfacemesh.SurfaceMesh`: Meshed object.
        )delim"
    );


    pygamer.def("readPQR_molsurf", &readPQR_molsurf,
        py::arg("filename"),
        R"delim(
            Read a PQR file into a mesh

            Args:
                filename (:py:class:`str`): PQR file to read
            Returns:
                :py:class:`surfacemesh.SurfaceMesh`: Meshed object
        )delim"
    );


    pygamer.def("readPQR_gauss", &readPQR_gauss,
        py::arg("filename"),
        py::arg("blobbyness") = -0.2,
        py::arg("isovalue") = 2.5,
        R"delim(
            Read a PQR file into a mesh

            Args:
                filename (:py:class:`str`): PQR file to read.
                blobbyness (:py:class:`float`): Blobbiness of the Gaussian.
                isovalue (:py:class:`float`): The isocontour value to mesh.

            Returns:
                :py:class:`surfacemesh.SurfaceMesh`: Meshed object.
        )delim"
    );


    pygamer.def("writeOFF", py::overload_cast<const std::string&, const SurfaceMesh&>(&writeOFF),
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in OFF format

            Args:
                filename (:py:class:`str`): Filename to write to.
                mesh (:py:class:`surfacemesh.SurfaceMesh`): Mesh of interest.
        )delim"
    );


    pygamer.def("writeOFF", py::overload_cast<const std::string&, const TetMesh&>(&writeOFF),
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in OFF format

            Args:
                filename (:py:class:`str`): Filename to write to.
                mesh (:py:class:`tetmesh.TetMesh`): Mesh of interest.
        )delim"
    );


    pygamer.def("readOBJ", &readOBJ,
        py::arg("filename"),
        R"delim(
            Read OBJ file to mesh

            Args:
                filename (:py:class:`str`): Filename to read from
            Returns:
                :py:class:`surfacemesh.SurfaceMesh`: Mesh of interest
        )delim"
    );


    pygamer.def("writeOBJ", &writeOBJ,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in OBJ format

            Args:
                filename (:py:class:`str`): Filename to write to
                mesh (:py:class:`surfacemesh.SurfaceMesh`): Mesh of interest
        )delim"
    );


    pygamer.def("writeVTK", &writeVTK,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in VTK format

            Args:
                filename (:py:class:`str`): Filename to write to
                mesh (:py:class:`tetmesh.TetMesh`): Mesh of interest
        )delim"
    );


    pygamer.def("writeDolfin", &writeDolfin,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in Dolfin XML format

            Args:
                filename (:py:class:`str`): Filename to write to
                mesh (:py:class:`tetmesh.TetMesh`): Mesh of interest
        )delim"
    );


    pygamer.def("writeTriangle", &writeTriangle,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Write mesh to file in Triangle format

            Args:
                filename (:py:class:`str`): Filename to write to
                mesh (:py:class:`tetmesh.TetMesh`): Mesh of interest
        )delim"
    );

    pygamer.def("printQualityInfo", &printQualityInfo,
        py::arg("filename"), py::arg("mesh"),
        R"delim(
            Print mesh quality info to file

            Args:
                filename (:py:class:`str`): Filename to write to
                mesh (:py:class:`surfacemesh.SurfaceMesh`): Mesh of interest
        )delim"
    );


    pygamer.def("readDolfin", &readDolfin,
        py::arg("filename"),
        R"delim(
            Read Dolfin format mesh into mesh

            Args:
                filename (:py:class:`str`): Filename to write to

            Returns:
                :py:class:`tetmesh.TetMesh`: Tetrahedral mesh
        )delim"
    );

    pygamer.def("makeTetMesh", &makeTetMesh,
        py::call_guard<py::scoped_ostream_redirect,
                py::scoped_estream_redirect>(),
        R"delim(
            Call tetgen to make a TetMesh

            Args:
                meshes (:py:class:`list`(:py:class:`surfacemesh.SurfaceMesh`): List of meshes with filled metadata

            Returns:
                :py:class:`tetmesh.TetMesh`: Resulting tetrahedral mesh
        )delim"
    );

    pygamer.def("__version__",
        [](){
            extern const std::string gVERSION;
            return gVERSION;
        },
        R"delim(
            Get the full version of PyGAMer

            Returns:
                :py:class:`str`: Long version of PyGAMer
        )delim"
    );

    pygamer.def("getVersion",
        [](){
            extern const std::string gVERSION_SHORT;
            return gVERSION_SHORT;
        },
        R"delim(
            Get the short version of PyGAMer

            Returns:
                :py:class:`str`: Short version of PyGAMer
        )delim"
    );
}

} // end namespace gamer