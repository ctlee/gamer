/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2018
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

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

#include <gamer.h>

#define TETLIBRARY
#include <libraries/tetgen/tetgen.h>

#include <libraries/casc/casc>

#include "Vertex.h"
#include "SurfaceMesh.h"

#ifndef SWIG
namespace tetmesh
{
/**
 * @brief      Properties that Faces should have
 */
struct FaceProperties
{
    int  marker;   /**< @brief Marker */
    // bool selected; /**< @brief Selection flag */
};

/**
 * @brief      Face object
 */
struct Face : FaceProperties
{
    /// Default constructor
    Face() {}
    Face(int marker) : Face(FaceProperties{marker}) {}
    Face(FaceProperties prop) : FaceProperties(prop) {}
};

struct CellProperties
{
    int marker;
    // int grpID;
    // int material;
};

struct Cell : casc::Orientable, CellProperties
{
    Cell() {}
    Cell(int orient, int marker) : Cell(Orientable{orient}, CellProperties{marker}) {}
    Cell(Orientable orient, CellProperties prop)
        : Orientable(orient), CellProperties(prop)
    {}
};

struct Edge : Vertex
{
    Edge() {}
    Edge(Vertex v) : Vertex(v) {}
};

/**
 * @brief      Type for containing root metadata
 */
struct Global
{
    bool isHigher_order() const;

    bool higher_order;
};

struct TetVertex : Vertex
{
    using Vertex::Vertex; // Inherit constructor from Vertex
    TetVertex() {}
    TetVertex(Vertex v) : Vertex(v) {}

    double getError() const;

    static double error;


    bool operator< (const TetVertex& rhs)
    {
        return  error < rhs.getError();
    };

    bool operator> (const TetVertex& rhs)
    {
        return  error > rhs.getError();
    };

};

struct TetEdge : TetVertex {
    TetEdge() {}
    TetEdge(TetVertex v) : TetVertex(v) {}

};

/**
 * @brief      A helper struct containing the traits/types in the simplicial
 *             complex
 */
struct complex_traits
{
    /// The index type
    using KeyType = int;
    /// The types of each node
    using NodeTypes = util::type_holder<Global, TetVertex, Edge, Face, Cell>;
    /// The types of each edge
    using EdgeTypes = util::type_holder<casc::Orientable, casc::Orientable, casc::Orientable, casc::Orientable>;
};
} // END namespace tetmesh
#endif //SWIG


using TetMesh = casc::simplicial_complex<tetmesh::complex_traits>;


#ifndef SWIG
/// Forward class declaration
class tetgenio;

std::unique_ptr<TetMesh> tetgenioToTetMesh(tetgenio &tetio);

std::unique_ptr<TetMesh> makeTetMesh(
        const std::vector<SurfaceMesh*> &surfmeshes,
        std::string tetgen_params);
#endif //SWIG

template <std:: size_t level>
TetMesh::SimplexID<level> get_lowest_err();

void smoothMesh(TetMesh & mesh);
void writeVTK(const std::string& filename, const TetMesh &mesh);
void writeOFF(const std::string& filename, const TetMesh &mesh);
void writeDolfin(const std::string &filename, const TetMesh &mesh);
void writeTriangle(const std::string &filename, const TetMesh &mesh);

// void writeMCSF(const std::string &filename, const TetMesh &mesh);
//void writeDiffPack
//void writeCARP
