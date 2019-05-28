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

#include <sstream>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

#include <gamer.h>

#define TETLIBRARY
#include <tetgen.h>

#include <casc/casc>

#include "Vertex.h"
#include "SurfaceMesh.h"

namespace tetmesh
{
/**
 * @brief      Type for containing root metadata
 */
struct Global
{
    bool higher_order;
};

struct VertexProperties
{
    double error;

    VertexProperties() : error(-1) {}
    VertexProperties(double error) : error(error) {}
};

struct TetVertex : Vertex, VertexProperties
{
    TetVertex(): TetVertex(Vertex(),VertexProperties()) {}
    template<typename... Args>
    TetVertex(Args&&... args) : TetVertex(Vertex(std::forward<Args>(args)...)) {}
    TetVertex(Vertex v) : TetVertex(v, VertexProperties(-1)) {}
    TetVertex(Vertex v, VertexProperties p) : Vertex(v), VertexProperties(p) {}

    /**
     * @brief      Operator<< overload
     *
     * @param      output  stream to print to
     * @param[in]  v       Vertex to print
     *
     * @return     the stream
     */
    friend std::ostream& operator<<(std::ostream& output, const TetVertex& v){
        output  << "tetmesh::TetVertex(x:" << v[0]
                << ",y:" << v[1]
                << ",z:" << v[2]
                << ";m:" << v.marker
                << ";sel:" << v.selected
                << ";err:" << v.error
                << ")";
        return output;
    }

    /**
     * @brief      Returns a string representation of the object.
     *
     * @return     String representation of the object.
     */
    std::string to_string() const{
        std::ostringstream output;
        output  << *this;
        return output.str();
    }
};


struct Edge : Vertex
{
    using Vertex::Vertex;

    /**
     * @brief      Operator<< overload
     *
     * @param      output  stream to print to
     * @param[in]  v       Vertex to print
     *
     * @return     the stream
     */
    friend std::ostream& operator<<(std::ostream& output, const Edge& v){
        output  << "tetmesh::Edge(x:" << v[0]
                << ",y:" << v[1]
                << ",z:" << v[2]
                << ";m:" << v.marker
                << ";sel:" << v.selected
                << ")";
        return output;
    }

    /**
     * @brief      Returns a string representation of the object.
     *
     * @return     String representation of the object.
     */
    std::string to_string() const{
        std::ostringstream output;
        output  << *this;
        return output.str();
    }
};

/**
 * @brief      Properties that Faces should have
 */
struct FaceProperties
{
    int  marker;   /**< @brief Marker */
    bool selected;
};

/**
 * @brief      Face object
 */
struct Face : FaceProperties
{
    /// Default constructor
    Face() : Face(FaceProperties{-1, false}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  marker    The marker
     */
    Face(int marker, bool selected) : Face(FaceProperties{marker, selected}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  prop    Properties of a face
     */
    Face(FaceProperties prop) : FaceProperties(prop) {}

    friend std::ostream& operator<<(std::ostream& output, const Face& f){
        output  << "tetmesh::Face("
                << "m:" << f.marker
                << ";sel:" << f.selected << ")";
        return output;
    }

    std::string to_string() const{
        std::ostringstream output;
        output  << *this;
        return output.str();
    }
};

struct CellProperties
{
    int marker;
    bool selected;
};

struct Cell : casc::Orientable, CellProperties
{
    Cell() : Cell(-1,false) {}
    Cell(int marker, bool selected) : Cell(0, marker, selected) {}
    Cell(int orient, int marker, bool selected) : Cell(Orientable{orient}, CellProperties{marker, selected}) {}
    Cell(Orientable orient, CellProperties prop)
        : Orientable(orient), CellProperties(prop)
    {}

    friend std::ostream& operator<<(std::ostream& output, const Cell& c){
        output  << "tetmesh::Cell("
                << "m:" << c.marker
                << ";sel:" << std::boolalpha << c.selected
                << ";o:" << c.orientation << ")";
        return output;
    }

    std::string to_string() const{
        std::ostringstream output;
        output  << *this;
        return output.str();
    }
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


using TetMesh = casc::simplicial_complex<tetmesh::complex_traits>;


/// Forward class declaration
class tetgenio;

std::unique_ptr<TetMesh> tetgenioToTetMesh(tetgenio &tetio);

std::unique_ptr<TetMesh> makeTetMesh(
        const std::vector<SurfaceMesh*> &surfmeshes,
        std::string tetgen_params);


void smoothMesh(TetMesh & mesh);
void writeVTK(const std::string& filename, const TetMesh &mesh);
void writeOFF(const std::string& filename, const TetMesh &mesh);
void writeDolfin(const std::string &filename, const TetMesh &mesh);
void writeTriangle(const std::string &filename, const TetMesh &mesh);

//void writeMCSF(const std::string &filename, const TetMesh &mesh);
//void writeDiffPack
//void writeCARP

std::unique_ptr<TetMesh> readDolfin(const std::string&filename);
