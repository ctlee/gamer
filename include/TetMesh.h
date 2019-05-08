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

struct ErrorProperties
{
    double error;
};

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
struct Face : FaceProperties, ErrorProperties
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

struct Cell : casc::Orientable, CellProperties, ErrorProperties
{
    Cell() {}
    Cell(int orient, int marker) : Cell(Orientable{orient}, CellProperties{marker}) {}
    Cell(Orientable orient, CellProperties prop)
        : Orientable(orient), CellProperties(prop)
    {}
};

struct Edge : Vertex, ErrorProperties
{
    Edge() {}
    Edge(Vertex v) : Vertex(v) {}

    bool operator< (const Edge& rhs)
    {
        return  error < rhs.error;
    };

    bool operator> (const Edge& rhs)
    {
        return  error > rhs.error;
    };
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

    double error;


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

// Core decimation operators
void vertexRemoval(tetmesh::TetVertex v) {
    // Remove v

    // Retriangulate hole
}

template <typename Complex, template <typename> class Callback>
void edgeCollapse(TetMesh & mesh, TetMesh::SimplexID<2> edge, double vertexLoc, Callback<Complex> &&clbk)
{
    auto name = mesh.get_name(edge);
    auto v = *mesh.get_simplex_down(edge, name[0])
             - *mesh.get_simplex_down(edge, name[1]);
    Vector pos = mesh.get_simplex_down(edge, name[0]).data().position - v*vertexLoc;

    using SimplexMap = typename casc::SimplexMap<Complex>;

    casc::SimplexMap<Complex> simplexMap;
    int np = casc::decimateFirstHalf(mesh, edge, simplexMap);
    typename casc::decimation_detail::SimplexDataSet<Complex>::type rv;

    run_user_callback(mesh, simplexMap,  std::forward<Callback<Complex>>(clbk), rv);

    tetmesh::TetVertex vert;
    for (auto vpair : std::get<1>(rv)) {
        auto vnames = vpair.first;
        for (auto name : vnames) {
            if (np == name) {
                vert = vpair.second;
                goto endloop;
            }

        }
    }
    endloop:
    vert.position = pos;
    casc::decimateBackHalf(mesh, simplexMap, rv);
}


TetMesh::SimplexID<2> getLowestErr(TetMesh &complex);

template <typename Complex, template <typename> class Callback>
void decimated(TetMesh & mesh, double threshold, Callback<Complex> &&cbk){
    if (threshold < 1) {
        throw std::invalid_argument("Threshold to decimate to must be at least 1");
    }

    // Exit, mesh already sufficiently decimated
    if (threshold > mesh.size<2>()) {
        std::cout << "Decimation finished to requested threshold: " << threshold <<"; final mesh size: "
        << mesh.size<2>() << std::endl;
        return;
    }

    std::vector<TetMesh::SimplexID<2>> sortedEdges;
    // initialize vector of edges
    for (TetMesh::SimplexID<2> e : mesh.get_level_id<2>()) {
        sortedEdges.push_back(e);
    }

    auto edgeErrorCmp = [](const TetMesh::SimplexID<2> e1, const TetMesh::SimplexID<2> e2) -> bool{
        return *e1 < *e2;
    }

    // sort edges in order of error
    std::sort_heap(sortedEdges.begin(), sortedEdges.end(), edgeErrorCmp);

    // remove invalid edges (if edge with both vertices decimated already decimated, will cause segfault)
    for (auto it = sortedEdges.begin(); it != sortedEdges.end();) {
        auto edge = *it;
        auto name =  mesh.get_name(edge);

        if (encountered[name[0]] || encountered[name[1]]) {
            it = sortedEdges.erase(it);
        } else {
            ++it;
        }
        encountered[name[0]] = true;
        encountered[name[1]] = true;
    }

    for (int i =0 ; i < 100; i++ ) {
        std::cout << i << ": " << sortedEdges[i] << std::endl;
    }


    int initialLevelSimplexCount = mesh.size<2>();
    int numToRemove = initialLevelSimplexCount - threshold;

    sortedEdges.erase(sortedEdges.begin() + 8);
    sortedEdges.erase(sortedEdges.begin() + 14);
    for(auto it = sortedEdges.begin(); it != sortedEdges.end(); it++){
        if (numToRemove == 0) {
            break;
        }

        if (*it != NULL) {
            numToRemove--;
            std::cout << "collapse #" << numToRemove << ": edge" << *it << std::endl;
            edgeCollapse(mesh, *it, 0.5, Callback<TetMesh>());
            std::cout << "collapsed" << std::endl;
        }

    }
}



void markDecimatedEdge(TetMesh & mesh, std::vector<std::pair<std::pair<int,int>, double>> list, TetMesh::SimplexID<2> edge,
                       double loc) {
    auto name =  mesh.get_name(edge);
    list.push_back(std::pair<std::pair<int,int>,double>(std::pair<int,int>(name[0], name[1]), loc));
}

void writeRestrictionMatrix(const std::string &filename, int origSize,
        std::vector<std::pair<std::pair<int,int>, double>> list) {

    std::ofstream fout(filename);
    if(!fout.is_open())
    {
        std::cerr   << "File '" << filename
                    << "' could not be writen to." << std::endl;
        exit(1);
    }

    /*
     * //Write identity matrix
    for (int i = 0; i < (origSize-list.size()); i++) {
        for (int j = 0; j < origSize; i++) {
            if (i == j) {
                fout << "1 ";
            } else {
                fout << "0 ";
            }
        }
        fout << std::endl;
    }
    */


    auto iter = list.begin();
    for (int i = (origSize-list.size()); i < origSize; i++) {
        int a = std::get<0>(std::get<0>(*iter));
        int b = std::get<1>(std::get<0>(*iter));
        double pos = std::get<1>(*iter);


    }
}

void smoothMesh(TetMesh & mesh);
void writeVTK(const std::string& filename, const TetMesh &mesh);
void writeOFF(const std::string& filename, const TetMesh &mesh);
void writeDolfin(const std::string &filename, const TetMesh &mesh);
void writeTriangle(const std::string &filename, const TetMesh &mesh);

// void writeMCSF(const std::string &filename, const TetMesh &mesh);
//void writeDiffPack
//void writeCARP
