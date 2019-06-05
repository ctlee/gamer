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

    bool operator< (const Edge& rhs) const
    {
        return  error < rhs.error;
    };

    bool operator> (const Edge& rhs) const
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


    const bool operator< (const TetVertex& rhs)
    {
        return  error < rhs.error;
    };

    const bool operator> (const TetVertex& rhs)
    {
        return  error > rhs.error;
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

typedef double (*PenaltyFunction)(TetMesh::SimplexID<2>, TetMesh&);

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


// Misc Error computation and propagation functions
template <std::size_t level>
void propagateHelper(TetMesh & mesh){
    if (level <= 1) return;

    for (auto node : mesh.get_level_id<level>()){
        auto parents = mesh.up(node);
        double tmperror;
        for (auto parent : parents){
            tmperror += parent->error;
        }
        double error = tmperror / parents.size();
        node->error = error;
    }
    propagateHelper<level-1>(mesh);
}





/*
 * Compute errors of each vertex/edge based on inputted functions
 * vertex Loc is location of new vertex relative to edge endpoints
 */
double computePenalty(TetMesh::SimplexID<2> & e, TetMesh& mesh, double vertexLoc,
                      std::vector<PenaltyFunction> penaltyList){
    // Scalarization of penalty function; using weighted sum method with 1/n as
    // coefficient for each individual function
    if (penaltyList.size() == 0) {
        return 0;
    }
    double weight = (1/penaltyList.size());

    double penalty = 0;
    for (auto f : penaltyList) {
        penalty += weight * f(e, mesh);
        if (std::isnan(penalty)) {
            return INFINITY;
        }
    }
    return penalty;
}


void propagateError(TetMesh & mesh, std::vector<PenaltyFunction> penalty){
    for (auto edge : mesh.get_level_id<2>()) {
        (*edge).error = computePenalty(edge, mesh, 0, penalty);
    }
    //propagateHelper<2>(mesh);
}

// Specific decimation operations
template <typename Complex, template <typename> class Callback>
int halfEdgeCollapse(TetMesh & mesh, TetMesh::SimplexID<2> e, std::vector<PenaltyFunction> penaltyList,
        Callback<Complex> &&clbk) {
    double penaltyA = computePenalty(e, mesh, 0, penaltyList);
    double penaltyB = computePenalty(e, mesh, 1, penaltyList);
    if (penaltyA < penaltyB) {
        edgeCollapse(mesh, e, 0, Callback<TetMesh>());
        return 0;
    } else {
        edgeCollapse(mesh, e, 1, Callback<TetMesh>());
        return 1;
    }
}


// Penalty/Constraint functions
double isBoundaryEdge(TetMesh::SimplexID<2> edge, TetMesh & mesh) {
    std::set<TetMesh::SimplexID<3>> faces= mesh.up(edge);
    for (auto face : faces) {
        std::set<TetMesh::SimplexID<4>> parentTets = mesh.up(face);
        if (parentTets.size() == 1) {
            return INFINITY;
        }
    }
    return 0;
}

bool isBoundaryVertex(TetMesh::SimplexID<1> v, TetMesh & mesh) {
    std::set<TetMesh::SimplexID<2>> edges = mesh.up(v);
    for (auto edge : edges) {
        if (isBoundaryEdge(edge, mesh)) {
            return true;
        }
    }
    return false;
}

double edgeLength(TetMesh::SimplexID<2> edge, TetMesh & mesh) {
    auto name =  mesh.get_name(edge);
    auto v = *mesh.get_simplex_down(edge, name[0])
             - *mesh.get_simplex_down(edge, name[1]);
    return std::sqrt(v|v);
}

double scalarInfoValue(tetmesh::TetEdge e) {}

double volumePreservation(TetMesh::SimplexID<2> edge, TetMesh & mesh) {}

double inwardCollapse(TetMesh::SimplexID<2> edge, TetMesh & mesh) {
    auto name =  mesh.get_name(edge);
    auto v1 = mesh.get_simplex_down(edge, name[0]);
    auto v2 = mesh.get_simplex_down(edge, name[1]);
    auto b1 = isBoundaryVertex(v1, mesh);
    auto b2 = isBoundaryVertex(v2, mesh);
    if (b1 && b2) {
        return INFINITY;
    } else if (!b1 && !b2) {
        return INFINITY;
    } else {
        return 0;
    }
}

// Test if a potential edge collapse preserves substructes
double substructurePreservation(TetMesh::SimplexID<2> edge, TetMesh & mesh) {}

// Test if a potential edge collapse results in inverted tetrahedra
double simplexInversionCheck(TetMesh::SimplexID<2> edge, TetMesh & mesh) {}

// Error from error quadrics as described in garland paper
double quadricError(TetMesh::SimplexID<2> edge, TetMesh & mesh) {}

void markDecimatedEdge(TetMesh & mesh, std::vector<std::pair<std::pair<int,int>, double>>& list, TetMesh::SimplexID<2> edge,
                       double loc) {
    auto name =  mesh.get_name(edge);
    list.push_back(std::pair<std::pair<int,int>,double>(std::pair<int,int>(name[0], name[1]), loc));
}

// Restriction matrix with convention as in http://www.cs.huji.ac.il/~csip/CSIP2007-MG.pdf
void writeRestrictionMatrix(const std::string &filename, int origSize,
        std::vector<std::pair<std::pair<int,int>, double>> list, std::vector<bool> encountered) {

    std::ofstream fout(filename);
    if(!fout.is_open())
    {
        std::cerr   << "File '" << filename
                    << "' could not be writen to." << std::endl;
        exit(1);
    }

    std::vector<std::vector<double>> matrix(origSize, std::vector<double> (origSize, 0));
    auto iter = list.begin();

    for (int i = 0; i < (origSize-list.size()); i++) {
        if (!encountered[i]) {
            matrix[i][i] = 1;
        }
    }

    for (int i = (origSize-list.size()); i < matrix.size(); i++) {
        int a = std::get<0>(std::get<0>(*iter));
        int b = std::get<1>(std::get<0>(*iter));
        double pos = std::get<1>(*iter);
        matrix[i][a] = pos;
        matrix[i][b] = 1 - pos;
        iter++;
    }

    // Delete rows that are all 0
    for (auto itr = matrix.begin(); itr != matrix.end();) {
        bool allzeros = true;
        for (auto col : *itr) {
            if (col != 0) {
                allzeros = false;
                break;
            }
        }
        if (allzeros) {
            itr = matrix.erase(itr);
        } else {
            ++itr;
        }

    }

    // Print to file as csv
    for (auto row : matrix) {
        for (auto element : row) {
            fout << element << ",";
        }
        fout << std::endl;
    }
}

template <typename Complex, template <typename> class Callback>
void decimation(TetMesh & mesh, double threshold, Callback<Complex> &&cbk){
    std::vector<PenaltyFunction> samplePenaltyFunction = {isBoundaryEdge, inwardCollapse};

    int origVertexNum = mesh.size<1>();
    int numToRemove = origVertexNum - origVertexNum*threshold;

    if (numToRemove < 1) {
        throw std::invalid_argument("Number of vertices to remove to must be at least 1");
    }

    // Exit, mesh already sufficiently decimated
    if (threshold*origVertexNum > mesh.size<2>()) {
        std::cout << "Decimation finished to requested threshold: " << threshold <<"; final mesh size: "
        << mesh.size<2>() << std::endl;
        return;
    }

    propagateError(mesh, samplePenaltyFunction);

    std::vector<TetMesh::SimplexID<2>> sortedEdges;
    // initialize vector of edges
    for (TetMesh::SimplexID<2> e : mesh.get_level_id<2>()) {
        sortedEdges.push_back(e);
    }

    auto edgeErrorCmp = [](const TetMesh::SimplexID<2> e1, const TetMesh::SimplexID<2> e2) -> bool{
        return *e1 < *e2;
    };

    // sort edges in order of error
    std::sort_heap(sortedEdges.begin(), sortedEdges.end(), edgeErrorCmp);
    // remove invalid edges (if edge with both vertices decimated already decimated, will cause segfault)
    std::vector<bool> encountered(mesh.size<1>(), false);
    int numValidEdges = 0;
    for (auto it = sortedEdges.begin(); it != sortedEdges.end();) {
        if (numValidEdges > numToRemove) {
            sortedEdges.erase(it, sortedEdges.end());
            break;
        }

        auto edge = *it;

        /*
        if ((*edge).error == INFINITY) {
            sortedEdges.erase(it, sortedEdges.end());
            std::cout << "No more valid edges for decimation, all remaining edges are protected by constraints." << std::endl;
            break;
        }
        */
        auto name =  mesh.get_name(edge);

        if (encountered[name[0]] || encountered[name[1]]) {
            it = sortedEdges.erase(it);
        } else {
            ++it;
            numValidEdges++;
        }
        encountered[name[0]] = true;
        encountered[name[1]] = true;
    }

    std::vector<std::pair<std::pair<int,int>, double>> decimatedList;

    for(auto it = sortedEdges.begin(); it != sortedEdges.end(); it++){
        auto edge = *it;
        std::cout << "collapsing edge: " << *it << std::endl;
        //markDecimatedEdge(mesh, decimatedList, edge, .5);
        if (isBoundaryEdge(*it, mesh)) {
            std::cout << "Warning: boundary edges being collapsed, raise threshold" << std::endl;
            break;
        }
        edgeCollapse(mesh, *it, .5, Callback<TetMesh>());
        //int pos = halfEdgeCollapse(mesh, *it, samplePenaltyFunction, Callback<TetMesh>());
    }
    std::cout << "Decimation finished, final size: " << mesh.size<1>() << " vertices." << std::endl;
    //writeRestrictionMatrix("restrictionM.csv", origVertexNum, decimatedList, encountered);
}


void smoothMesh(TetMesh & mesh);
void writeVTK(const std::string& filename, const TetMesh &mesh);
void writeOFF(const std::string& filename, const TetMesh &mesh);
std::unique_ptr<TetMesh> readDolfin(const std::string&filename);
void writeDolfin(const std::string &filename, const TetMesh &mesh);
void writeTriangle(const std::string &filename, const TetMesh &mesh);

// void writeMCSF(const std::string &filename, const TetMesh &mesh);
//void writeDiffPack
//void writeCARP
