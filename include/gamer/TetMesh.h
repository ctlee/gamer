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

#include <gamer/gamer.h>

#define TETLIBRARY
#include <tetgen.h>

#include <casc/casc>

#include "gamer/Vertex.h"
#include "gamer/SurfaceMesh.h"

namespace tetmesh
{

struct ErrorProperties
{
    double error;
};

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

    bool operator< (const TetVertex& rhs) const
    {
        return  error < rhs.error;
    };

    bool operator> (const TetVertex& rhs) const
    {
        return  error > rhs.error;
    };
};


struct Edge : Vertex, ErrorProperties
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
struct Face : FaceProperties, ErrorProperties
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

struct Cell : casc::Orientable, CellProperties, ErrorProperties
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

typedef double (*PenaltyFunction)(TetMesh::SimplexID<2>, TetMesh&);

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
                      std::vector<PenaltyFunction> penaltyList);

void propagateError(TetMesh & mesh, std::vector<PenaltyFunction> penalty);

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
double isBoundaryEdge(TetMesh::SimplexID<2> edge, TetMesh & mesh);

bool isBoundaryVertex(TetMesh::SimplexID<1> v, TetMesh & mesh);

double edgeLength(TetMesh::SimplexID<2> edge, TetMesh & mesh);

// double scalarInfoValue(tetmesh::TetEdge e) {}

// double volumePreservation(TetMesh::SimplexID<2> edge, TetMesh & mesh) {}

double inwardCollapse(TetMesh::SimplexID<2> edge, TetMesh & mesh);

// // Test if a potential edge collapse preserves substructes
// double substructurePreservation(TetMesh::SimplexID<2> edge, TetMesh & mesh);

// // Test if a potential edge collapse results in inverted tetrahedra
// double simplexInversionCheck(TetMesh::SimplexID<2> edge, TetMesh & mesh);

// // Error from error quadrics as described in garland paper
// double quadricError(TetMesh::SimplexID<2> edge, TetMesh & mesh);

void markDecimatedEdge(TetMesh & mesh, std::vector<std::pair<std::pair<int,int>, double>>& list, TetMesh::SimplexID<2> edge,
                       double loc);

// Restriction matrix with convention as in http://www.cs.huji.ac.il/~csip/CSIP2007-MG.pdf
void writeRestrictionMatrix(const std::string &filename, int origSize,
        std::vector<std::pair<std::pair<int,int>, double>> list, std::vector<bool> encountered);

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
