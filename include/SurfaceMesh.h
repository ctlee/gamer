#pragma once

#include <iostream>
#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>
#include <string>
#include <unordered_set>
#include <utility>
#include "util.h"
#include "SimplicialComplex.h"
#include "Orientable.h"
#include "Vertex.h"
/**
 * @brief      Properties that Faces should have
 */
struct FaceProperties
{
    int  marker;   /**< @brief Marker */
    bool selected; /**< @brief selection flag */
};

/**
 * @brief      Face object
 */
struct Face : Orientable, FaceProperties
{
    Face() {}
    Face(Orientable orient, FaceProperties prop)
        : Orientable(orient), FaceProperties(prop)
    {}
};

/**
 * @brief      Type for containing root metadata
 */
struct Global
{
    bool  closed;                /**< @brief is the surface mesh closed or not */
    int   _marker;               /**< @brief doman marker, to be used when tetrahedralizing */
    float volume_constraint;     /**< @brief volume constraint of the tetrahedralized domain */
    bool  use_volume_constraint; /**< @brief flag that determines if the volume constraint is used */
    float min[3];                /**< @brief minimal coordinate of nodes */
    float max[3];                /**< @brief maximal coordinate of nodes */
    float avglen;                /**< @brief average edge length */
    bool hole;                   /**< @brief flag that determines if the mesh is a hole or not */
};

/**
 * @brief      A helper struct containing the traits/types in the simplicial
 *             complex
 */
struct complex_traits
{
    using KeyType = int;                                                    /**< @brief the index type */
    using NodeTypes = util::type_holder<Global,Vertex,void,Face>;           /**< @brief the types of each Node */
    using EdgeTypes = util::type_holder<Orientable,Orientable,Orientable>;  /**< @brief the types of each Edge */
};

// This alias is for legacy purposes...
using SurfaceMesh_ASC = simplicial_complex<complex_traits>;
using ASC = simplicial_complex<complex_traits>; // Alias for the lazy
using SurfaceMesh = simplicial_complex<complex_traits>;

template <std::size_t dimension>
auto getTangentH(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::NodeID<SurfaceMesh::topLevel> curr)
{
    return (*curr).orientation;
}

template <std::size_t level, std::size_t dimension>
auto getTangentH(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::NodeID<level> curr)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;
    auto cover = mesh.get_cover(curr);
    for(auto alpha : cover)
    {
        auto edge = *mesh.get_edge_up(curr, alpha);
        const auto& v = (*mesh.get_node_up({alpha})).position;
        auto next = mesh.get_node_up(curr,alpha);
        rval += edge.orientation * (v-origin) * getTangentH(mesh, origin, next);
    }

    return rval/cover.size();
}

template <std::size_t dimension>
auto getTangentF(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, 
        SurfaceMesh::NodeID<SurfaceMesh::topLevel> curr, std::set<SurfaceMesh::KeyType>& cover)
{
    return (*curr).orientation;
}

template <std::size_t level, std::size_t dimension>
auto getTangentF(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, 
        SurfaceMesh::NodeID<level> curr, std::set<SurfaceMesh::KeyType>& cover)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;
    for(auto alpha : cover)
    {
        auto edge = *mesh.get_edge_up(curr, alpha);
        const auto& v = (*mesh.get_node_up({alpha})).position;
        auto next = mesh.get_node_up(curr,alpha); 
        auto coverup = cover;
        coverup.erase(alpha);
        rval += edge.orientation * (v-origin) * getTangentF(mesh, origin, next, coverup);
    }
    return rval/cover.size();
}

/**
 * @brief      Reads off.
 *
 * @param[in]  filename  The filename
 *
 * @return     { description_of_the_return_value }
 */
std::pair<SurfaceMesh*, bool> readOFF(const std::string& filename);

/**
 * @brief      Writes off.
 *
 * @param[in]  filename  The filename
 * @param[in]  mesh      The mesh
 */
void writeOFF(const std::string& filename, SurfaceMesh& mesh);

void print(const SurfaceMesh& mesh);

void generateHistogram(const SurfaceMesh& mesh);
bool smoothMesh(const SurfaceMesh& mesh, std::size_t minAngle, std::size_t maxAngle, std::size_t maxIter, bool preserveRidges);
void edgeFlip(SurfaceMesh& mesh, SurfaceMesh::NodeID<2> edgeID);
std::vector<SurfaceMesh::NodeID<2>> selectFlipEdges(SurfaceMesh& mesh, bool preserveRidges, 
        const std::function<bool(SurfaceMesh&, SurfaceMesh::NodeID<2>&)> &checkFlip);
bool checkFlipAngle(const SurfaceMesh& mesh, const SurfaceMesh::NodeID<2>& edgeID);
bool checkFlipValence(const SurfaceMesh& mesh, const SurfaceMesh::NodeID<2>& edgeID);

void barycenterVertexSmooth(SurfaceMesh& mesh, SurfaceMesh::NodeID<1> vertexID);
void weightedVertexSmooth(SurfaceMesh& mesh, SurfaceMesh::NodeID<1> vertexID);

int getValence(const SurfaceMesh& mesh, const SurfaceMesh::NodeID<1> vertexID);

tensor<double,3,2> getTangent(const SurfaceMesh& mesh, SurfaceMesh::NodeID<1> vertexID);
tensor<double,3,2> getTangent(const SurfaceMesh& mesh, SurfaceMesh::NodeID<3> faceID);
Vector getNormalFromTangent(const tensor<double,3,2> tangent);

// These exist for the a potential python interface
void translate(SurfaceMesh& mesh, Vector v);
void translate(SurfaceMesh& mesh, double dx, double dy, double dz);
void scale(SurfaceMesh& mesh, Vector v);
void scale(SurfaceMesh& mesh, double sx, double sy , double sz);
void scale(SurfaceMesh& mesh, double s);

struct LocalStructureTensorVisitor
{
    tensor<double,3,2> lst;

    LocalStructureTensorVisitor(){
        lst = tensor<double,3,2>();
    }

    template <std::size_t level> 
    bool visit(const SurfaceMesh& F, SurfaceMesh::NodeID<level> s)
    {
        auto tan = getTangent(F, s);
        auto norm = getNormalFromTangent(tan);
        lst += norm*norm;
        //std::cout << norm*norm << " " << lst << std::endl; 
        return true;
    }
};

Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double,3,2> mat);