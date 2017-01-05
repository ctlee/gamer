#pragma once

#include <iostream>
#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include "util.h"
#include "SimplicialComplex.h"
#include "SimplicialComplexVisitors.h"
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
auto getTangentH(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::SimplexID<SurfaceMesh::topLevel> curr)
{
    return (*curr).orientation;
}

template <std::size_t level, std::size_t dimension>
auto getTangentH(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::SimplexID<level> curr)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;
    auto cover = mesh.get_cover(curr);
    for(auto alpha : cover)
    {
        auto edge = *mesh.get_edge_up(curr, alpha);
        const auto& v = (*mesh.get_simplex_up({alpha})).position;
        auto next = mesh.get_simplex_up(curr,alpha);
        rval += edge.orientation * (v-origin) * getTangentH(mesh, origin, next);
    }

    return rval/cover.size();
}

template <std::size_t dimension>
auto getTangentF(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, 
        SurfaceMesh::SimplexID<SurfaceMesh::topLevel> curr, std::set<SurfaceMesh::KeyType>& cover)
{
    return (*curr).orientation;
}

template <std::size_t level, std::size_t dimension>
auto getTangentF(const SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, 
        SurfaceMesh::SimplexID<level> curr, std::set<SurfaceMesh::KeyType>& cover)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;
    for(auto alpha : cover)
    {
        auto edge = *mesh.get_edge_up(curr, alpha);
        const auto& v = (*mesh.get_simplex_up({alpha})).position;
        auto next = mesh.get_simplex_up(curr,alpha); 
        auto coverup = cover;
        coverup.erase(alpha);
        rval += edge.orientation * (v-origin) * getTangentF(mesh, origin, next, coverup);
    }
    return rval/cover.size();
}

/**
 * READERS AND WRITERS
 */
// Geomview OFF
std::unique_ptr<SurfaceMesh> readOFF(const std::string& filename);
void writeOFF(const std::string& filename, const SurfaceMesh& mesh);

// Wavefront OBJ
std::unique_ptr<SurfaceMesh> readOBJ(const std::string& filename);
void writeOBJ(const std::string& filename, const SurfaceMesh& mesh);

void print(const SurfaceMesh& mesh);
void generateHistogram(const SurfaceMesh& mesh);
int getValence(const SurfaceMesh& mesh, const SurfaceMesh::SimplexID<1> vertexID);
double getArea(const SurfaceMesh& mesh);
double getArea(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<3> faceID);
double getVolume(const SurfaceMesh& mesh);

void edgeFlip(SurfaceMesh& mesh, SurfaceMesh::SimplexID<2> edgeID);
std::vector<SurfaceMesh::SimplexID<2>> selectFlipEdges(const SurfaceMesh& mesh, bool preserveRidges, 
        std::function<bool(const SurfaceMesh&, SurfaceMesh::SimplexID<2>&)> &checkFlip);
bool checkFlipAngle(const SurfaceMesh& mesh, const SurfaceMesh::SimplexID<2>& edgeID);
bool checkFlipValence(const SurfaceMesh& mesh, const SurfaceMesh::SimplexID<2>& edgeID);

void barycenterVertexSmooth(SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID);
void normalSmooth(SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID);

tensor<double,3,2> getTangent(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID);
tensor<double,3,2> getTangent(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<3> faceID);
Vector getNormalFromTangent(const tensor<double,3,2> tangent);
Vector getNormal(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID);
Vector getNormal(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<3> faceID);

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
    bool visit(const SurfaceMesh& F, SurfaceMesh::SimplexID<level> s)
    {
        // auto tan = getTangent(F, s);
        // auto norm = getNormalFromTangent(tan);
        auto norm = getNormal(F,s);
        lst += norm*norm;
        return true;
    }
};

Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double,3,2> mat);

template <std::size_t rings>
void weightedVertexSmooth(SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID){
    auto centerName = mesh.get_name(vertexID)[0]; 
    auto& center = *vertexID; // get the vertex data

    double sumWeights = 0;
    Vector newPos;

    /**
     * Compute the following sum to get the new position
     * \bar{x} = \frac{1}{\sum_{i=1}^{N_2}(\alpha_i+1)}\sum_{i=1}^{N_2}(alpha_i + 1) x_i
     */
    for(auto edge : mesh.up(vertexID)){
        // Get the vertex connected by edge
        auto edgeName = mesh.get_name(edge);
        auto shared = *mesh.get_simplex_up({(edgeName[0] == centerName) ? edgeName[1] : edgeName[0]});

        // Get the vertices connected to adjacent edge
        auto up = mesh.get_cover(edge);
        auto prev = *mesh.get_simplex_up({up[0]});
        auto next = *mesh.get_simplex_up({up[1]}); 
        
        auto pS = prev - shared;
        pS /= std::sqrt(pS|pS);
        auto nS = next - shared;
        nS /= std::sqrt(nS|nS);
        auto bisector = (pS + nS)/2;
        bisector /= std::sqrt(bisector|bisector);

        // Get a reference vecter to shared which lies on the plane of interest.
        auto disp = center - shared;
        Eigen::Map<Eigen::Vector3d> disp_e(disp.data());

        //auto tanNorm = getNormalFromTangent(pS^nS);
        auto tanNorm = cross(pS, nS);

        // Get the perpendicular plane made up of plane normal of bisector
        //auto perpPlane = tanNorm^bisector;
        //auto perpNorm = getNormalFromTangent(perpPlane);
        auto perpNorm = cross(tanNorm, bisector);
        perpNorm /= std::sqrt(perpNorm|perpNorm);
        auto perpProj = perpNorm*perpNorm; // tensor product
      
        // Compute perpendicular component 
        Vector perp;
        Eigen::Map<Eigen::Matrix3d> perpProj_e(perpProj.data());
        Eigen::Map<Eigen::Vector3d> perp_e(perp.data());
        perp_e = perpProj_e*disp_e;

        auto alpha = (pS|nS)+1; // keep the dot product positive
        sumWeights += alpha;
        newPos += alpha*(center.position - perp);
    }
    newPos /= sumWeights;
    /**
     * Scale by PCA
     * \bar{x} = x + \sum_{k=1}^3 \frac{1}{1+\lambda_k}((\bar{x} - x)\cdot \vec{e_k})\vec{e_k}
     */
    auto v = LocalStructureTensorVisitor();
    visit_neighbors_up<rings>(v, mesh, vertexID); // TODO: how should we set this?
    auto eigen_result = getEigenvalues(v.lst);
    // std::cout << "The eigenvalues of A are:\n" << eigen_result.eigenvalues() << std::endl;
    // std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
    //      << "corresponding to these eigenvalues:\n"
    //      << eigen_result.eigenvectors() << std::endl;

    newPos -= center.position;
    Eigen::Map<Eigen::Vector3d> newPos_e(newPos.data());
    
    auto w = ((eigen_result.eigenvectors().transpose()*newPos_e).array() // dot product
             / (eigen_result.eigenvalues().array()+1)).matrix();         // elementwise-division
    newPos_e = eigen_result.eigenvectors()*w; // matrix product
    center.position += newPos;
}