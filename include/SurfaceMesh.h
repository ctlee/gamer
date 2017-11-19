/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2017
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
#include <libraries/casc/include/SimplicialComplex.h>
#include <libraries/casc/include/CASCTraversals.h>
#include <libraries/casc/include/util.h>
#include <libraries/casc/include/Orientable.h>
#include <libraries/casc/include/stringutil.h>

#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include "Vertex.h"

// The number of rings to use to compute local structure tensor
#define RINGS 2

/**
 * @brief      Properties that Faces should have
 */
struct FaceProperties
{
    int  marker;   /**< @brief Marker */
    bool selected; /**< @brief Selection flag */
};

/**
 * @brief      Face object
 */
struct Face : casc::Orientable, FaceProperties
{
    /// Default constructor
    Face() {}
    /**
     * @brief      Constructor
     *
     * @param[in]  orient  Orientable object
     * @param[in]  prop    Properties of a face
     */
    Face(Orientable orient, FaceProperties prop)
        : Orientable(orient), FaceProperties(prop)
    {}
};

/**
 * @brief      Type for containing root metadata
 */
struct Global
{
    /// Is the SurfaceMesh closed or not.
    bool  closed;
    /// Domain marker to be used when tetrahedralizing.
    int   _marker;
    /// Volume constraint of the tetrahedralized domain.
    float volume_constraint;
    /// flag that determines if the volume constraint is used.
    bool  use_volume_constraint;
    /// Minimum coordinate of vertices
    float min[3];
    /// Max coordinate of vertices
    float max[3];
    /// Average edge length
    float avglen;
    /// Flag that determines if the mesh has a hole or not
    bool  hole;
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
    using NodeTypes = util::type_holder<Global, Vertex, void, Face>;
    /// The types of each edge
    using EdgeTypes = util::type_holder<casc::Orientable, casc::Orientable, casc::Orientable>;
};

// This alias is for legacy purposes...
using SurfaceMesh_ASC = casc::simplicial_complex<complex_traits>;
using ASC = casc::simplicial_complex<complex_traits>; // Alias for the lazy
using SurfaceMesh = casc::simplicial_complex<complex_traits>;

/**
 * @brief      Reads in a GeomView OFF file.
 *
 * @param[in]  filename  The filename.
 *
 * @return     Returns a unique_ptr to the SurfaceMesh.
 */
std::unique_ptr<SurfaceMesh> readOFF(const std::string &filename);

/**
 * @brief      Write the SurfaceMesh to file in OFF format.
 *
 * @param[in]  filename  The filename to write to.
 * @param[in]  mesh      SurfaceMesh of interest.
 */
void writeOFF(const std::string &filename, const SurfaceMesh &mesh);

// Wavefront OBJ
std::unique_ptr<SurfaceMesh> readOBJ(const std::string &filename);
void writeOBJ(const std::string &filename, const SurfaceMesh &mesh);

void print(const SurfaceMesh &mesh);
void generateHistogram(const SurfaceMesh &mesh);
std::tuple<double, double, int, int> getMinMaxAngles(const SurfaceMesh& mesh, 
    int maxMinAngle, int minMaxAngle);
double getArea(const SurfaceMesh &mesh);
double getArea(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID);
double getVolume(const SurfaceMesh &mesh);
int getValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID);


/**
 * @brief      Terminal case
 *
 * @param[in]  mesh       The mesh
 * @param[in]  origin     The origin
 * @param[in]  curr       The curr
 *
 * @tparam     dimension  { description }
 *
 * @return     The tangent h.
 */
template <std::size_t dimension>
auto getTangentH(const SurfaceMesh &mesh,
                 const tensor<double, dimension, 1> &origin,
                 SurfaceMesh::SimplexID<SurfaceMesh::topLevel> curr)
{
    return (*curr).orientation;
}

/**
 * @brief      Vertex tangent by computing the average wedge product of incident
 *             faces.
 *
 * @param[in]  mesh       The mesh
 * @param[in]  origin     The origin
 * @param[in]  curr       The curr
 *
 * @tparam     level      { description }
 * @tparam     dimension  { description }
 *
 * @return     The tangent h.
 */
template <std::size_t level, std::size_t dimension>
auto getTangentH(const SurfaceMesh &mesh,
                 const tensor<double, dimension, 1> &origin,
                 SurfaceMesh::SimplexID<level> curr)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;
    auto cover = mesh.get_cover(curr);
    for (auto alpha : cover)
    {
        auto        edge = *mesh.get_edge_up(curr, alpha);
        const auto &v = (*mesh.get_simplex_up({alpha})).position;
        auto        next = mesh.get_simplex_up(curr, alpha);
        rval += edge.orientation * (v-origin) * getTangentH(mesh, origin, next);
    }
    return rval/cover.size();
}

/**
 * @brief      Compute the tangent of a face as the sum of wedge products.
 *
 * This is the terminal case. The wedge products must be scaled by the
 *
 * @param[in]  mesh       SurfaceMesh of interest.
 * @param[in]  origin     The Vector position of the first vertex.
 * @param[in]  curr       Current simplex to compute on.
 * @param      cover      Set of coboundary simplices.
 *
 * @tparam     level      Simplex dimension of curr.
 * @tparam     dimension  Dimension of the embedding.
 *
 * @return     Returns a 2-tensor corresponding to the tangent plane.
 */
template <std::size_t dimension>
auto getTangentF(const SurfaceMesh &mesh,
                 const tensor<double, dimension, 1> &origin,
                 SurfaceMesh::SimplexID<SurfaceMesh::topLevel> curr,
                 std::set<SurfaceMesh::KeyType> &cover)
{
    return (*curr).orientation;
}

/**
 * @brief      Compute the tangent of a face as the sum of wedge products.
 *
 * @param[in]  mesh       SurfaceMesh of interest.
 * @param[in]  origin     The Vector position of the first vertex.
 * @param[in]  curr       Current simplex to compute on.
 * @param      cover      Set of coboundary simplices.
 *
 * @tparam     level      Simplex dimension of curr.
 * @tparam     dimension  Dimension of the embedding.
 *
 * @return     Returns a 2-tensor corresponding to the tangent plane.
 */
template <std::size_t level, std::size_t dimension>
auto getTangentF(const SurfaceMesh &mesh,
                 const tensor<double, dimension, 1> &origin,
                 SurfaceMesh::SimplexID<level> curr,
                 std::set<SurfaceMesh::KeyType> &cover)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;
    for (auto alpha : cover)
    {
        auto        edge = *mesh.get_edge_up(curr, alpha);
        const auto &v = (*mesh.get_simplex_up({alpha})).position;
        auto        next = mesh.get_simplex_up(curr, alpha);
        auto        coverup = cover;
        coverup.erase(alpha);
        rval += edge.orientation * (v-origin) * getTangentF(mesh, origin, next, coverup);
    }
    return rval/cover.size();
}

/**
 * @brief      Comupute the tangent to a vertex.
 *
 * The vertex normal is defined as the mean of incident triangle normals as
 * follows,
 * \f$ \mathbf{n} = \frac{1}{N} \sum_{i=0}^{N, \textrm{incident triangle}}
 * n_t\f$.
 * where \f$\mathbf{n}\f$ is the vertex normal, \f$N\f$ is the number of
 * incident
 * faces, and \f$n_i\f$ is the normal of incident triangle \f$i\f$.
 *
 * @param[in]  mesh      SurfaceMesh of interest.
 * @param[in]  vertexID  SimplexID of the vertex to get the tangent of.
 *
 * @return     Returns a 2-tensor representing the tangent plane.
 */
tensor<double, 3, 2> getTangent(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID);

/**
 * @brief      Compute the tangent of a face.
 *
 * This function gets the tangent by computing the sum of oriented wedge
 * products.
 *
 * @param[in]  mesh    SurfaceMesh of interest.
 * @param[in]  faceID  SimplexID of the face to get the tangent of.
 *
 * @return     Returns a 2-tensor representing the tangent plane.
 */
tensor<double, 3, 2> getTangent(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID);

/**
 * @brief      Gets the normal vector from the tangent.
 *
 * @param[in]  tangent  2-tensor representing the tangent plane.
 *
 * @return     The normal Vector.
 */
Vector getNormalFromTangent(const tensor<double, 3, 2> tangent);

/**
 * @brief      Gets the normal of a vertex.
 *
 * The vertex normal is defined as the mean of incident triangle normals as
 * follows,
 * \f$\mathbf{n} = \frac{1}{N} \sum_{i=0}^{N, \textrm{incident triangle}}n_t\f$.
 * where \f$\mathbf{n}\f$ is the vertex normal, \f$N\f$ is the number of 
 * incident faces, and \f$n_i\f$ is the normal of incident triangle \f$i\f$.
 *
 * @param[in]  mesh      SurfaceMesh of interest.
 * @param[in]  vertexID  SimplexID of the vertex to get the normal of.
 *
 * @return     Returns the Vector normal to the vertex.
 */
Vector getNormal(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID);

/**
 * @brief      Compute the normal of a face.
 *
 * This function computes the normal as the cross product with respect to the
 * orientation. See also getTangent() and getNormalFromTangent().
 *
 * @param[in]  mesh    SurfaceMesh of interest.
 * @param[in]  faceID  SimplexID of the face to get the tangent of.
 *
 * @return     Returns a Vector normal to the face.
 */
Vector getNormal(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID);

// These exist for the a potential python interface
void translate(SurfaceMesh &mesh, Vector v);
void translate(SurfaceMesh &mesh, double dx, double dy, double dz);
void scale(SurfaceMesh &mesh, Vector v);
void scale(SurfaceMesh &mesh, double sx, double sy, double sz);
void scale(SurfaceMesh &mesh, double s);
std::pair<Vector, double> getCenterRadius(SurfaceMesh &mesh);
void centeralize(SurfaceMesh &mesh);

tensor<double,3,2> computeLocalStructureTensor(
        const SurfaceMesh &mesh, 
        const SurfaceMesh::SimplexID<1> vertexID, 
        const int rings);

bool smoothMesh(SurfaceMesh &mesh, int maxMinAngle, int minMaxAngle, int maxIter, bool preserveRidges);

void coarse(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight);

template <typename InsertIter>
void triangulateHole(SurfaceMesh &mesh, 
        std::vector<SurfaceMesh::SimplexID<1>> &boundary,
        const Face &fdata,
        InsertIter iter){
    // Terminal case
    if(boundary.size() == 3){
        // create the face
        auto a = mesh.get_name(boundary[0])[0];
        auto b = mesh.get_name(boundary[1])[0];
        auto c = mesh.get_name(boundary[2])[0];
        mesh.insert({a,b,c}, fdata);
        return;
    }

    // Construct a sorted vector of pairs... (valence, vertexID)
    std::vector<std::pair<int, SurfaceMesh::SimplexID<1>>> list;
    for(auto vertexID : boundary){
        list.push_back(std::make_pair(getValence(mesh, vertexID), vertexID));
    }
    std::sort(list.begin(), list.end(), [](
                const std::pair<int, SurfaceMesh::SimplexID<1>> &lhs, 
                const std::pair<int, SurfaceMesh::SimplexID<1>> &rhs){
            return lhs.first < rhs.first;
        });

    SurfaceMesh::SimplexID<1> v1, v2;
    v1 = list[0].second;

    // Find v1 and rotate so that it is first for easy splitting later
    auto v1it =  std::find(boundary.begin(), boundary.end(), v1);
    std::rotate(boundary.begin(), v1it, boundary.end());

    // Get the next lowest valence vertex
    for(auto it = ++list.begin(); it != list.end(); ++it){
        v2 = (*it).second;
        // Check that it is not already connected to v1
        if(v2 != boundary[1] && v2 != boundary.back()){
            break;
        }
    }
    // Insert new edge
    auto a = mesh.get_name(v1)[0];
    auto b = mesh.get_name(v2)[0];
    mesh.insert({a,b});
    // TODO: (0) Get the simplex from insert.
    *iter++ = mesh.get_simplex_up({a,b});

    auto v2it =  std::find(boundary.begin(), boundary.end(), v2);
    std::vector<SurfaceMesh::SimplexID<1>> other;
    other.push_back(v2);
    std::move(v2it+1, boundary.end(), std::back_inserter(other));
    boundary.erase(v2it+1, boundary.end());
    
    other.push_back(v1);
   
    // Recurse to fill sub-holes 
    triangulateHole(mesh, boundary, fdata, iter);
    triangulateHole(mesh, other, fdata, iter);
}

template <class K>
struct orientHoleHelper {};

template <std::size_t k>
struct orientHoleHelper<std::integral_constant<std::size_t, k>>{
    template <typename Iterator>
    static void apply(SurfaceMesh & mesh, 
            const std::set<int> &&names, 
            Iterator begin,
            Iterator end){
        std::vector<SurfaceMesh::SimplexID<k+1>> next;
        for(auto curr = begin; curr != end; ++curr){
            for(auto a : mesh.get_cover(*curr)){
                auto find = names.find(a);
                if(find != names.end()){
                    next.push_back(mesh.get_simplex_up(*curr, a));

                    int orient = 1;
                    for(auto b : mesh.get_name(*curr)){
                        if(a > b){
                            if(a > b)
                            {
                                orient *= -1;
                            }
                            else{
                                break;
                            }
                        }
                    }
                    (*mesh.get_edge_up(*curr,a)).orientation = orient;
                }
            }
        }
        orientHoleHelper<std::integral_constant<std::size_t,k+1>>::apply(mesh, std::move(names), next.begin(), next.end());
    } 
};

template <>
struct orientHoleHelper<std::integral_constant<std::size_t, SurfaceMesh::topLevel>>{
    template <typename Iterator>
    static void apply(SurfaceMesh &mesh, const std::set<int> &&names, Iterator begin, Iterator end){} 
};

template<typename Iterator>
bool computeHoleOrientation(SurfaceMesh &mesh, Iterator begin, Iterator end){
    bool orientable = true;
    for(auto currIT = begin; currIT != end; ++currIT){
        SurfaceMesh::SimplexID<2> curr = *currIT;
        auto w = mesh.get_cover(curr);

        if(w.size() == 1)
        {
            //std::cout << curr << ":" << w[0] << " ~ Boundary" << std::endl;
        }
        else if(w.size() == 2)
        {
            auto& edge0 = *mesh.get_edge_up(curr, w[0]);
            auto& edge1 = *mesh.get_edge_up(curr, w[1]);

            auto& node0 = *mesh.get_simplex_up(curr, w[0]);
            auto& node1 = *mesh.get_simplex_up(curr, w[1]);

            if(node0.orientation == 0)
            {
                if(node1.orientation == 0)
                {
                    node0.orientation = 1;
                    node1.orientation = -edge1.orientation * edge0.orientation * node0.orientation;
                }
                else
                {
                    node0.orientation = -edge0.orientation * edge1.orientation * node1.orientation;
                }
            }
            else
            {
                if(node1.orientation == 0)
                {
                    node1.orientation = -edge1.orientation * edge0.orientation * node0.orientation;
                }
                else
                {
                    if(edge0.orientation*node0.orientation + edge1.orientation*node1.orientation != 0)
                    {
                        orientable = false;
                        std::cout << "+++++" << std::endl;
                        std::cout << edge0.orientation << " : " << node0.orientation << std::endl;
                        std::cout << edge1.orientation << " : " << node1.orientation << std::endl;

                        std::cout << " : "
                                  << edge0.orientation*node0.orientation + edge1.orientation*node1.orientation
                                  << std::endl;
                        std::cout << "-----"
                                  << std::endl;
                        std::cout << "Non-Orientable: "
                                  << edge0.orientation*node0.orientation + edge1.orientation*node1.orientation
                                  << std::endl;
                    }
                }
            }
        }
    }
    return orientable;
}

/**
 * @brief      Compute the eigenvalues of a 3x3 matrix.
 *
 * @param[in]  mat   3x3 matrix to compute eigenvalues of.
 *
 * @return     Returns an Eigen::SelfAdjointEigenSolver containing the results
 */
Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double, 3, 2> mat);

/**
 * @brief      Smooth the vertex according to Section 2.2.2 of GAMer paper.
 *
 * @param[out] mesh      SurfaceMesh of interest.
 * @param[in]  vertexID  SimplexID of the vertex to move.
 * @param[in]  rings     The number of neighbor rings to use to compute the LST.
 */
void weightedVertexSmooth(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, int rings);

/**
 * @brief      Perform an edge flip operation
 *
 * @param      mesh    SurfaceMesh of interest.
 * @param[in]  edgeID  SimplexID of the edge to flip.
 */
void edgeFlip(SurfaceMesh &mesh, SurfaceMesh::SimplexID<2> edgeID);

/**
 * @brief      Select edges which are good candidates for flipping
 *
 * @param[in]  mesh            SurfaceMesh to operate on.
 * @param[in]  preserveRidges  Whether or not to try preserving ridges.
 * @param[in]  checkFlip       Functor specifying flip criteria
 * @param[in]  iter            Inserter for container of edges to flip
 *
 * @tparam     Inserter        Typename of the inserter.
 */
template <class Inserter>
void selectFlipEdges(const SurfaceMesh &mesh, 
                     bool preserveRidges,
                     std::function<bool(const SurfaceMesh &, const SurfaceMesh::SimplexID<2> &)> && checkFlip,
                     Inserter iter){
    casc::NodeSet<SurfaceMesh::SimplexID<2> > ignoredEdges;

    for (auto edgeID : mesh.get_level_id<2>())
    {
        if (!ignoredEdges.count(edgeID))
        {
            auto up = mesh.get_cover(edgeID);
            // The mesh is not a surface mesh...
            if (up.size() > 2)
            {
                // std::cerr << "This edge participates in more than 2 faces. "
                //           << "Returning..." << std::endl;
                continue;
            }
            else if (up.size() < 2) // Edge is a boundary
            {
                // std::cerr << "This edge participates in fewer than 2 faces. "
                //           << "Returning..." << std::endl;
                continue;
            }
            auto name = mesh.get_name(edgeID);
            std::pair<Vertex, Vertex> shared;
            shared.first  = *mesh.get_simplex_up({name[0]});
            shared.second = *mesh.get_simplex_up({name[1]});


            std::pair<Vertex, Vertex> notShared;
            notShared.first  = *mesh.get_simplex_up({up[0]});
            notShared.second = *mesh.get_simplex_up({up[1]});

            // Check if the edge is a part of a tetrahedron.
            if (mesh.exists<2>({up[0], up[1]}))
            {
                // std::cerr << "Found a tetrahedron cannot edge flip." 
                //           << std::endl;
                continue;
            }

            // Check if we're on a ridge. This prevents folding also.
            if (preserveRidges)
            {
                auto a   = getNormal(mesh, mesh.get_simplex_up(edgeID, up[0]));
                auto b   = getNormal(mesh, mesh.get_simplex_up(edgeID, up[1]));
                auto val = angle(a, b);
                if (val > 60)
                {
                    continue;
                }
            }

            // Check if the triangles make a wedge shape. Flipping this can
            // cause knife faces.
            auto f1   = (shared.first - notShared.first)^(shared.second - notShared.first);
            auto f2   = (shared.first - notShared.second)^(shared.second - notShared.second);
            auto f3   = (notShared.first - shared.first)^(notShared.second - shared.first);
            auto f4   = (notShared.first - shared.second)^(notShared.second - shared.second);
            auto area = std::pow(std::sqrt(f1|f1) + std::sqrt(f2|f2), 2);
            auto areaFlip = std::pow(std::sqrt(f3|f3) + std::sqrt(f4|f4), 2);
            if (areaFlip/area > 1.01) continue;

            // Check the flip using user function
            if (checkFlip(mesh, edgeID))
            {
                *iter++ = edgeID;   // Insert into edges to flip
                // The local topology will be changed. Don't flip neighbors 
                // which share a common face. 
                neighbors(mesh, edgeID, std::inserter(ignoredEdges, ignoredEdges.end()));
            }
        }
    }
}
/**
 * @brief      Check if we should flip an edge to improve the angles.
 *             
 * @param[in]  mesh    SurfaceMesh to manipulate.
 * @param[in]  edgeID  SimplexID of the edge to consider.
 *
 * @return     Returns true if flippiing the edge will improve angles, false 
 *             otherwise.
 */
bool checkFlipAngle(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID);

/**
 * @brief      Check if we should flip the edge to improve the valence.
 *
 * @param[in]  mesh    SurfaceMesh to manipulate.
 * @param[in]  edgeID  SimplexID of the edge to consider.
 *
 * @return     Returns true if flipping the edge will improve the valence, false
 *             otherwise.
 */
bool checkFlipValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID);

/**
 * @brief      Traditional barycenter smooth. 
 *
 * @param      mesh      SurfaceMesh to manipulate.
 * @param[in]  vertexID  SimplexID of the vertex to move.
 */
void barycenterVertexSmooth(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID);

void normalSmooth(SurfaceMesh &mesh);
void normalSmoothH(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID);

/**
 * @brief      Refine the mesh by quadrisection of faces
 *
 * @param      mesh  The mesh
 */
std::unique_ptr<SurfaceMesh> refineMesh(const SurfaceMesh &mesh);
