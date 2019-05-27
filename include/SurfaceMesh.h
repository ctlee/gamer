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
#include <stdexcept>
#include <string>
#include <memory>
#include <unordered_set>
#include <utility>

#include <casc/casc>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Vertex.h"

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
    Face() : Face(Orientable{0}, FaceProperties{0, false}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  marker    The marker
     * @param[in]  selected  The selected
     */
    Face(int marker, bool selected) : Face(Orientable{0}, FaceProperties{marker, selected}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  marker    The marker
     * @param[in]  selected  The selected
     */
    Face(int orient, int marker, bool selected) : Face(Orientable{orient}, FaceProperties{marker, selected}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  orient  Orientable object
     * @param[in]  prop    Properties of a face
     */
    Face(Orientable orient, FaceProperties prop)
        : Orientable(orient), FaceProperties(prop)
    {}

    friend std::ostream& operator<<(std::ostream& output, const Face& f){
        output  << "Face(o:" << f.orientation
                << ";m:" << f.marker
                << ";sel:" << std::boolalpha << f.selected << ")";
        return output;
    }

    std::string to_string() const{
        std::ostringstream output;
        output  << "Face(orient:" << orientation
                << ";m:" << marker
                << ";sel:" << std::boolalpha << selected << ")";
        return output.str();
    }

};

struct Edge{
    bool selected;
    // Default constructor
    Edge() : Edge(0) {}
    Edge(bool select) : selected(select) {}
};

/**
 * @brief      Type for containing root metadata
 */
struct Global
{
    /// Domain marker to be used when tetrahedralizing.
    int   marker;
    /// Volume constraint of the tetrahedralized domain.
    float volumeConstraint;
    /// flag that determines if the volume constraint is used.
    bool  useVolumeConstraint;
    /// Flag that determines if the mesh represents a hole or not
    bool  ishole;

    // Default constructor
    Global(int marker=-1, float volumeConstraint=-1, bool useVolumeConstraint=false, bool ishole=false) :
        marker(marker), volumeConstraint(volumeConstraint), useVolumeConstraint(useVolumeConstraint), ishole(ishole) {}
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
    using NodeTypes = util::type_holder<Global, Vertex, Edge, Face>;
    /// The types of each edge
    using EdgeTypes = util::type_holder<casc::Orientable, casc::Orientable, casc::Orientable>;
};

using SurfaceMesh = casc::simplicial_complex<complex_traits>;

/**
 * @brief      Compute the tangent to a vertex.
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



namespace surfacemesh_detail{
/**
 * @brief      Gets the mean edge length.
 *
 * @param      mesh  Surface mesh of interest
 *
 * @return     The mean edge length.
 */
double getMeanEdgeLength(SurfaceMesh& mesh);

/**
 * @brief      Remove a vertex form mesh and triangulate the resulting hole.
 *
 * @param      mesh      The mesh
 * @param[in]  vertexID  The vertex id
 */
void decimateVertex(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID);

tensor<double,3,2> computeLocalStructureTensor(
        const SurfaceMesh &mesh,
        const SurfaceMesh::SimplexID<1> vertexID,
        const int rings);

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
 * @brief      General template
 *
 * @tparam     K     { description }
 */
template <class K>
struct orientHoleHelper {};

template <std::size_t k>
struct orientHoleHelper<std::integral_constant<std::size_t, k>>
{
    /**
     * @brief      Set the orientation of a filled hole
     *
     * @param      mesh       Simplicial complex
     * @param[in]  names      Names of participant vertices
     * @param[in]  begin      Begin iterator of simplices to traverse
     * @param[in]  end        Past the end iterator of simplices to traverse
     *
     * @tparam     Iterator   Typename of the iterator
     */
    template <typename Iterator>
    static void apply(SurfaceMesh & mesh,
            const std::set<int> &&names,
            Iterator begin,
            Iterator end)
    {
        std::vector<SurfaceMesh::SimplexID<k+1>> next;
        for(auto curr = begin; curr != end; ++curr)
        {
            auto currSimplexID = *curr;
            for(auto a : mesh.get_cover(currSimplexID))
            {
                // Look for key a in names
                auto find = names.find(a);
                if(find != names.end())
                {
                    next.push_back(mesh.get_simplex_up(currSimplexID, a));

                    int orient = 1;
                    for(auto b : mesh.get_name(currSimplexID)){
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
                    //std::cout << casc::to_string(mesh.get_name(*curr)) << " " << a << " " << orient << std::endl;
                    (*mesh.get_edge_up(currSimplexID,a)).orientation = orient;
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

/**
 * @brief      Calculates the orientation of simplices of a hole
 *
 * @param      mesh       The mesh
 * @param[in]  <unnamed>  { parameter_description }
 *
 * @return     The hole orientation.
 */
bool computeHoleOrientation(SurfaceMesh &mesh, const std::vector<SurfaceMesh::SimplexID<2> > &edgeList);

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
 * @brief      Traditional barycenter smooth.
 *
 * @param      mesh      SurfaceMesh to manipulate.
 * @param[in]  vertexID  SimplexID of the vertex to move.
 */
void barycenterVertexSmooth(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID);

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
        if ((*edgeID).selected == true){
            if (!ignoredEdges.count(edgeID))
            {
                auto up = mesh.get_cover(edgeID);
                // The mesh is not a surface mesh...
                if (up.size() > 2)
                {
                    // std::cerr << "This edge participates in more than 2 faces. "
                    //           << "Returning..." << std::endl;
                    throw std::runtime_error("SurfaceMesh is not pseudomanifold. Found an edge connected to more than 2 faces.");
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

                // TODO: (5) Change to check link condition
                // Check if the triangles make a wedge shape. Flipping this can
                // cause knife faces.
                auto f1   = (shared.first - notShared.first)^(shared.second - notShared.first);
                auto f2   = (shared.first - notShared.second)^(shared.second - notShared.second);
                auto f3   = (notShared.first - shared.first)^(notShared.second - shared.first);
                auto f4   = (notShared.first - shared.second)^(notShared.second - shared.second);
                auto area = std::pow(std::sqrt(f1|f1) + std::sqrt(f2|f2), 2);
                auto areaFlip = std::pow(std::sqrt(f3|f3) + std::sqrt(f4|f4), 2);
                double edgeFlipCriterion = 1.001;
                if (areaFlip/area > edgeFlipCriterion) continue; // If area changes by a lot then continue

                // Check the flip using user function
                if (checkFlip(mesh, edgeID))
                {
                    *iter++ = edgeID;   // Insert into edges to flip

                    // The local topology will be changed by edge flip.
                    // Don't flip edges which share a common face.
                    std::set<SurfaceMesh::SimplexID<2> > tmpIgnored;
                    kneighbors(mesh, edgeID, 3, tmpIgnored);
                    ignoredEdges.insert(tmpIgnored.begin(), tmpIgnored.end());
                    // neighbors(mesh, edgeID, std::inserter(ignoredEdges, ignoredEdges.end()));
                }
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

void normalSmoothH(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, double k);

/**
 * @brief      Traverse mesh and find holes
 *
 * @param[in]  mesh      SurfaceMesh
 * @param      holeList  List of list of edge rings to store holes in
 */
void findHoles(const SurfaceMesh &mesh,
    std::vector<std::vector<SurfaceMesh::SimplexID<2>>>& holeList);

/**
 * @brief      Recursive DFS of edges defining hole ring.
 *
 * @param[in]  mesh                SurfaceMesh of interest
 * @param      unvisitedBdryEdges  Set of unvisited boundary edges
 * @param      bdryRing            Stack of ordered boundary edges
 *
 * @return     True if hole ring found. False otherwise.
 */
bool findHoleHelper(const SurfaceMesh &mesh,
        std::set<SurfaceMesh::SimplexID<2>>& unvisitedBdryEdges,
        std::vector<SurfaceMesh::SimplexID<1>>& visitedVerts,
        std::vector<SurfaceMesh::SimplexID<2>>& bdryRing);

void edgeRingToVertices(const SurfaceMesh &mesh,
        std::vector<SurfaceMesh::SimplexID<2>>& edgeRing,
        std::back_insert_iterator<std::vector<SurfaceMesh::SimplexID<1>>> iter);

void triangulateHoleHelper(SurfaceMesh &mesh,
        std::vector<SurfaceMesh::SimplexID<1>> &boundary,
        const Face &fdata,
        std::back_insert_iterator<std::vector<SurfaceMesh::SimplexID<2>>> iter);

/**
 * @brief      Recursively triangulate hole by connecting lowest valence
 *             vertices.
 *
 * @param      mesh      Surface mesh
 * @param      boundary  List of boundary edges
 * @param[in]  fdata     Data to store on each face
 * @param[in]  iter      Back inserter to store new edges and boundary edges
 */
void triangulateHole(SurfaceMesh &mesh,
        std::vector<SurfaceMesh::SimplexID<1>> &sortedVerts,
        const Face &fdata,
        std::vector<SurfaceMesh::SimplexID<2>>  &edgeList);

} // end namespace surfacemesh_detail

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
void printQualityInfo(const std::string &filename, const SurfaceMesh &mesh);
void generateHistogram(const SurfaceMesh &mesh);
std::tuple<double, double, int, int> getMinMaxAngles(const SurfaceMesh& mesh,
    double maxMinAngle, double minMaxAngle);
double getArea(const SurfaceMesh &mesh);
double getArea(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID);
double getArea(Vertex a, Vertex b, Vertex c);
double getVolume(const SurfaceMesh &mesh);
bool hasHole(const SurfaceMesh &mesh);

std::size_t getValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID);

double getMeanCurvature(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID);
double getGaussianCurvature(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID);


// These exist for the a potential python interface
void translate(SurfaceMesh &mesh, Vector v);
void translate(SurfaceMesh &mesh, double dx, double dy, double dz);
void scale(SurfaceMesh &mesh, Vector v);
void scale(SurfaceMesh &mesh, double sx, double sy, double sz);
void scale(SurfaceMesh &mesh, double s);

std::pair<Vector, double> getCenterRadius(SurfaceMesh &mesh);
void centeralize(SurfaceMesh &mesh);

/**
 * @brief      Smooth the surface mesh
 *
 * @param      mesh            The mesh
 * @param[in]  maxIter         The maximum iterator
 * @param[in]  preserveRidges  The preserve ridges
 * @param[in]  verbose         The verbose
 */
void smoothMesh(SurfaceMesh &mesh, int maxIter, bool preserveRidges, bool verbose);

/**
 * @brief      Coarsens the mesh by selecting vertices first.
 *
 * @param      mesh         The mesh
 * @param[in]  coarseRate   The coarse rate
 * @param[in]  flatRate     The flat rate
 * @param[in]  denseWeight  The dense weight
 */
void coarse(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight);

/**
 * @brief      Coarsen the mesh in one loop. Faster but uses more memory.
 *
 * @param      mesh         The mesh
 * @param[in]  coarseRate   The coarse rate
 * @param[in]  flatRate     The flat rate
 * @param[in]  denseWeight  The dense weight
 */
void coarseIT(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight);

void normalSmooth(SurfaceMesh &mesh);

void fillHoles(SurfaceMesh &mesh);

void flipNormals(SurfaceMesh &mesh);

/**
 * @brief      Refine the mesh by quadrisection of faces
 *
 * @param      mesh  The mesh
 */
std::unique_ptr<SurfaceMesh> refineMesh(const SurfaceMesh &mesh);


/**
 * @brief      Create a triangulated octahedron
 *
 * @param[in]  order  Number of times to subdivide octahedron
 *
 * @return     Pointer to resulting mesh
 */
std::unique_ptr<SurfaceMesh> sphere(int order);


/**
 * @brief      Create a triangulated cube
 *
 * @param[in]  order  Number of times to subdivide cube
 *
 * @return     Pointer to resulting mesh
 */
std::unique_ptr<SurfaceMesh> cube(int order);

