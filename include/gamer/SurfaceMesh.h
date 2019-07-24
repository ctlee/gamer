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

/**
 * @file SurfaceMesh.h
 * @brief Surface mesh functionality and definition
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

#include "gamer/Vertex.h"


/// Namespace for all things gamer
namespace gamer
{

/**
 * @brief      Type for containing root metadata
 */
struct SMGlobal
{
    /// Domain marker to be used when tetrahedralizing.
    int   marker;
    /// Volume constraint of the tetrahedralized domain.
    float volumeConstraint;
    /// flag that determines if the volume constraint is used.
    bool  useVolumeConstraint;
    /// Flag that determines if the mesh represents a hole or not
    bool  ishole;

    /**
     * @brief      Default constructor
     *
     * @param[in]  marker               Global marker value
     * @param[in]  volumeConstraint     Value of volume constraint to use
     * @param[in]  useVolumeConstraint  Whether or not a volume constraint
     *                                  should be applied when tetrahedralizing
     * @param[in]  ishole               Is this domain a hole?
     */
    SMGlobal(int marker = -1, float volumeConstraint = -1, bool useVolumeConstraint = false, bool ishole = false) :
        marker(marker), volumeConstraint(volumeConstraint), useVolumeConstraint(useVolumeConstraint), ishole(ishole) {}
};

struct SMVertex : Vertex {
    /// Cached normal vector
    Vector normal;
    using Vertex::Vertex;
};

/**
 * @brief      Edge data
 */
struct SMEdge{
    /// Selection status of the edge
    bool selected;

    /// Default constructor constructs unselected edge
    SMEdge() : SMEdge(0) {}

    /**
     * @brief      Overload constructor constructs edge with
     *
     * @param[in]  select  Selection status
     */
    SMEdge(bool select) : selected(select) {}
};

/**
 * @brief      Properties that Faces should have
 */
struct SMFaceProperties
{
    int  marker;   /**< @brief Marker */
    bool selected; /**< @brief Selection flag */
    Vector normal; /**< @brief cached normal of face */

    /**
     * @brief      Face properties constructor
     *
     * @param[in]  marker    The marker
     * @param[in]  selected  The selected
     */
    SMFaceProperties(int marker, bool selected) :
        marker(marker), selected(selected) {}
};

/**
 * @brief      SMFace object
 */
struct SMFace : casc::Orientable, SMFaceProperties
{
    /// Default constructor
    SMFace() : SMFace(Orientable{0}, SMFaceProperties{-1, false}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  marker    Marker value
     * @param[in]  selected  Selection status
     */
    SMFace(int marker, bool selected) : SMFace(Orientable{0}, SMFaceProperties{marker, selected}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  orient    Orientation of the simplex
     * @param[in]  marker    Marker value
     * @param[in]  selected  Selection status
     */
    SMFace(int orient, int marker, bool selected) : SMFace(Orientable{orient}, SMFaceProperties{marker, selected}) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  orient  Orientable object
     * @param[in]  prop    Properties of a face
     */
    SMFace(Orientable orient, SMFaceProperties prop)
        : Orientable(orient), SMFaceProperties(prop)
    {}

    /**
     * @brief      Print operator overload
     *
     * @param      output  The output
     * @param[in]  f       Face to print
     *
     * @return     Output stream
     */
    friend std::ostream &operator<<(std::ostream &output, const SMFace &f)
    {
        output << "SMFace("
               << "m:" << f.marker
               << ";sel:" << std::boolalpha << f.selected
               << ";o:" << f.orientation << ")";
        return output;
    }

    /**
     * @brief      Returns a string representation of the object.
     *
     * @return     String representation of the object.
     */
    std::string to_string() const
    {
        std::ostringstream output;
        output << *this;
        return output.str();
    }

};


/// @cond detail
namespace surfmesh_detail
{
/**
 * @brief      A helper struct containing the traits/types in the simplicial
 *             complex
 */
struct surfmesh_traits
{
    /// The index type
    using KeyType = int;
    /// The types of each node
    using NodeTypes = util::type_holder<SMGlobal, SMVertex, SMEdge, SMFace>;
    /// The types of each edge
    using EdgeTypes = util::type_holder<casc::Orientable, casc::Orientable, casc::Orientable>;
};
} // end namespace surfmesh_detail
/// @endcond

/// Surface Mesh Object
using SurfaceMesh = casc::simplicial_complex<surfmesh_detail::surfmesh_traits>;

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

/// @cond detail
/// Namespace for surface mesh detail functions
namespace surfacemesh_detail
{
/**
 * @brief      Remove a vertex from mesh and triangulate the resulting hole.
 *
 * @param      mesh      The mesh
 * @param[in]  vertexID  The vertex id
 */
void decimateVertex(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, std::size_t rings = 2);

/**
 * @brief      Computes the local structure tensor
 *
 * @param[in]  mesh      Surface mesh of interest
 * @param[in]  vertexID  Vertex of interest
 * @param[in]  rings     Number of neighborhood rings to consider
 *
 * @return     The local structure tensor.
 */
tensor<double, 3, 2> computeLocalStructureTensor(
    const SurfaceMesh              &mesh,
    const SurfaceMesh::SimplexID<1> vertexID,
    const int                       rings);


/**
 * @brief      Computes the local structure tensor from cached normals
 *
 * @param[in]  mesh      Surface mesh of interest
 * @param[in]  vertexID  Vertex of interest
 * @param[in]  rings     Number of neighborhood rings to consider
 *
 * @return     The local structure tensor.
 */
tensor<double, 3, 2> computeLSTFromCache(
    const SurfaceMesh              &mesh,
    const SurfaceMesh::SimplexID<1> vertexID,
    const int                       rings);


/**
 * @brief      Terminal case
 *
 * @param[in]  mesh       The mesh
 * @param[in]  origin     The origin
 * @param[in]  curr       The curr
 *
 * @tparam     dimension  Vector space dimension
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
 * @tparam     level      Level of the complex currently traversed
 * @tparam     dimension  Vector dimension
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
 * @brief      General template for function to initialize the local orientation
 *             of a mesh.
 *
 * @tparam     K     Class template intended for storing level
 */
template <class K>
struct initLocalOrientation {};

/**
 * @brief      Specialization to initialize the local orientation of a mesh
 *
 * @tparam     k     The current level traversed
 */
template <std::size_t k>
struct initLocalOrientation<std::integral_constant<std::size_t, k> >
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
    static void apply(SurfaceMesh          &mesh,
                      const std::set<int> &&names,
                      Iterator              begin,
                      Iterator              end)
    {
        std::vector<SurfaceMesh::SimplexID<k+1> > next;
        for (auto curr = begin; curr != end; ++curr)
        {
            auto currSimplexID = *curr;
            for (auto a : mesh.get_cover(currSimplexID))
            {
                // Look for key a in names
                auto find = names.find(a);
                if (find != names.end())
                {
                    next.push_back(mesh.get_simplex_up(currSimplexID, a));

                    int orient = 1;
                    for (auto b : mesh.get_name(currSimplexID))
                    {
                        if (a > b)
                        {
                            if (a > b)
                            {
                                orient *= -1;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    (*mesh.get_edge_up(currSimplexID, a)).orientation = orient;
                }
            }
        }
        initLocalOrientation<std::integral_constant<std::size_t, k+1> >::apply(mesh, std::move(names), next.begin(), next.end());
    }
};

/**
 * @brief      Terminal case. The top level does not need to be initialized.
 */
template <>
struct initLocalOrientation<std::integral_constant<std::size_t, SurfaceMesh::topLevel> > {
    template <typename Iterator>
    static void apply(SurfaceMesh &mesh, const std::set<int> &&names, Iterator begin, Iterator end){}
};

/**
 * @brief      Calculates the orientation of faces around a local region
 *
 * @param      mesh      Surface mesh of interest
 * @param[in]  edgeList  List of edges around which the orientation fill is
 *                       restricted.
 *
 * @return     True if the hole was orientable, False if an inconsistent
 *             orientation was found.
 */
bool computeLocalOrientation(SurfaceMesh &mesh, const std::vector<SurfaceMesh::SimplexID<2> > &edgeList);

/**
 * @brief      Compute the eigenvalues of a 3x3 matrix.
 *
 * @param[in]  mat   3x3 matrix to compute eigenvalues of.
 *
 * @return     Returns an Eigen::SelfAdjointEigenSolver containing the results
 */
Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double, 3, 2> &mat);

/**
 * @brief      Smooth the vertex according to Section 2.2.2 of GAMer paper.
 *
 * @param[out] mesh      SurfaceMesh of interest.
 * @param[in]  vertexID  SimplexID of the vertex to move.
 * @param[in]  rings     The number of neighbor rings to use to compute the LST.
 */
void weightedVertexSmooth(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, int rings);

Vector weightedVertexSmoothCache(SurfaceMesh &mesh,
                               SurfaceMesh::SimplexID<1> vertexID,
                               std::size_t rings);

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

void edgeFlipCache(SurfaceMesh &mesh, SurfaceMesh::SimplexID<2> edgeID);

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
                     std::function<bool(const SurfaceMesh &, const SurfaceMesh::SimplexID<2> &)> &&checkFlip,
                     Inserter iter)
{
    casc::NodeSet<SurfaceMesh::SimplexID<2> > ignoredEdges;

    for (auto edgeID : mesh.get_level_id<2>())
    {
        if ((*edgeID).selected == true)
        {
            if (!ignoredEdges.count(edgeID))
            {
                auto up = mesh.get_cover(edgeID);
                // The mesh is not a surface mesh...
                if (up.size() > 2)
                {
                    // std::cerr << "This edge participates in more than 2
                    // faces. "
                    //           << "Returning..." << std::endl;
                    throw std::runtime_error("SurfaceMesh is not pseudomanifold. Found an edge connected to more than 2 faces.");
                }
                else if (up.size() < 2) // Edge is a boundary
                {
                    // std::cerr << "This edge participates in fewer than 2
                    // faces. "
                    //           << "Returning..." << std::endl;
                    continue;
                }

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

                // Check the flip using user function
                if (checkFlip(mesh, edgeID))
                {
                    *iter++ = edgeID;   // Insert into edges to flip

                    // The local topology will be changed by edge flip.
                    // Don't flip edges which share a common face.
                    std::set<SurfaceMesh::SimplexID<2> > tmpIgnored;
                    kneighbors(mesh, edgeID, 3, tmpIgnored);
                    ignoredEdges.insert(tmpIgnored.begin(), tmpIgnored.end());

                    // Local neighborhood append. Larger neighborhood selected
                    // above appears to work better...
                    // neighbors(mesh, edgeID, std::inserter(ignoredEdges,
                    // ignoredEdges.end()));
                }
            }
        }
    }
}

/**
 * @brief      Check if an edge should be flipped
 *
 * @param[in]  mesh            SurfaceMesh to operate on.
 * @param[in]  preserveRidges  Whether or not to try preserving ridges.
 * @param[in]  edgeID          The edge id
 * @param[in]  checkFlip       Functor specifying flip criteria
 *
 * @return     True if edge should be flipped
 */
bool checkEdgeFlip(const SurfaceMesh &mesh,
                     bool preserveRidges,
                     SurfaceMesh::SimplexID<2> edgeID,
                     std::function<bool(const SurfaceMesh &, const SurfaceMesh::SimplexID<2> &)> &&checkFlip
                );
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
int checkFlipValenceExcess(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID);

/**
 * @brief      Apply normal smoothing to a region around a vertex
 *
 * @param      mesh      The mesh
 * @param[in]  vertexID  Vertex of interest
 * @param[in]  k         Anisotropic smoothing factor
 */
void normalSmoothH(SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID, const double k);

/**
 * @brief      Traverse mesh and find holes
 *
 * @param[in]  mesh      SurfaceMesh
 * @param      holeList  List of list of edge rings to store holes in
 */
void findHoles(const SurfaceMesh                                     &mesh,
               std::vector<std::vector<SurfaceMesh::SimplexID<2> > > &holeList);

/**
 * @brief      Sort a list of boundary edges into ring order.
 *
 *             Recursive DFS of edges to find the ring ordering. Note that this
 *             function assumes that the existing content in @p visitedVerts and
 *             @p BdryRings are sorted in ring order. Also, @p
 *             unvisitedBdryEdges can contain edges from other holes. These
 *             other edges will remain in @p univistedBdryEdges post operation.
 *             When first called, if \p visitedVerts and \p bdryRing are empty,
 *             the first \p univistedBdryEdge will be used to seed the hole.
 *
 * @param[in]  mesh                SurfaceMesh of interest
 * @param[in]  unvisitedBdryEdges  Set of unvisited boundary edges
 * @param      visitedVerts        Visited vertices in order
 * @param[out] bdryRing            Stack of ordered boundary edges
 *
 * @return     True if hole ring found. False otherwise.
 */
bool orderBoundaryEdgeRing(const SurfaceMesh                       &mesh,
                           std::set<SurfaceMesh::SimplexID<2> >    &unvisitedBdryEdges,
                           std::vector<SurfaceMesh::SimplexID<1> > &visitedVerts,
                           std::vector<SurfaceMesh::SimplexID<2> > &bdryRing);

void edgeRingToVertices(const SurfaceMesh                                                  &mesh,
                        std::vector<SurfaceMesh::SimplexID<2> >                            &edgeRing,
                        std::back_insert_iterator<std::vector<SurfaceMesh::SimplexID<1> > > iter);

void triangulateHoleHelper(SurfaceMesh                                                        &mesh,
                           std::vector<SurfaceMesh::SimplexID<1> >                            &boundary,
                           const SMFace                                                       &fdata,
                           std::back_insert_iterator<std::vector<SurfaceMesh::SimplexID<2> > > iter);

/**
 * @brief      Recursively triangulate hole by connecting lowest valence
 *             vertices.
 *
 * @param      mesh      Surface mesh
 * @param      sortedVerts  List of boundary edges
 * @param[in]  fdata     Data to store on each face
 * @param[in]  edgeList  Back inserter to store new edges and boundary edges
 */
void triangulateHole(SurfaceMesh                             &mesh,
                     std::vector<SurfaceMesh::SimplexID<1> > &sortedVerts,
                     const SMFace                            &fdata,
                     std::vector<SurfaceMesh::SimplexID<2> > &edgeList);

template <typename Complex>
struct CopyHelper
{
    using SimplexSet = typename casc::SimplexSet<Complex>;
    using KeyType = typename Complex::KeyType;

    template <std::size_t k>
    static void apply(Complex          &before,
                      Complex          &after,
                      const SimplexSet &S)
    {
        for (auto sID : casc::get<k>(S))
        {
            auto name = before.get_name(sID);
            auto data = *sID;
            after.insert(name, data);
        }
    }
};
} // end namespace surfacemesh_detail
/// @endcond

/**
 * @brief      Reads in a GeomView OFF file
 *
 * @param[in]  filename  Filename of file of interest
 *
 * @return     Returns a unique_ptr to the SurfaceMesh
 */
std::unique_ptr<SurfaceMesh> readOFF(const std::string &filename);


/**
 * @brief      Write the SurfaceMesh to file in OFF format.
 *
 * @param[in]  filename  The filename to write to.
 * @param[in]  mesh      SurfaceMesh of interest.
 */
void writeOFF(const std::string &filename, const SurfaceMesh &mesh);

/**
 * @brief      Reads an obj file.
 *
 * @param[in]  filename  The filename
 *
 * @return     Unique pointer to SurfaceMesh
 */
std::unique_ptr<SurfaceMesh> readOBJ(const std::string &filename);


/**
 * @brief      Writes a mesh to obj file format.
 *
 * @param[in]  filename  The filename to write out to
 * @param[in]  mesh      Surface mesh to output
 */
void writeOBJ(const std::string &filename, const SurfaceMesh &mesh);

/**
 * @brief      Pretty print the mesh.
 *
 * @param[in]  mesh  The mesh to print
 */
void print(const SurfaceMesh &mesh);

/**
 * @brief      { function_description }
 *
 * @param[in]  filename  The filename
 * @param[in]  mesh      The mesh
 */
void printQualityInfo(const std::string &filename, const SurfaceMesh &mesh);

/**
 * @brief      { function_description }
 *
 * @param[in]  mesh  The mesh
 */
void generateHistogram(const SurfaceMesh &mesh);

/**
 * @brief      Gets the minimum maximum angles.
 *
 * @param[in]  mesh         The mesh
 * @param[in]  maxMinAngle  The maximum minimum angle
 * @param[in]  minMaxAngle  The minimum maximum angle
 *
 * @return     The minimum maximum angles.
 */
std::tuple<double, double, int, int> getMinMaxAngles(const SurfaceMesh &mesh,
                                                     double maxMinAngle, double minMaxAngle);

/**
 * @brief      Gets the area.
 *
 * @param[in]  mesh  The mesh
 *
 * @return     The area.
 */
double getArea(const SurfaceMesh &mesh);

/**
 * @brief      Gets the area.
 *
 * @param[in]  mesh    The mesh
 * @param[in]  faceID  The face id
 *
 * @return     The area.
 */
double getArea(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID);

/**
 * @brief      Gets the area.
 *
 * @param[in]  a     { parameter_description }
 * @param[in]  b     { parameter_description }
 * @param[in]  c     { parameter_description }
 *
 * @return     The area.
 */
double getArea(Vertex a, Vertex b, Vertex c);

REAL getArea(std::array<Vertex, 3> t);

/**
 * @brief      Gets the volume.
 *
 * @param[in]  mesh  The mesh
 *
 * @return     The volume.
 */
double getVolume(const SurfaceMesh &mesh);

/**
 * @brief      Determines if a surface mesh contains holes.
 *
 * @param[in]  mesh  The mesh
 *
 * @return     True if has hole, False otherwise.
 */
bool hasHole(const SurfaceMesh &mesh);

/**
 * @brief      Gets the number of edges connected to a vertex.
 *
 * @param[in]  mesh      The mesh
 * @param[in]  vertexID  The vertex id
 *
 * @return     The valence.
 */
inline std::size_t getValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID){
    return mesh.get_cover(vertexID).size();
}

/**
 * @brief      Gets the mean curvature.
 *
 * @param[in]  mesh      The mesh
 * @param[in]  vertexID  The vertex id
 *
 * @return     The mean curvature.
 */
double getMeanCurvature(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID);
/**
 * @brief      Gets the gaussian curvature.
 *
 * @param[in]  mesh      The mesh
 * @param[in]  vertexID  The vertex id
 *
 * @return     The gaussian curvature.
 */
double getGaussianCurvature(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID);


std::tuple<REAL*, REAL*, REAL*, REAL*, std::map<typename SurfaceMesh::KeyType,typename SurfaceMesh::KeyType>>
computeCurvatures(const SurfaceMesh &mesh);

//
// @param      mesh  The mesh
// @param[in]  v     Displacement vector
//
void translate(SurfaceMesh &mesh, Vector v);
/**
 * @brief      Translate the mesh
 *
 * @param      mesh  The mesh
 * @param[in]  dx    Distance to move in x direction
 * @param[in]  dy    Distance to move in y direction
 * @param[in]  dz    Distance to move in z direction
 */
void translate(SurfaceMesh &mesh, double dx, double dy, double dz);
/**
 * @brief      Scale mesh anisotropically
 *
 * @param      mesh  The mesh
 * @param[in]  v     Vector with components representing the anisotropic scaling
 *                   factors.
 */
void scale(SurfaceMesh &mesh, Vector v);
/**
 * @brief      Scale a mesh anisotropically
 *
 * @param      mesh  The mesh
 * @param[in]  sx    Scale factor for x axis
 * @param[in]  sy    Scale factor for y axis
 * @param[in]  sz    Scale factor for z axis
 */
void scale(SurfaceMesh &mesh, double sx, double sy, double sz);
/**
 * @brief      Scale a mesh isotropically
 *
 * @param      mesh  The mesh
 * @param[in]  s     Scale factor
 */
void scale(SurfaceMesh &mesh, double s);

/**
 * @brief      Compute the center and radius of the mesh
 *
 * @param      mesh  The mesh
 *
 * @return     The center and radius.
 */
std::pair<Vector, double> getCenterRadius(SurfaceMesh &mesh);
/**
 * @brief      Center the mesh on its center of mass
 *
 * @param      mesh  The mesh
 */
void centeralize(SurfaceMesh &mesh);

/**
 * @brief      Smooth the surface mesh
 *
 * @param      mesh            The mesh
 * @param[in]  maxIter         The maximum iterator
 * @param[in]  preserveRidges  The preserve ridges
 * @param[in]  verbose         The verbose
 */
void smoothMesh(SurfaceMesh &mesh, int maxIter, bool preserveRidges, bool verbose, std::size_t rings = 2);

/**
 * @brief      Coarsens the mesh
 *
 * @param      mesh         The mesh
 * @param[in]  coarseRate   The coarse rate
 * @param[in]  flatRate     The flat rate
 * @param[in]  denseWeight  The dense weight
 */
void coarse(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight, std::size_t rings = 2);

/**
 * @brief      Coarsens the mesh by selecting vertices first.
 *
 * @param      mesh       The mesh
 * @param[in]  threshold  The threshold
 * @param[in]  weight     The weight
 * @param[in]  rings      The rings
 */
void coarse_dense(SurfaceMesh &mesh, REAL threshold, REAL weight, std::size_t rings = 2);

/**
 * @brief      Coarsens flat regions by LST analysis
 *
 * @param      mesh       The mesh
 * @param[in]  threshold  The threshold
 * @param[in]  weight     The weight
 * @param[in]  rings      The rings
 */
void coarse_flat(SurfaceMesh &mesh, REAL threshold, REAL weight, std::size_t rings = 2);

/**
 * @brief      Perform smoothing of the mesh normals
 *
 * @param      mesh  The mesh
 */
void normalSmooth(SurfaceMesh &mesh);

/**
 * @brief      Fill holes in the mesh
 *
 * @param      mesh  The mesh
 */
void fillHoles(SurfaceMesh &mesh);

/**
 * @brief      Flip the normals of the mesh
 *
 * @param      mesh  The mesh
 */
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

/**
 * @brief      Split connected surfaces from a single mesh.
 *
 *             Note that this function creates new meshes from surfaces
 *
 * @param      mesh  The mesh
 *
 * @return     Vector of surface meshes
 */
std::vector<std::unique_ptr<SurfaceMesh> > splitSurfaces(SurfaceMesh &mesh);

void cacheNormals(SurfaceMesh &mesh);
} // end namespace gamer
