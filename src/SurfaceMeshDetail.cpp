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


#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>

#include <casc/casc>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gamer/SurfaceMesh.h"
#include "gamer/Vertex.h"

/// Namespace for all things gamer
namespace gamer
{
namespace surfacemesh_detail
{
tensor<double, 3, 2> computeLocalStructureTensor(const SurfaceMesh              &mesh,
                                                 const SurfaceMesh::SimplexID<1> vertexID,
                                                 const int                       rings)
{
    // Set of neighbors
    std::set<SurfaceMesh::SimplexID<1> > nbors;
    // Get list of neighbors
    casc::kneighbors_up(mesh, vertexID, rings, nbors);
    // local structure tensor
    tensor<double, 3, 2> lst = tensor<double, 3, 2>();
    for (SurfaceMesh::SimplexID<1> nid : nbors)
    {
        auto norm = getNormal(mesh, nid);           // Get Vector normal
        auto mag  = std::sqrt(norm|norm);
        if (mag <= 0)
        {
            continue;
        }
        norm /= mag;               // normalize
        lst  += norm*norm;         // tensor product
    }

    // Print the LST nicely
    // std::cout << "LST:\n" << std::endl;
    // for(int i = 0; i < 3; ++i){
    //     for(int j = 0; j < 3; ++j)
    //         std::cout << lst.get(i,j) << " ";
    //     std::cout << "\n";
    // }
    return lst;
}

tensor<double, 3, 2> computeLSTFromCache(const SurfaceMesh              &mesh,
                                                 const SurfaceMesh::SimplexID<1> vertexID,
                                                 const int                       rings)
{
    // Set of neighbors
    std::set<SurfaceMesh::SimplexID<1> > nbors;
    // Get list of neighbors
    casc::kneighbors_up(mesh, vertexID, rings, nbors);
    // local structure tensor
    tensor<double, 3, 2> lst = tensor<double, 3, 2>();
    for (SurfaceMesh::SimplexID<1> vID : nbors)
    {
        auto norm = (*vID).normal;           // Get Vector normal
        lst  += norm*norm;         // tensor product
    }
    return lst;
}

void decimateVertex(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, std::size_t rings)
{
    // TODO: (10) Come up with a better scheme
    // Pick an arbitrary face's data
    auto fdata = **mesh.up(std::move(mesh.up(vertexID))).begin();
    fdata.orientation = 0; // Reset the orientation accordingly

    // Compute and backup ring of vertices prior to vertex removal
    std::set<SurfaceMesh::SimplexID<1> > boundary;
    casc::neighbors_up(mesh, vertexID, std::inserter(boundary, boundary.end()));
    std::set<SurfaceMesh::SimplexID<1> > backupBoundary(boundary);

    // Remove the vertex
    mesh.remove(vertexID);

    // Sort vertices into 'ring' order
    std::vector<SurfaceMesh::SimplexID<1> > sortedVerts;
    std::set<int> bNames; // boundary names
    std::vector<SurfaceMesh::SimplexID<2> > edgeList;

    auto it = boundary.begin();
    int  firstName = mesh.get_name(*it)[0];
    while (boundary.size() > 0)
    {
        std::vector<SurfaceMesh::SimplexID<1> > nbors;
        auto currID = *it; // Get SimplexID

        int  currName = mesh.get_name(currID)[0];
        bNames.insert(currName);

        bool success = false; // Flag to track success
        // Push current into list of sorted vertices.
        std::move(it, std::next(it), std::back_inserter(sortedVerts));
        boundary.erase(it); // Erase current from boundary

        // If nothing is left in the boundary, check that the current
        // vertex has an edge to the first vertex.
        if (boundary.size() == 0)
        {
            if (mesh.exists({currName, firstName}))
            {
                edgeList.push_back(mesh.get_simplex_up(currID, firstName));
                break; // Break out of while loop
            }
        }
        else
        {
            // Get neighbors and search for next vertex
            casc::neighbors_up(mesh, currID, std::back_inserter(nbors));
            for (auto nbor : nbors)
            {
                auto result = boundary.find(nbor);
                if (result != boundary.end())
                {
                    // Check that the edge is a boundary
                    auto tmp = mesh.get_simplex_up(*result, currName);
                    if (mesh.get_cover(tmp).size() == 1)
                    {
                        it = result;
                        edgeList.push_back(tmp);
                        success = true;
                        break; // leave for loop
                    }
                }
            }
        }
        // The ring isn't really a ring
        if (!success)
        {
            throw std::runtime_error("ERROR(coarse): Hole ring is not closed. "
                                     "Please contact the developers with this error.");
        }
    }

    triangulateHole(mesh, sortedVerts, fdata, edgeList);

    // Smooth vertices around the filled hole
    for (auto v : backupBoundary)
    {
        weightedVertexSmooth(mesh, v, rings);
    }
}

void triangulateHoleHelper(SurfaceMesh                                                        &mesh,
                           std::vector<SurfaceMesh::SimplexID<1> >                            &boundary,
                           const SMFace                                                       &fdata,
                           std::back_insert_iterator<std::vector<SurfaceMesh::SimplexID<2> > > edgeListIT)
{
    // Terminal case
    if (boundary.size() == 3)
    {
        // create the face
        auto a = mesh.get_name(boundary[0])[0];
        auto b = mesh.get_name(boundary[1])[0];
        auto c = mesh.get_name(boundary[2])[0];
        mesh.insert({a, b, c}, fdata);
        return;
    }

    // Construct a sorted vector of pairs... (valence, vertexID) by valence
    std::vector<std::pair<int, SurfaceMesh::SimplexID<1> > > list;
    for (auto vertexID : boundary)
    {
        list.push_back(std::make_pair(getValence(mesh, vertexID), vertexID));
    }
    std::sort(list.begin(), list.end(), [](
                  const std::pair<int, SurfaceMesh::SimplexID<1> > &lhs,
                  const std::pair<int, SurfaceMesh::SimplexID<1> > &rhs){
                return lhs.first < rhs.first;
            });

    SurfaceMesh::SimplexID<1> v1, v2;
    v1 = list[0].second;

    // Find v1 and rotate so that it is first for easy splitting later
    auto v1it =  std::find(boundary.begin(), boundary.end(), v1);
    std::rotate(boundary.begin(), v1it, boundary.end());

    // Get the next lowest valence vertex
    for (auto it = ++list.begin(); it != list.end(); ++it)
    {
        v2 = (*it).second;
        // Check that it is not already connected to v1
        if (v2 != boundary[1] && v2 != boundary.back())
        {
            break;
        }
    }
    // Insert new edge
    auto a = mesh.get_name(v1)[0];
    auto b = mesh.get_name(v2)[0];
    mesh.insert({a, b});
    // TODO: (0) Get the simplex from insert.
    *edgeListIT++ = mesh.get_simplex_up({a, b});

    auto v2it =  std::find(boundary.begin(), boundary.end(), v2);
    std::vector<SurfaceMesh::SimplexID<1> > other;
    other.push_back(v2);
    std::move(v2it+1, boundary.end(), std::back_inserter(other));
    boundary.erase(v2it+1, boundary.end());

    other.push_back(v1);

    // Recurse to fill sub-holes
    triangulateHoleHelper(mesh, boundary, fdata, edgeListIT);
    triangulateHoleHelper(mesh, other, fdata, edgeListIT);
}

void triangulateHole(SurfaceMesh                             &mesh,
                     std::vector<SurfaceMesh::SimplexID<1> > &sortedVerts,
                     const SMFace                            &fdata,
                     std::vector<SurfaceMesh::SimplexID<2> > &holeEdges)
{
    // backup sortedVerts
    std::vector<SurfaceMesh::SimplexID<1> > backupVerts = sortedVerts;

    // Triangulate the hole
    triangulateHoleHelper(mesh, sortedVerts, fdata, std::back_inserter(holeEdges));

    std::set<int> keys;
    for (auto vertexID : backupVerts)
    {
        keys.insert(mesh.get_name(vertexID)[0]);
    }

    initLocalOrientation<std::integral_constant<size_t, 1> >::apply(mesh, std::move(keys), backupVerts.begin(), backupVerts.end());

    if (!computeLocalOrientation(mesh, holeEdges))
    {
        throw std::runtime_error("ERROR(triangulateHole): Mesh became non-orientable");
    }
}

bool computeLocalOrientation(SurfaceMesh                                   &mesh,
                             const std::vector<SurfaceMesh::SimplexID<2> > &edgeList)
{
    std::deque<SurfaceMesh::SimplexID<2> > frontier;
    std::set<SurfaceMesh::SimplexID<2> >   visited;
    bool orientable = true;

    for (auto outer : edgeList)
    {
        if (visited.find(outer) == visited.end())
        {
            frontier.push_back(outer);

            while (!frontier.empty())
            {
                auto curr = frontier.front();

                if (visited.find(curr) == visited.end())
                {
                    visited.insert(curr);
                    auto w = mesh.get_cover(curr);

                    if (w.size() == 1)
                    {
                        //std::cout << curr << ":" << w[0] << " ~ Boundary" <<
                        // std::endl;
                    }
                    else if (w.size() == 2)
                    {
                        auto &edge0 = *mesh.get_edge_up(curr, w[0]);
                        auto &edge1 = *mesh.get_edge_up(curr, w[1]);

                        auto &node0 = *mesh.get_simplex_up(curr, w[0]);
                        auto &node1 = *mesh.get_simplex_up(curr, w[1]);

                        if (node0.orientation == 0)
                        {
                            if (node1.orientation == 0)
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
                            if (node1.orientation == 0)
                            {
                                node1.orientation = -edge1.orientation * edge0.orientation * node0.orientation;
                            }
                            else
                            {
                                if (edge0.orientation*node0.orientation + edge1.orientation*node1.orientation != 0)
                                {
                                    orientable = false;
                                }
                            }
                        }
                        std::vector<SurfaceMesh::SimplexID<2> > tmp;
                        neighbors_up(mesh, curr, std::back_inserter(tmp));
                        for (auto e : tmp)
                        {
                            if (std::find(edgeList.begin(), edgeList.end(), e) != edgeList.end())
                                frontier.push_back(e);
                        }
                    }
                    else
                    {
                        std::cerr << "ERROR(computeLocalOrientation): Found an edge"
                                  << " connected to " << w.size() << " faces. The SurfaceMesh "
                                  << "is no longer a surface mesh." << std::endl;
                        return false;
                    }
                }
                frontier.pop_front();
            }
        }
    }
    return orientable;
}

Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double, 3, 2> &mat)
{
    Eigen::Map<Eigen::Matrix3d>                    emat(mat.data());
    // TODO: (99) How much optimization can we get from having a persistent
    // eigensolver?
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(emat);
    if (eigensolver.info() != Eigen::Success){
        std::stringstream ss;
        ss << "getEigenvalues has encountered Eigen error " << eigensolver.info();
        throw std::runtime_error(ss.str());
    }
    return eigensolver;
}

void weightedVertexSmooth(SurfaceMesh              &mesh,
                          SurfaceMesh::SimplexID<1> vertexID,
                          int                       rings)
{
    auto   centerName = mesh.get_name(vertexID)[0];
    auto  &center = *vertexID; // get the vertex data

    double sumWeights = 0;
    Vector newPos;

    // Compute the following sum to get the new position
    // \bar{x} = \frac{1}{\sum_{i=1}^{N_2}(\alpha_i+1)}\sum_{i=1}^{N_2}(\alpha_i
    // + 1) x_i
    for (auto edge : mesh.up(vertexID))
    {
        // Get the vertex connected by edge
        auto   edgeName = mesh.get_name(edge);
        Vertex shared   = *mesh.get_simplex_up({(edgeName[0] == centerName) ? edgeName[1] : edgeName[0]});

        // Get the vertices connected to adjacent edge by face
        auto up = mesh.get_cover(edge);
        if (up.size() != 2)
        {
            // Vertex is on a boundary or otherwise...
            // Let's not move it
            return;
        }

        Vertex prev = *mesh.get_simplex_up({up[0]});
        Vertex next = *mesh.get_simplex_up({up[1]});

        // std::cout << mesh.get_simplex_up({up[0]}) << " "
        //             << mesh.get_simplex_up({up[1]}) << " "
        //             << mesh.get_simplex_up({(edgeName[0] == centerName) ?
        // edgeName[1] : edgeName[0]}) << std::endl;

        // std::cout << prev << " " << next << " " << shared << std::endl;
        Vector pS, nS;
        try
        {
            pS = prev - shared;
            normalize(pS);
            nS = next - shared;
            normalize(nS);
        }
        catch (std::exception &e)
        {
            throw std::runtime_error("ERROR: Zero length edge found. "
                                     "weightedVertexSmooth expects no zero length edges.");
        }

        // Bisector of the 'rhombus'
        Vector bisector = pS + nS;
        Vector perpNorm;
        double alpha = (pS|nS) + 1;

        // Check if vertices are colinear
        if (length(bisector) == 0 || alpha < 1e-6 || fabs(alpha-2) < 1e-6)
        {
            perpNorm = pS;
            normalize(perpNorm);
        }
        else
        {
            normalize(bisector);
            // Normal of tangent plane
            Vector tanNorm = cross(pS, nS);
            // Get the perpendicular plane made up of plane normal of
            // bisector
            perpNorm = cross(tanNorm, bisector);
            normalize(perpNorm);
        }

        // Get a reference vector to shared which lies on the plane of interest.
        Vector disp = center - shared;
        Eigen::Map<Eigen::Vector3d> disp_e(disp.data());

        // Perpendicular projector
        auto   perpProj = perpNorm*perpNorm; // tensor product
        // Compute perpendicular component
        Vector perp;
        Eigen::Map<Eigen::Matrix3d> perpProj_e(perpProj.data());
        Eigen::Map<Eigen::Vector3d> perp_e(perp.data());
        perp_e = perpProj_e*disp_e; // matrix (3x3) * vector = vector

        sumWeights += alpha;

        newPos += alpha*(center.position - perp);
    }
    newPos /= sumWeights;
    /**
     * Project to move onto LST eigenvector space and scale by eigenvalues.
     *
     * \bar{x} = x + \sum_{k=1}^3 \frac{1}{1+\lambda_k}((\bar{x} - x)\cdot
     * \vec{e_k})\vec{e_k}
     */
    auto lst = surfacemesh_detail::computeLocalStructureTensor(mesh, vertexID, rings);

    auto eigen_result = getEigenvalues(lst);
    // std::cout << "Eigenvalues(LST): "
    //      << eigen_result.eigenvalues().transpose() << std::endl;
    // std::cout << "Here's a matrix whose columns are eigenvectors of LST \n"
    //      << "corresponding to these eigenvalues:\n"
    //      << eigen_result.eigenvectors() << std::endl;

    newPos -= center.position;  // Vector of old position to new position
    Eigen::Map<Eigen::Vector3d> newPos_e(newPos.data());

    // dot product followed by elementwise-division EQN 4. w is a scale factor.
    auto w = (
        (eigen_result.eigenvectors().transpose()*newPos_e).array()
        / (eigen_result.eigenvalues().array()+1)
        ).matrix();                           // vector 3x1
    newPos_e = eigen_result.eigenvectors()*w; // matrix 3x3 * vector = vector
    center.position += newPos;
}

Vector weightedVertexSmoothCache(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, std::size_t rings){
    auto   centerName = mesh.get_name(vertexID)[0];
    auto  &center = *vertexID; // get the vertex data

    double sumWeights = 0;
    Vector newPos;

    // Compute the following sum to get the new position
    // \bar{x} = \frac{1}{\sum_{i=1}^{N_2}(\alpha_i+1)}\sum_{i=1}^{N_2}(\alpha_i
    // + 1) x_i
    for (auto edge : mesh.up(vertexID))
    {
        // Get the vertex connected by edge
        auto   edgeName = mesh.get_name(edge);
        Vertex shared   = *mesh.get_simplex_up({(edgeName[0] == centerName) ? edgeName[1] : edgeName[0]});

        // Get the vertices connected to adjacent edge by face
        auto cover = mesh.get_cover(edge);
        if (cover.size() != 2)
        {
            // Vertex is on a boundary... Don't move it
            return Vector({0,0,0});
        }
        std::vector<SurfaceMesh::SimplexID<3>> faces;
        mesh.up(edge, std::back_inserter(faces));

        Vertex prev = *mesh.get_simplex_down(faces[0], edgeName);
        Vertex next = *mesh.get_simplex_down(faces[1], edgeName);

        // std::cout << mesh.get_simplex_up({up[0]}) << " "
        //             << mesh.get_simplex_up({up[1]}) << " "
        //             << mesh.get_simplex_up({(edgeName[0] == centerName) ?
        // edgeName[1] : edgeName[0]}) << std::endl;

        // std::cout << prev << " " << next << " " << shared << std::endl;
        Vector pS, nS;
        try
        {
            pS = prev - shared;
            normalize(pS);
            nS = next - shared;
            normalize(nS);
        }
        catch (std::exception &e)
        {
            throw std::runtime_error("ERROR: Zero length edge found. "
                                     "weightedVertexSmooth expects no zero length edges.");
        }

        // Bisector of the 'rhombus'
        Vector bisector = pS + nS;
        Vector perpNorm;
        double alpha = (pS|nS) + 1;

        // Check if vertices are colinear
        if (length(bisector) == 0 || alpha < 1e-6 || fabs(alpha-2) < 1e-6)
        {
            perpNorm = pS;
            normalize(perpNorm);
        }
        else
        {
            normalize(bisector);
            // Normal of tangent plane
            Vector tanNorm = cross(pS, nS);
            // Get the perpendicular plane made up of plane normal of
            // bisector
            perpNorm = cross(tanNorm, bisector);
            normalize(perpNorm);
        }

        // Get a reference vector to shared which lies on the plane of interest.
        Vector disp = center - shared;
        Eigen::Map<Eigen::Vector3d> disp_e(disp.data());

        // Perpendicular projector
        auto   perpProj = perpNorm*perpNorm; // tensor product
        // Compute perpendicular component
        Vector perp;
        Eigen::Map<Eigen::Matrix3d> perpProj_e(perpProj.data());
        Eigen::Map<Eigen::Vector3d> perp_e(perp.data());
        perp_e = perpProj_e*disp_e; // matrix (3x3) * vector = vector

        sumWeights += alpha;

        newPos += alpha*(center.position - perp);
    }
    newPos /= sumWeights;
    /**
     * Project to move onto LST eigenvector space and scale by eigenvalues.
     *
     * \bar{x} = x + \sum_{k=1}^3 \frac{1}{1+\lambda_k}((\bar{x} - x)\cdot
     * \vec{e_k})\vec{e_k}
     */
    auto lst = surfacemesh_detail::computeLocalStructureTensor(mesh, vertexID, rings);

    auto eigen_result = getEigenvalues(lst);
    // std::cout << "Eigenvalues(LST): "
    //      << eigen_result.eigenvalues().transpose() << std::endl;
    // std::cout << "Here's a matrix whose columns are eigenvectors of LST \n"
    //      << "corresponding to these eigenvalues:\n"
    //      << eigen_result.eigenvectors() << std::endl;

    newPos -= center.position;  // Vector of old position to new position
    Eigen::Map<Eigen::Vector3d> newPos_e(newPos.data());

    // dot product followed by elementwise-division EQN 4. w is a scale factor.
    auto w = (
        (eigen_result.eigenvectors().transpose()*newPos_e).array()
        / (eigen_result.eigenvalues().array()+1)
        ).matrix();                           // vector 3x1
    newPos_e = eigen_result.eigenvectors()*w; // matrix 3x3 * vector = vector
    // center.position += newPos;
    return newPos;
}


/**
 * @brief      Perona-Malik normal based smoothing algorithm
 *
 * @param      mesh      The mesh
 * @param[in]  vertexID  The vertex id
 */
void normalSmoothH(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID, double k)
{
    auto            name = mesh.get_name(vertexID)[0];
    double          areaSum = 0;

    auto            p = (*vertexID).position;
    Eigen::Vector4d pos_e;
    pos_e << p[0], p[1], p[2], 1;   // 4D for affine3D
    Eigen::Vector4d newPos_e;
    newPos_e << 0, 0, 0, 0;

    // For each incident face get the average normal
    auto incidentFaces = mesh.up(mesh.up(vertexID));
    for (auto faceID : incidentFaces)
    {
        // Get the area of incident face
        double area = getArea(mesh, faceID);
        areaSum += area;

        auto normal = getNormal(mesh, faceID);
        normal /= std::sqrt(normal|normal);

        // get the incident incident faces
        std::vector<SurfaceMesh::SimplexID<3> > faces;
        neighbors(mesh, faceID, std::back_inserter(faces));
        Vector avgNorm;

        double sumWeight = 0;

        for (auto face : faces)
        {
            auto   inorm = getNormal(mesh, face);
            inorm /= std::sqrt(inorm|inorm);
            double weight = exp(k*(normal|inorm));
            avgNorm   += weight*inorm;
            sumWeight += weight;
        }
        avgNorm /= sumWeight;

        // Compute the edge (axis) to rotate about.
        auto edge = mesh.get_simplex_down(faceID, name);
        auto edgeName = mesh.get_name(edge);
        auto a = *mesh.get_simplex_up({edgeName[0]});
        auto b = *mesh.get_simplex_up({edgeName[1]});

        auto ab = a-b;
        ab /= std::sqrt(ab|ab); // Eigen AngleAxis requires unit vector

        // Angle between normals in radians. This is the angle to rotate the
        // normal by.
        double angle = std::copysign(std::acos(normal|avgNorm), dot(cross(normal, avgNorm), ab));
        // Catch for floating point dot product issues
        if (std::isnan(angle))
        {
            angle = 0;
        }
        //Vector rotAxis = (*faceID).orientation * ab; // We don't need this
        // because the angle is relative.
        Vector rotAxis = ab;

        // build the transformation
        Eigen::Map<Eigen::Vector3d> rotAxis_e(rotAxis.data());
        Eigen::Map<Eigen::Vector3d> center_e(b.position.data());
        Eigen::Affine3d             A = Eigen::Translation3d(center_e) * Eigen::AngleAxisd(angle, rotAxis_e) * Eigen::Translation3d(-center_e);

        // Weight the new position by the area of the current face
        newPos_e += area*(A*pos_e);
    }
    newPos_e /= areaSum; // weighted by area

    (*vertexID).position = Vertex({newPos_e[0], newPos_e[1], newPos_e[2]});
}

void edgeFlip(SurfaceMesh &mesh, SurfaceMesh::SimplexID<2> edgeID)
{
    // Assuming that the mesh is manifold and edge has been vetted for flipping
    auto name = mesh.get_name(edgeID);
    auto up   = mesh.get_cover(edgeID);

    std::array<SurfaceMesh::SimplexID<3>, 2> faces;
    mesh.up(edgeID, faces.begin());

    // Assume flipped edges are always selected
    auto fdata = SMFace(0, true);
    if ((*faces[0]).marker == (*faces[1]).marker)
    {
        fdata.marker = (*faces[0]).marker;
    }

    mesh.remove<2>({name[0], name[1]});
    mesh.insert<3>({name[0], up[0], up[1]}, fdata);
    mesh.insert<3>({name[1], up[0], up[1]}, fdata);

    edgeID = mesh.get_simplex_up({up[0], up[1]});
    (*edgeID).selected = true;


    std::vector<SurfaceMesh::SimplexID<1>> verts;
    mesh.down(edgeID, std::back_inserter(verts));
    initLocalOrientation<std::integral_constant<size_t, 1> >::apply(mesh, std::set<int>({up[0], up[1], name[0], name[1]}), verts.begin(), verts.end());

    std::vector<SurfaceMesh::SimplexID<2>> nbors;
    casc::neighbors(mesh, edgeID, std::back_inserter(nbors));
    nbors.push_back(edgeID);
    computeLocalOrientation(mesh, nbors);
}

void edgeFlipCache(SurfaceMesh &mesh, SurfaceMesh::SimplexID<2> edgeID)
{
    // Assuming that the mesh is manifold and edge has been vetted for flipping
    auto name = mesh.get_name(edgeID);
    auto up   = mesh.get_cover(edgeID);

    std::array<SurfaceMesh::SimplexID<3>, 2> faces;
    mesh.up(edgeID, faces.begin());

    // Assume flipped edges are always selected
    auto fdata = SMFace(0, true);
    if ((*faces[0]).marker == (*faces[1]).marker)
    {
        fdata.marker = (*faces[0]).marker;
    }

    std::array<SurfaceMesh::SimplexID<3>, 2> newFaces;
    mesh.remove<2>({name[0], name[1]});
    newFaces[0] = mesh.insert<3>({name[0], up[0], up[1]}, fdata);
    newFaces[1] = mesh.insert<3>({name[1], up[0], up[1]}, fdata);

    edgeID = mesh.get_simplex_up({up[0], up[1]});
    (*edgeID).selected = true;

    std::vector<SurfaceMesh::SimplexID<1>> verts;
    mesh.down(edgeID, std::back_inserter(verts));
    initLocalOrientation<std::integral_constant<size_t, 1> >::apply(mesh, std::set<int>({up[0], up[1], name[0], name[1]}), verts.begin(), verts.end());

    std::vector<SurfaceMesh::SimplexID<2>> nbors;
    casc::neighbors(mesh, edgeID, std::back_inserter(nbors));
    nbors.push_back(edgeID);
    computeLocalOrientation(mesh, nbors);

    for(auto fID : newFaces){
        auto norm = getNormal(mesh, fID);
        normalize(norm);
        (*fID).normal = norm;
    }

    for(auto vID : verts){
        Vector norm;
        auto   faces = mesh.up(mesh.up(vID));
        for (auto faceID : faces)
        {
            norm += (*faceID).normal;
        }
        normalize(norm);
        (*vID).normal = norm;
    }
}

bool checkFlipAngle(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID)
{
    auto getMinAngle = [](const Vertex &a, const Vertex &b, const Vertex &c){
            double              minAngle = 999; // dummy for now
            double              tmp;
            std::array<const Vertex, 3> triangle = {a, b, c};
            std::array<std::size_t, 3> sigma = {0,1,2};
            for(int i = 0; i < 3; ++i){
                try
                {
                    tmp = angleDeg(triangle[sigma[0]], triangle[sigma[1]], triangle[sigma[2]]);
                }
                catch (std::runtime_error &e)
                {
                    throw std::runtime_error("Angle is undefined for face with zero area. Try running degenerate dissolve in Blender and ensure manifoldness.");
                }
                if (tmp < minAngle)
                minAngle = tmp;

                std::rotate(sigma.begin(), sigma.begin() + 1, sigma.end());
            }
            return minAngle;
        };

    auto name = mesh.get_name(edgeID);
    std::pair<Vertex, Vertex> shared;

    shared.first  = *mesh.get_simplex_down(edgeID, name[0]);
    shared.second = *mesh.get_simplex_down(edgeID, name[1]);

    // TODO: (50) Benchmark for performance increase from going up to face then back down...
    std::pair<Vertex, Vertex> notShared;
    auto up   = mesh.get_cover(edgeID);
    notShared.first  = *mesh.get_simplex_up({up[0]});
    notShared.second = *mesh.get_simplex_up({up[1]});

    // TODO: (5) Change to check link condition
    // Check if the triangles make a wedge shape. Flipping this can
    // cause knife faces.
    auto   f1   = (shared.first - notShared.first)^(shared.second - notShared.first);
    auto   f2   = (shared.first - notShared.second)^(shared.second - notShared.second);
    auto   f3   = (notShared.first - shared.first)^(notShared.second - shared.first);
    auto   f4   = (notShared.first - shared.second)^(notShared.second - shared.second);
    auto   area = std::pow(std::sqrt(f1|f1) + std::sqrt(f2|f2), 2);
    auto   areaFlip = std::pow(std::sqrt(f3|f3) + std::sqrt(f4|f4), 2);
    double edgeFlipCriterion = 1.001;
    // If area changes by a lot continue
    if (areaFlip/area > edgeFlipCriterion)
        return false;

    int excess = checkFlipValenceExcess(mesh, edgeID);
    if (excess < -3){
        // Force flip since it improves valence by a lot
        return true;
    }
    else if (excess > 3){
        // Flipping is too bad...
        return false;
    }

    // Go through all angle combinations
    double tmp;
    double minAngle = getMinAngle(shared.first, shared.second, notShared.first);
    tmp = getMinAngle(shared.first, shared.second, notShared.second);
    if (tmp < minAngle)
        minAngle = tmp;

    double minAngleFlip = getMinAngle(notShared.first, notShared.second, shared.first);
    tmp = getMinAngle(notShared.first, notShared.second, shared.second);
    if (tmp < minAngleFlip)
        minAngleFlip = tmp;

    if (minAngleFlip > minAngle){
        return true;
    }
    return false;
}

int checkFlipValenceExcess(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID)
{
    auto name = mesh.get_name(edgeID);
    std::pair<SurfaceMesh::SimplexID<1>, SurfaceMesh::SimplexID<1> > shared;
    shared.first  = mesh.get_simplex_down(edgeID, name[0]);
    shared.second = mesh.get_simplex_down(edgeID, name[1]);
    std::pair<SurfaceMesh::SimplexID<1>, SurfaceMesh::SimplexID<1> > notShared;
    auto up = mesh.get_cover(edgeID);
    notShared.first  = mesh.get_simplex_up({up[0]});
    notShared.second = mesh.get_simplex_up({up[1]});
    std::array<int, 4> valence;

    // assuming there are no boundaries
    valence[0] = getValence(mesh, shared.first)-6;
    valence[1] = getValence(mesh, shared.second)-6;
    valence[2] = getValence(mesh, notShared.first)-6;
    valence[3] = getValence(mesh, notShared.second)-6;

    // std::cout << "Before: " << casc::to_string(valence) << std::endl;
    int excess = 0;
    for (int i = 0; i < 4; i++)
    {
        excess += std::pow(valence[i], 2);
    }
    // simulate the flip
    valence[0] -= 1;
    valence[1] -= 1;
    valence[2] += 1;
    valence[3] += 1;

    // std::cout << "After: " << casc::to_string(valence) << std::endl;
    int flipExcess = 0;
    for (int i = 0; i < 4; i++)
    {
        flipExcess += std::pow(valence[i], 2);
    }
    return flipExcess - excess;
}

// Angle Mesh improvement by projecting barycenter on tangent plane...
void barycenterVertexSmooth(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    // get the neighbors
    std::vector<SurfaceMesh::SimplexID<1> > vertices;
    neighbors_up(mesh, vertexID, std::back_inserter(vertices));
    // compute the average position
    Vector avgPos;
    for (auto vertex : vertices)
    {
        avgPos += (*vertex).position;
    }
    avgPos /= vertices.size();

    auto disp = avgPos - (*vertexID).position;
    // Project onto tangent plane
    // A||B = Bx(AxB/|B|)/|B|
    // A_|_B = A.B*B/|B|^2
    //auto norm = getNormalFromTangent(getTangent(mesh, vertexID));
    auto norm = getNormal(mesh, vertexID);
    norm /= std::sqrt(norm|norm); // normalize
    auto perpProj = norm*norm;    // tensor product

    Eigen::Map<Eigen::Vector3d> disp_e(disp.data());

    // Compute perpendicular component
    // Vector perp;
    // Eigen::Map<Eigen::Matrix3d> perpProj_e(perpProj.data());
    // Eigen::Map<Eigen::Vector3d> perp_e(perp.data());
    // perp_e = perpProj_e*disp_e;

    // Compute the parallel component
    Vector               parallel;
    tensor<double, 3, 2> identity{{
        1, 0, 0, 0, 1, 0, 0, 0, 1
    }};
    auto llproj = identity-perpProj; // perpendicular projector
    Eigen::Map<Eigen::Matrix3d> llproj_e(llproj.data());
    Eigen::Map<Eigen::Vector3d> parallel_e(parallel.data());
    parallel_e = llproj_e*disp_e;

    (*vertexID).position = (*vertexID).position + parallel;
}

void findHoles(const SurfaceMesh                                     &mesh,
               std::vector<std::vector<SurfaceMesh::SimplexID<2> > > &holeList)
{
    // TODO: (25) Update this to return a pair of edge ring and vertex ring
    std::set<SurfaceMesh::SimplexID<2> > bdryEdges;

    // Collect all of the boundary edges into boundarySet
    for (auto edgeID : mesh.get_level_id<2>())
    {
        auto cover = mesh.get_cover(edgeID);
        if (cover.size() == 1)
        {
            bdryEdges.insert(edgeID);
        }
    }

    // Connect the edges into rings...
    while (bdryEdges.size() > 0)
    {
        // Container to store ring
        std::vector<SurfaceMesh::SimplexID<2> > bdryRing;
        // // Pop first edge from remaining edges
        // auto it = bdryEdges.begin();
        // auto firstEdge = *it;
        // bdryRing.push_back(firstEdge);
        // bdryEdges.erase(it);
        std::vector<SurfaceMesh::SimplexID<1> > visitedVerts;
        // mesh.down(firstEdge, std::back_inserter(visitedVerts));
        // Try to complete the ring
        if (!orderBoundaryEdgeRing(mesh, bdryEdges, visitedVerts, bdryRing))
        {
            throw std::runtime_error("Couldn't connect ring");
        }
        holeList.push_back(bdryRing);
    }
}

bool orderBoundaryEdgeRing(const SurfaceMesh                       &mesh,
                           std::set<SurfaceMesh::SimplexID<2> >    &unvisitedBdryEdges,
                           std::vector<SurfaceMesh::SimplexID<1> > &visitedVerts,
                           std::vector<SurfaceMesh::SimplexID<2> > &bdryRing)
{
    if (visitedVerts.empty() && bdryRing.empty())
    {
        // Pop first edge from unvisited edges
        auto it = unvisitedBdryEdges.begin();
        auto firstEdge = *it;
        // Push the edge into the ring
        bdryRing.push_back(firstEdge);
        unvisitedBdryEdges.erase(it);
        // Push vertices of edge into visited vertices
        mesh.down(firstEdge, std::back_inserter(visitedVerts));
    }

    auto firstEdge = bdryRing.front();
    auto currEdge  = bdryRing.back();
    auto currVert  = visitedVerts.back();

    std::vector<SurfaceMesh::SimplexID<2> > nbors;
    neighbors(mesh, currEdge, std::back_inserter(nbors));

    // Look for connected unvisited boundary edge
    for (auto nbor : nbors)
    {
        auto result = unvisitedBdryEdges.find(nbor);
        if (result != unvisitedBdryEdges.end())
        {
            auto nextEdge = *result;
            std::array<SurfaceMesh::SimplexID<1>, 2> verts;
            mesh.down(nextEdge, verts.begin());
            int found = 0;
            for (; found < 2; ++found)
            {
                if (verts[found] == currVert)
                    break;
            }

            SurfaceMesh::SimplexID<1> nextVert;

            if (found == 0)
            {
                nextVert = verts[1];
            }
            else if (found == 1)
            {
                nextVert = verts[0];
            }
            else
            {
                continue;
            }

            // If this closes the ring...
            if (nextVert == visitedVerts.front())
            {
                bdryRing.push_back(nextEdge);
                unvisitedBdryEdges.erase(result);
                return true;
            }

            // Check if we have visited the next vertex already
            auto visited = std::find(visitedVerts.begin(), visitedVerts.end(), nextVert);
            if (visited == visitedVerts.end())
            {
                visitedVerts.push_back(nextVert);
                bdryRing.push_back(nextEdge);
                unvisitedBdryEdges.erase(result);

                if (orderBoundaryEdgeRing(mesh, unvisitedBdryEdges, visitedVerts, bdryRing))
                {
                    return true;
                }
            }
        }
    }

    // Ring isn't complete. Backtrack...
    unvisitedBdryEdges.insert(currEdge);
    visitedVerts.pop_back();
    bdryRing.pop_back();
    return false;
}

void edgeRingToVertices(const SurfaceMesh                                                  &mesh,
                        std::vector<SurfaceMesh::SimplexID<2> >                            &edgeRing,
                        std::back_insert_iterator<std::vector<SurfaceMesh::SimplexID<1> > > iter)
{
    auto edgeRingIT = edgeRing.begin();
    std::array<SurfaceMesh::SimplexID<1>, 2> first, next;
    SurfaceMesh::SimplexID<1>                prevVertex;

    mesh.down(*edgeRingIT++, first.begin());
    mesh.down(*edgeRingIT, next.begin());

    // Set the first edge
    if (first[0] != next[0] && first[0] != next[1])
    {
        *iter++ = first[0];
        *iter++ = first[1];
        prevVertex = first[1];
    }
    else
    {
        *iter++ = first[1];
        *iter++ = first[0];
        prevVertex = first[0];
    }

    // Go through remaining Edges
    for (; edgeRingIT != edgeRing.end()-1; ++edgeRingIT)
    {
        mesh.down(*edgeRingIT, next.begin());
        if (next[0] == prevVertex)
        {
            *iter++ = next[1];
            prevVertex = next[1];
        }
        else
        {
            *iter++ = next[0];
            prevVertex = next[0];
        }
    }
}

bool checkEdgeFlip(const SurfaceMesh &mesh,
                     bool preserveRidges,
                     SurfaceMesh::SimplexID<2> edgeID,
                     std::function<bool(const SurfaceMesh &, const SurfaceMesh::SimplexID<2> &)> &&checkFlip
                )
{
    if ((*edgeID).selected == false)
        return false;

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
        return false;
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
        return false;
    }

    // Check if we're on a ridge. This prevents folding also.
    if (preserveRidges)
    {
        auto a   = getNormal(mesh, mesh.get_simplex_up(edgeID, up[0]));
        auto b   = getNormal(mesh, mesh.get_simplex_up(edgeID, up[1]));
        auto val = angleDeg(a, b);
        if (val > 60)
        {
            return false;
        }
    }

    // TODO: (5) Implement a better algorithm to prevent folding
    // Check if the triangles make a wedge shape. Flipping this can
    // cause knife faces.
    auto   f1   = (shared.first - notShared.first)^(shared.second - notShared.first);
    auto   f2   = (shared.first - notShared.second)^(shared.second - notShared.second);
    auto   f3   = (notShared.first - shared.first)^(notShared.second - shared.first);
    auto   f4   = (notShared.first - shared.second)^(notShared.second - shared.second);
    auto   area = std::pow(std::sqrt(f1|f1) + std::sqrt(f2|f2), 2);
    auto   areaFlip = std::pow(std::sqrt(f3|f3) + std::sqrt(f4|f4), 2);
    double edgeFlipCriterion = 1.001;
    // If area changes by a lot continue
    if (areaFlip/area > edgeFlipCriterion)
        return false;

    // Check the flip using user function
    return checkFlip(mesh, edgeID);
}
} // end namespace surfacemesh_detail
} // end namespace gamer
