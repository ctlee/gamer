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

#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <strstream>
#include <vector>
#include <casc/casc>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gamer/EigenDiagonalization.h"
#include "gamer/SurfaceMesh.h"
#include "gamer/Vertex.h"

/// Namespace for all things gamer
namespace gamer
{
void print(const SurfaceMesh& mesh)
{
    std::cout << "Level: 1" << std::endl;
    for (auto node : mesh.get_level_id<1>())
    {
        std::cout << "    " << *node << std::endl;
    }
    std::cout << "Level: 2" << std::endl;
    for (auto node : mesh.get_level_id<2>())
    {
        auto name = mesh.get_name(node);
        std::cout << "    " << casc::to_string(name) << std::endl;
    }
    std::cout << "Level: 3" << std::endl;
    for (auto node : mesh.get_level_id<3>())
    {
        auto name = mesh.get_name(node);
        std::cout << "    " << casc::to_string(name) << std::endl;
    }
}

void generateHistogram(const SurfaceMesh& mesh)
{
    // compute angle distribution
    std::array<double, 18> histogram;
    histogram.fill(0);
    for (auto face : mesh.get_level_id<3>())
    {
        auto vertexIDs = mesh.get_name(face);
        // Unpack the ID's for convenience
        Vertex a = *mesh.get_simplex_up<1>({vertexIDs[0]});
        Vertex b = *mesh.get_simplex_up<1>({vertexIDs[1]});
        Vertex c = *mesh.get_simplex_up<1>({vertexIDs[2]});

        auto binAngle = [&](double angle) -> std::size_t {
                            return static_cast<std::size_t>(std::floor(angle/10));
                        };
        histogram[binAngle(angleDeg(a, b, c))]++;
        histogram[binAngle(angleDeg(b, a, c))]++;
        histogram[binAngle(angleDeg(c, a, b))]++;
    }

    std::size_t factor = mesh.size<3>()*3;
    std::for_each(histogram.begin(), histogram.end(), [&factor](double& n){
            n = 100.0*n/factor;
        });

    std::cout << "Angle Distribution:" << std::endl;
    for (std::size_t x = 0; x < 18; x++)
        std::cout << x*10 << "-" << (x+1)*10 << ": " << std::setprecision(2)
                  << std::fixed << histogram[x] << std::endl;
    std::cout << std::endl << std::endl;

    // compute the edge length distribution
    std::cout << "Edge Length Distribution:" << std::endl;
    std::vector<double> lengths;
    for (auto edge : mesh.get_level_id<2>())
    {
        auto vertexIDs = mesh.down(edge);
        auto t1  = *vertexIDs.cbegin();
        auto t2  = *(++vertexIDs.cbegin());
        auto v1  = *t1;
        auto v2  = *t2;
        double len = length(v2-v1);
        lengths.push_back(len);
    }
    std::sort(lengths.begin(), lengths.end());

    std::array<double, 20> histogramLength;
    histogramLength.fill(0);
    double interval = (lengths.back() - lengths.front())/20;
    double low = lengths.front();

    if (interval <= 0.0000001) // floating point roundoff prevention
    {
        std::cout << lengths.front() << ": " << 100 << std::endl << std::endl;
    }
    else
    {
        for (auto length : lengths)
        {
            histogramLength[std::floor((length-low)/interval)]++;
        }

        factor = mesh.size<2>();
        std::for_each(histogramLength.begin(), histogramLength.end(), [&factor](double& n){
                n = 100.0*n/factor;
            });

        for (std::size_t x = 0; x < 20; x++)
            std::cout << x*interval << "-" << (x+1)*interval << ": " << std::setprecision(2)
                      << std::fixed << histogramLength[x] << std::endl;
        std::cout << std::endl << std::endl;
    }

    // Compute the valence distribution
    std::array<double, 20> histogramValence;
    histogramValence.fill(0);

    for (auto vertexID : mesh.get_level_id<1>())
    {
        // TODO bounds checking here...
        histogramValence[getValence(mesh, vertexID)]++;
    }

    factor = mesh.size<1>();
    // std::for_each(histogramValence.begin(), histogramValence.end(),
    // [&factor](double& n){
    //         n = 100.0*n/factor;});
    std::cout << "Valence distribution:" << std::endl;
    for (int x = 0; x < 20; x++)
        std::cout << x << ": " << histogramValence[x] << std::endl;
    std::cout << std::endl << std::endl;
}

void printQualityInfo(const std::string& filename, const SurfaceMesh& mesh)
{

    std::stringstream anglefilename;
    anglefilename << filename << ".angle";

    std::ofstream angleOut(anglefilename.str());
    if (!angleOut.is_open())
    {
        std::stringstream ss;
        ss << "Couldn't open file " << filename << ".angle";
        throw std::runtime_error(ss.str());
    }

    std::stringstream areafilename;
    areafilename << filename << ".area";
    std::ofstream areaOut(areafilename.str());
    if (!areaOut.is_open())
    {
        std::stringstream ss;
        ss << "Couldn't open file " << filename << ".area";
        throw std::runtime_error(ss.str());
    }

    for (auto faceID : mesh.get_level_id<3>())
    {
        auto name = mesh.get_name(faceID);
        auto a = *mesh.get_simplex_up({name[0]});
        auto b = *mesh.get_simplex_up({name[1]});
        auto c = *mesh.get_simplex_up({name[2]});

        areaOut << getArea(a, b, c) << "\n";
        // Print angles in degrees
        angleOut << angleDeg(a, b, c) << "\n"
                 << angleDeg(b, a, c) << "\n"
                 << angleDeg(a, c, b) << "\n";
    }
    areaOut.close();
    angleOut.close();
}

std::tuple<double, double, int, int> getMinMaxAngles(
    const SurfaceMesh& mesh,
    double maxMinAngle,
    double minMaxAngle)
{
    double minAngle = 360;
    double maxAngle = 0;
    int small = 0;
    int large = 0;

    // for each triangle
    for (auto nid : mesh.get_level_id<3>())
    {
        auto name = mesh.get_name(nid);
        auto a = *mesh.get_simplex_up({name[0]});
        auto b = *mesh.get_simplex_up({name[1]});
        auto c = *mesh.get_simplex_up({name[2]});
        std::array<double, 3> angles;
        try
        {
            angles[0] = angleDeg(a, b, c);
            angles[1] = angleDeg(b, a, c);
            angles[2] = angleDeg(a, c, b);
        }
        catch (std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
            throw std::runtime_error("ERROR(getMinMaxAngles): Cannot compute angles of face with zero area.");
        }

        for (double angle : angles)
        {
            if (angle < minAngle)
            {
                minAngle = angle;
            }
            if (angle > maxAngle)
            {
                maxAngle = angle;
            }
            if (angle < maxMinAngle)
            {
                ++small;
            }
            if (angle > minMaxAngle)
            {
                ++large;
            }
        }
    }
    return std::make_tuple(minAngle, maxAngle, small, large);
}

double getArea(const SurfaceMesh& mesh)
{
    double area = 0.0;
    for (auto faceID : mesh.get_level_id<3>())
        area += getArea(mesh, faceID);
    return area;
}

double getArea(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<3> faceID)
{
    auto name = mesh.get_name(faceID);
    auto a = *mesh.get_simplex_up({name[0]});
    auto b = *mesh.get_simplex_up({name[1]});
    auto c = *mesh.get_simplex_up({name[2]});
    return getArea(a, b, c);
}

double getArea(Vertex a, Vertex b, Vertex c)
{
    auto wedge = (b-a)^(b-c);
    return std::sqrt(wedge|wedge)/2;
}

/**
 * http://research.microsoft.com/en-us/um/people/chazhang/publications/icip01_ChaZhang.pdf
 * http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
 */
double getVolume(const SurfaceMesh& mesh)
{
    bool orientError = false;
    double volume = 0;
    for (auto faceID : mesh.get_level_id<3>())
    {
        auto name = mesh.get_name(faceID);
        auto a = (*mesh.get_simplex_up({name[0]})).position;
        auto b = (*mesh.get_simplex_up({name[1]})).position;
        auto c = (*mesh.get_simplex_up({name[2]})).position;

        Vector norm;
        double tmp = 0.0;
        if ((*faceID).orientation == 1)
        {
            // a->b->c
            norm = cross(b, c);
            tmp  = dot(a, norm);

        }
        else if ((*faceID).orientation == -1)
        {
            // c->b->a
            norm = cross(b, a);
            tmp  = dot(c, norm);
        }
        else
        {
            orientError = true;
        }

        // * Probably less efficient way #1
        // tmp = dot(a, getNormal(mesh, faceID));

        // * Less efficient way #2
        // norm = getNormalFromTangent(getTangent(mesh, faceID));
        // auto wedge = a^b^c;
        // tmp = std::sqrt(wedge|wedge);
        // auto sgn = dot((a+b+c)/3, norm);
        // if(sgn <= 0) {
        //     tmp = -1*tmp;
        // }
        volume += tmp;
    }
    if (orientError)
    {
        std::cerr << "ERROR getVolume(): Orientation undefined for one or more "
                  << "simplices. Did you call compute_orientation()?" << std::endl;
    }
    return volume/6;
}

void translate(SurfaceMesh& mesh, Vector v)
{
    for (auto& vertex : mesh.get_level<1>())
        vertex += v;
}

void translate(SurfaceMesh& mesh, double dx, double dy, double dz)
{
    Vector v = Vector({dx, dy, dz});
    translate(mesh, v);
}

void scale(SurfaceMesh& mesh, Vector v)
{
    for (auto& vertex : mesh.get_level<1>())
    {
        vertex[0] *= v[0];
        vertex[1] *= v[1];
        vertex[2] *= v[2];
    }
}

void scale(SurfaceMesh& mesh, double sx, double sy, double sz)
{
    Vector v = Vector();
    v[0] = sx; v[1] = sy; v[2] = sz;
    scale(mesh, v);
}

void scale(SurfaceMesh& mesh, double s)
{
    Vector v = Vector();
    v[0] = s; v[1] = s; v[2] = s;
    scale(mesh, v);
}

std::pair<Vector, double> getCenterRadius(SurfaceMesh& mesh)
{
    Vector center = Vector();
    double radius = 0;
    if (mesh.size<1>() > 0)
    {
        for (auto vertexID : mesh.get_level_id<1>())
        {
            center += *vertexID;
        }
        center /= static_cast<double>(mesh.size<1>());

        for (auto vertexID : mesh.get_level_id<1>())
        {
            Vector tmp  = *vertexID - center;
            double dist = std::sqrt(tmp|tmp);
            if (dist > radius)
                radius = dist;
        }
    }
    return std::make_pair(center, radius);
}

void centeralize(SurfaceMesh& mesh)
{
    Vector center;
    double radius;
    std::tie(center, radius) = getCenterRadius(mesh);
    translate(mesh, -center);
}

void normalSmooth(SurfaceMesh& mesh)
{
    for (auto nid : mesh.get_level_id<1>())
    {
        surfacemesh_detail::normalSmoothH(mesh, nid, 2);
    }
    double min, max;
    int nSmall, nLarge;
    std::tie(min, max, nSmall, nLarge) = getMinMaxAngles(mesh, 15, 150);
    std::cout << "  Min Angle: " << min << ", Max Angle: " << max
              << ", # Small Angles: " << nSmall
              << ", # Large Angles: " << nLarge << std::endl;
}

tensor<double, 3, 2> getTangent(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    return surfacemesh_detail::getTangentH(mesh, (*vertexID).position, vertexID);
}

tensor<double, 3, 2> getTangent(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<3> faceID)
{
    auto cover = mesh.get_name(faceID);
    auto vertexID = mesh.get_simplex_up({cover[0]});
    std::set<SurfaceMesh::KeyType> next(std::begin(cover)+1, std::end(cover));
    return surfacemesh_detail::getTangentF(mesh, (*vertexID).position, vertexID, next);
}

Vector getNormalFromTangent(const tensor<double, 3, 2> tangent)
{
    Vector xp;
    // upper right...
    xp[0] = -tangent.get(1, 2);  // x
    xp[1] = tangent.get(0, 2);   // y
    xp[2] = -tangent.get(0, 1);  // z
    // lower left...
    // xp[0] = tangent.get(2,1);
    // xp[1] = -tangent.get(2,0);
    // xp[2] = tangent.get(1,0);
    return xp;
}

Vector getNormal(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    Vector norm;
    auto faces = mesh.up(mesh.up(vertexID));
    for (auto faceID : faces)
    {
        norm += getNormal(mesh, faceID);
    }
    // norm /= faces.size();
    return norm;
}

Vector getNormal(const SurfaceMesh& mesh, SurfaceMesh::SimplexID<3> faceID)
{
    Vector norm;
    auto name = mesh.get_name(faceID);

    // Get the three vertices
    auto a = *mesh.get_simplex_up({name[0]});
    auto b = *mesh.get_simplex_up({name[1]});
    auto c = *mesh.get_simplex_up({name[2]});

    if ((*faceID).orientation == 1)
    {
        norm = cross(c-b, a-b);
    }
    else if ((*faceID).orientation == -1)
    {
        norm = cross(a-b, c-b);
    }
    else
    {
        // std::cerr << "ERROR(getNormal): Orientation undefined, cannot compute
        // "
        // << "normal. Did you call compute_orientation()?" << std::endl;
        throw std::runtime_error("ERROR(getNormal): Orientation undefined, cannot compute normal. Did you call compute_orientation()?");
    }
    return norm;
}


void smoothMesh(SurfaceMesh& mesh, int maxIter, bool preserveRidges, std::size_t rings, bool verbose)
{
    double maxMinAngle = 15;
    double minMaxAngle = 165;
    double minAngle, maxAngle;
    int nSmall, nLarge;
    int nIter = 1;

    if (verbose)
    {
        std::tie(minAngle, maxAngle, nSmall, nLarge) = getMinMaxAngles(mesh, maxMinAngle, minMaxAngle);
        std::cout << "Initial Quality: Min Angle = " << minAngle << ", "
                  << "Max Angle = " << maxAngle << ", "
                  << "# smaller-than-" << maxMinAngle << " = " << nSmall << ", "
                  << "# larger-than-" << minMaxAngle << " = " << nLarge << std::endl;
    }

    std::vector<std::pair<SurfaceMesh::SimplexID<1>, Vector>> delta;

    // Cache normals before entering loop
    cacheNormals(mesh);
    for (int nIter = 1; nIter <= maxIter; ++nIter)
    {
        for (auto vertex : mesh.get_level_id<1>())
        {
            if ((*vertex).selected == true)
            {
                // surfacemesh_detail::weightedVertexSmooth(mesh, vertex,
                // rings);
                delta.push_back(std::make_pair(vertex, surfacemesh_detail::weightedVertexSmoothCache(mesh, vertex, rings)));
            }
            //barycenterVertexSmooth(mesh, vertex);
        }

        for (auto pair : delta) {
            *pair.first += pair.second;
        }
        delta.clear();
        cacheNormals(mesh);

        // ATOMIC EDGE FLIP
        std::vector<SurfaceMesh::SimplexID<2>> edgesToFlip;
        // Get set of good, non-interfering edges to flip according to the
        // Angle based criteria.
        surfacemesh_detail::selectFlipEdges(mesh,
                                            preserveRidges,
                                            surfacemesh_detail::checkFlipAngle,
                                            std::back_inserter(edgesToFlip));
        for (auto edgeID : edgesToFlip)
        {
            surfacemesh_detail::edgeFlipCache(mesh, edgeID);
        }

        if (verbose)
        {
            std::tie(minAngle, maxAngle, nSmall, nLarge) = getMinMaxAngles(mesh, maxMinAngle, minMaxAngle);
            std::cout << "Iteration " << nIter << ":" << std::endl;
            std::cout << "Min Angle = " << minAngle << ", "
                      << "Max Angle = " << maxAngle << ", "
                      << "# smaller-than-" << maxMinAngle << " = " << nSmall << ", "
                      << "# larger-than-" << minMaxAngle << " = " << nLarge << std::endl;
        }
    }
}

/**
 * @brief      Refine the mesh by quadrisection.
 *
 * Note that this function will delete all stored data on edges and faces. But
 * this can be easily fixed.
 *
 * @param      mesh  The mesh
 */
std::unique_ptr<SurfaceMesh> refineMesh(const SurfaceMesh& mesh)
{
    std::unique_ptr<SurfaceMesh> refinedMesh(new SurfaceMesh);

    // Copy over vertices to refinedMesh
    for (auto vertex : mesh.get_level_id<1>())
    {
        auto key = mesh.get_name(vertex);
        refinedMesh->insert(key, *vertex);
    }

    // Split edges and generate a map of names before to after
    std::map<std::array<int, 2>, int> edgeMap;

    for (auto edge : mesh.get_level_id<2>())
    {
        auto edgeName = mesh.get_name(edge);
        Vector v1 = (*mesh.get_simplex_up({edgeName[0]})).position;
        Vector v2 = (*mesh.get_simplex_up({edgeName[1]})).position;

        auto newVertex = refinedMesh->add_vertex(SMVertex(std::move(0.5*(v1+v2))));
        edgeMap.emplace(std::make_pair(edgeName, newVertex));
    }

    // Connect edges and face and copy data
    for (auto face : mesh.get_level_id<3>())
    {
        auto name = mesh.get_name(face);
        int a, b, c;

        // Skip checking if found
        auto it = edgeMap.find({name[0], name[1]});
        a = it->second;

        it = edgeMap.find({name[1], name[2]});
        b  = it->second;

        it = edgeMap.find({name[0], name[2]});
        c  = it->second;

        refinedMesh->insert({name[0], a});
        refinedMesh->insert({a, name[1]});

        refinedMesh->insert({name[1], b});
        refinedMesh->insert({b, name[2]});

        refinedMesh->insert({name[0], c});
        refinedMesh->insert({c, name[2]});

        refinedMesh->insert({a, b, c});
        refinedMesh->insert({name[0], a, c}, *face);
        refinedMesh->insert({name[1], a, b}, *face);
        refinedMesh->insert({name[2], b, c}, *face);
    }
    return refinedMesh;
}


void coarse(SurfaceMesh& mesh, double coarseRate, double flatRate, double denseWeight, std::size_t rings, bool verbose)
{
    // TODO: Check if all polygons are closed (0)

    // Compute the average edge length
    double avgLen = 0;
    if (denseWeight > 0)
    {
        for (auto edgeID : mesh.get_level_id<2>())
        {
            auto name =  mesh.get_name(edgeID);
            auto v = *mesh.get_simplex_down(edgeID, name[0])
                     - *mesh.get_simplex_down(edgeID, name[1]);
            avgLen += length(v);
        }
        avgLen /= static_cast<double>(mesh.size<2>());
    }

    double sparsenessRatio = 1;
    double flatnessRatio   = 1;

    auto range = mesh.get_level_id<1>();

    for (auto vertexIDIT = range.begin(); vertexIDIT != range.end();)
    {
        // Immediately cache vertexID and increment IT so destruction of
        // vertexID won't invalidate the iterator.
        auto vertexID = *vertexIDIT;
        ++vertexIDIT;

        if ((*vertexID).selected == false)
        {
            continue;
        }

        // Sparseness as coarsening criteria
        if (denseWeight > 0)
        {
            // Get max length of edges.
            auto edges  = mesh.up(vertexID);
            double maxLen = 0;
            for (auto edgeID : edges)
            {
                auto name =  mesh.get_name(edgeID);
                auto v = *mesh.get_simplex_down(edgeID, name[0])
                         - *mesh.get_simplex_down(edgeID, name[1]);
                double tmpLen = std::sqrt(v|v);
                if (tmpLen > maxLen)
                    maxLen = tmpLen;
            }
            sparsenessRatio = std::pow(maxLen/avgLen, denseWeight);
        }

        // Curvature as coarsening criteria
        if (flatRate > 0)
        {
            auto lst = surfacemesh_detail::computeLocalStructureTensor(mesh, vertexID, rings);

            EigenVector eigenvalues;
            EigenMatrix eigenvectors;

            EigenDiagonalizeTraits<REAL, 3>::diagonalizeSelfAdjointMatrix(lst, eigenvalues, eigenvectors);

            // The closer this ratio is to 0 the flatter the local region.
            flatnessRatio = std::pow(eigenvalues[1]/eigenvalues[2], flatRate);
        }

        // Add vertex to delete list
        if (sparsenessRatio * flatnessRatio < coarseRate)
        {
            surfacemesh_detail::decimateVertex(mesh, vertexID);
        }
    }
}

void coarse_dense(SurfaceMesh& mesh, REAL threshold, REAL weight, std::size_t rings, bool verbose)
{
    // Compute the average edge length
    REAL avgLen = 0;
    for (auto edgeID : mesh.get_level_id<2>())
    {
        auto name =  mesh.get_name(edgeID);
        auto v = *mesh.get_simplex_down(edgeID, name[0])
                 - *mesh.get_simplex_down(edgeID, name[1]);
        avgLen += length(v);
    }
    avgLen /= static_cast<REAL>(mesh.size<2>());

    REAL sparsenessRatio = 1;

    auto range = mesh.get_level_id<1>();
    for (auto vertexIDIT = range.begin(); vertexIDIT != range.end();)
    {
        // Immediately cache vertexID and increment IT so destruction of
        // vertexID won't invalidate the iterator.
        auto vertexID = *vertexIDIT;
        ++vertexIDIT;

        if ((*vertexID).selected == false)
        {
            continue;
        }

        // Get max length of edges.
        auto edges  = mesh.up(vertexID);
        REAL maxLen = 0;
        for (auto edgeID : edges)
        {
            auto name =  mesh.get_name(edgeID);
            auto v = *mesh.get_simplex_down(edgeID, name[0])
                     - *mesh.get_simplex_down(edgeID, name[1]);
            REAL tmpLen = std::sqrt(v|v);
            if (tmpLen > maxLen)
                maxLen = tmpLen;
        }
        sparsenessRatio = std::pow(maxLen/avgLen, weight);

        // Decimate if under the threshold
        if (sparsenessRatio < threshold)
        {
            surfacemesh_detail::decimateVertex(mesh, vertexID, rings);
        }
    }
}

void coarse_flat(SurfaceMesh& mesh, REAL threshold, REAL weight, std::size_t rings, bool verbose){
    REAL flatnessRatio = 1;

    auto range = mesh.get_level_id<1>();
    for (auto vertexIDIT = range.begin(); vertexIDIT != range.end();)
    {
        // Immediately cache vertexID and increment IT so destruction of
        // vertexID won't invalidate the iterator.
        auto vertexID = *vertexIDIT;
        ++vertexIDIT;

        if ((*vertexID).selected == false)
        {
            continue;
        }

        // Curvature as coarsening criteria
        auto lst = surfacemesh_detail::computeLocalStructureTensor(mesh, vertexID, rings);

        EigenVector eigenvalues;
        EigenMatrix eigenvectors;

        EigenDiagonalizeTraits<REAL, 3>::diagonalizeSelfAdjointMatrix(lst, eigenvalues, eigenvectors);

        // The closer this ratio is to 0 the flatter the local region.
        flatnessRatio = std::pow(eigenvalues[1]/eigenvalues[2], weight);

        // Add vertex to delete list
        if (flatnessRatio < threshold)
        {
            surfacemesh_detail::decimateVertex(mesh, vertexID);
        }
    }
}


void fillHoles(SurfaceMesh& mesh)
{
    std::vector<std::vector<SurfaceMesh::SimplexID<2>>> holeList;
    surfacemesh_detail::findHoles(mesh, holeList);

    for (auto& holeEdges : holeList)
    {
        std::vector<SurfaceMesh::SimplexID<1>> sortedVertices;
        surfacemesh_detail::edgeRingToVertices(mesh, holeEdges, std::back_inserter(sortedVertices));

        surfacemesh_detail::triangulateHole(mesh, sortedVertices, SMFace(), holeEdges);
    }
}

void flipNormals(SurfaceMesh& mesh)
{
    for (auto& fdata : mesh.get_level<3>())
    {
        fdata.orientation *= -1;
    }
}

bool hasHole(const SurfaceMesh& mesh)
{
    for (auto eID : mesh.get_level_id<2>())
    {
        auto cover = mesh.get_cover(eID);
        if (cover.size() != 2)
        {
            return true;
        }
    }
    return false;
}

std::unique_ptr<SurfaceMesh> sphere(int order)
{
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    mesh->insert({0}, SMVertex(0, 0, -1));
    mesh->insert({1}, SMVertex(0.723607, -0.525725, -0.44722));
    mesh->insert({2}, SMVertex(-0.276388, -0.850649, -0.44722));
    mesh->insert({3}, SMVertex(-0.894426, 0, -0.447216));
    mesh->insert({4}, SMVertex(-0.276388, 0.850649, -0.44722));
    mesh->insert({5}, SMVertex(0.723607, 0.525725, -0.44722));
    mesh->insert({6}, SMVertex(0.276388, -0.850649, 0.44722));
    mesh->insert({7}, SMVertex(-0.723607, -0.525725, 0.44722));
    mesh->insert({8}, SMVertex(-0.723607, 0.525725, 0.44722));
    mesh->insert({9}, SMVertex(0.276388, 0.850649, 0.44722));
    mesh->insert({10}, SMVertex(0.894426, 0, 0.447216));
    mesh->insert({11}, SMVertex(0, 0, 1));
    mesh->insert({12}, SMVertex(-0.162456, -0.499995, -0.850654));
    mesh->insert({13}, SMVertex(0.425323, -0.309011, -0.850654));
    mesh->insert({14}, SMVertex(0.262869, -0.809012, -0.525738));
    mesh->insert({15}, SMVertex(0.850648, 0, -0.525736));
    mesh->insert({16}, SMVertex(0.425323, 0.309011, -0.850654));
    mesh->insert({17}, SMVertex(-0.52573, 0, -0.850652));
    mesh->insert({18}, SMVertex(-0.688189, -0.499997, -0.525736));
    mesh->insert({19}, SMVertex(-0.162456, 0.499995, -0.850654));
    mesh->insert({20}, SMVertex(-0.688189, 0.499997, -0.525736));
    mesh->insert({21}, SMVertex(0.262869, 0.809012, -0.525738));
    mesh->insert({22}, SMVertex(0.951058, -0.309013, 0));
    mesh->insert({23}, SMVertex(0.951058, 0.309013, 0));
    mesh->insert({24}, SMVertex(0, -1, 0));
    mesh->insert({25}, SMVertex(0.587786, -0.809017, 0));
    mesh->insert({26}, SMVertex(-0.951058, -0.309013, 0));
    mesh->insert({27}, SMVertex(-0.587786, -0.809017, 0));
    mesh->insert({28}, SMVertex(-0.587786, 0.809017, 0));
    mesh->insert({29}, SMVertex(-0.951058, 0.309013, 0));
    mesh->insert({30}, SMVertex(0.587786, 0.809017, 0));
    mesh->insert({31}, SMVertex(0, 1, 0));
    mesh->insert({32}, SMVertex(0.688189, -0.499997, 0.525736));
    mesh->insert({33}, SMVertex(-0.262869, -0.809012, 0.525738));
    mesh->insert({34}, SMVertex(-0.850648, 0, 0.525736));
    mesh->insert({35}, SMVertex(-0.262869, 0.809012, 0.525738));
    mesh->insert({36}, SMVertex(0.688189, 0.499997, 0.525736));
    mesh->insert({37}, SMVertex(0.162456, -0.499995, 0.850654));
    mesh->insert({38}, SMVertex(0.52573, 0, 0.850652));
    mesh->insert({39}, SMVertex(-0.425323, -0.309011, 0.850654));
    mesh->insert({40}, SMVertex(-0.425323, 0.309011, 0.850654));
    mesh->insert({41}, SMVertex(0.162456, 0.499995, 0.850654));
    mesh->insert({0, 12, 13});
    mesh->insert({1, 13, 15});
    mesh->insert({0, 12, 17});
    mesh->insert({0, 17, 19});
    mesh->insert({0, 16, 19});
    mesh->insert({1, 15, 22});
    mesh->insert({2, 14, 24});
    mesh->insert({3, 18, 26});
    mesh->insert({4, 20, 28});
    mesh->insert({5, 21, 30});
    mesh->insert({1, 22, 25});
    mesh->insert({2, 24, 27});
    mesh->insert({3, 26, 29});
    mesh->insert({4, 28, 31});
    mesh->insert({5, 23, 30});
    mesh->insert({6, 32, 37});
    mesh->insert({7, 33, 39});
    mesh->insert({8, 34, 40});
    mesh->insert({9, 35, 41});
    mesh->insert({10, 36, 38});
    mesh->insert({11, 38, 41});
    mesh->insert({36, 38, 41});
    mesh->insert({9, 36, 41});
    mesh->insert({11, 40, 41});
    mesh->insert({35, 40, 41});
    mesh->insert({8, 35, 40});
    mesh->insert({11, 39, 40});
    mesh->insert({34, 39, 40});
    mesh->insert({7, 34, 39});
    mesh->insert({11, 37, 39});
    mesh->insert({33, 37, 39});
    mesh->insert({6, 33, 37});
    mesh->insert({11, 37, 38});
    mesh->insert({32, 37, 38});
    mesh->insert({10, 32, 38});
    mesh->insert({10, 23, 36});
    mesh->insert({23, 30, 36});
    mesh->insert({9, 30, 36});
    mesh->insert({9, 31, 35});
    mesh->insert({28, 31, 35});
    mesh->insert({8, 28, 35});
    mesh->insert({8, 29, 34});
    mesh->insert({26, 29, 34});
    mesh->insert({7, 26, 34});
    mesh->insert({7, 27, 33});
    mesh->insert({24, 27, 33});
    mesh->insert({6, 24, 33});
    mesh->insert({6, 25, 32});
    mesh->insert({22, 25, 32});
    mesh->insert({10, 22, 32});
    mesh->insert({9, 30, 31});
    mesh->insert({21, 30, 31});
    mesh->insert({4, 21, 31});
    mesh->insert({8, 28, 29});
    mesh->insert({20, 28, 29});
    mesh->insert({3, 20, 29});
    mesh->insert({7, 26, 27});
    mesh->insert({18, 26, 27});
    mesh->insert({2, 18, 27});
    mesh->insert({6, 24, 25});
    mesh->insert({14, 24, 25});
    mesh->insert({1, 14, 25});
    mesh->insert({10, 22, 23});
    mesh->insert({15, 22, 23});
    mesh->insert({5, 15, 23});
    mesh->insert({5, 16, 21});
    mesh->insert({16, 19, 21});
    mesh->insert({4, 19, 21});
    mesh->insert({4, 19, 20});
    mesh->insert({17, 19, 20});
    mesh->insert({3, 17, 20});
    mesh->insert({3, 17, 18});
    mesh->insert({12, 17, 18});
    mesh->insert({2, 12, 18});
    mesh->insert({5, 15, 16});
    mesh->insert({13, 15, 16});
    mesh->insert({0, 13, 16});
    mesh->insert({2, 12, 14});
    mesh->insert({12, 13, 14});
    mesh->insert({1, 13, 14});

    for (int i = 1; i < order; ++i)
    {
        mesh = refineMesh(*mesh);
    }

    casc::compute_orientation(*mesh);
    if (getVolume(*mesh) < 0)
    {
        for (auto& data : mesh->get_level<3>())
            data.orientation *= -1;
    }
    return mesh;
}

std::unique_ptr<SurfaceMesh> cube(int order)
{
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    mesh->insert({0}, SMVertex(-1, -1, -1));
    mesh->insert({1}, SMVertex(-1, 1, -1));
    mesh->insert({2}, SMVertex(1, 1, -1));
    mesh->insert({3}, SMVertex(1, -1, -1));

    mesh->insert({4}, SMVertex(-1, -1, 1));
    mesh->insert({5}, SMVertex(-1, 1, 1));
    mesh->insert({6}, SMVertex(1, 1, 1));
    mesh->insert({7}, SMVertex(1, -1, 1));

    mesh->insert({0, 1, 5});
    mesh->insert({0, 4, 5});

    mesh->insert({0, 3, 4});
    mesh->insert({3, 4, 7});

    mesh->insert({2, 3, 7});
    mesh->insert({2, 6, 7});

    mesh->insert({1, 6, 2});
    mesh->insert({1, 5, 6});

    mesh->insert({4, 5, 6});
    mesh->insert({4, 6, 7});

    mesh->insert({0, 1, 3});
    mesh->insert({1, 2, 3});

    for (int i = 1; i < order; ++i)
    {
        mesh = refineMesh(*mesh);
    }

    casc::compute_orientation(*mesh);
    if (getVolume(*mesh) < 0)
    {
        for (auto& data : mesh->get_level<3>())
            data.orientation *= -1;
    }
    return mesh;
}

std::vector<std::unique_ptr<SurfaceMesh>> splitSurfaces(SurfaceMesh& mesh)
{
    // Queue to store the next edges to visit
    std::deque<SurfaceMesh::SimplexID<2>>     frontier;
    // Set of visited edges
    std::set<SurfaceMesh::SimplexID<2>>       visited;

    std::vector<std::unique_ptr<SurfaceMesh>> meshes;

    using SimplexSet = typename casc::SimplexSet<SurfaceMesh>;
    SimplexSet surface;

    int connected_components = 0;
    for (auto edgeID : mesh.get_level_id<2>())
    {
        // This is a new edge if we haven't visited before
        if (visited.find(edgeID) == visited.end())
        {
            ++connected_components;
            frontier.push_back(edgeID);
            surface.insert(edgeID);

            // Traverse across connected edges.
            while (!frontier.empty())
            {
                SurfaceMesh::SimplexID<2> curr = frontier.front();
                if (visited.find(curr) == visited.end())
                {
                    visited.insert(curr);
                    surface.insert(curr);
                    neighbors_up(mesh, curr, std::back_inserter(frontier));
                }
                frontier.pop_front();
            }

            // Extract out the surface
            SimplexSet tmp;
            casc::getClosure(mesh, surface, tmp);
            surface.clear();
            casc::getStar(mesh, tmp, surface);

            std::unique_ptr<SurfaceMesh> newSMPtr(new SurfaceMesh);
            auto& newSM = *newSMPtr;

            util::int_for_each<std::size_t, typename SimplexSet::cLevelIndex>(surfacemesh_detail::CopyHelper<SurfaceMesh>(), mesh, newSM, surface);

            compute_orientation(newSM);
            if (getVolume(newSM) < 0)
            {
                flipNormals(newSM);
            }

            meshes.push_back(std::move(newSMPtr));
            surface.clear();
        }
    }

    return meshes;
}

void cacheNormals(SurfaceMesh& mesh){
    for (auto fID : mesh.get_level_id<3>()) {
        auto norm = getNormal(mesh, fID);
        normalize(norm);
        (*fID).normal = norm;
    }

    for (auto vID: mesh.get_level_id<1>()) {
        Vector norm;
        auto faces = mesh.up(mesh.up(vID));
        for (auto faceID : faces)
        {
            norm += (*faceID).normal;
        }
        normalize(norm);
        (*vID).normal = norm;
    }
}

// http://pub.ist.ac.at/~edels/Papers/1995-J-03-IncrementalBettiNumbers.pdf
std::tuple<bool, int, int, int> getBettiNumbers(SurfaceMesh& mesh){
    bool valid = true;
    int connected_components = mesh.size<1>();
    int holes = 0;
    int voids = 0;

    std::deque<SurfaceMesh::SimplexID<2>> frontier_edges;
    std::set<SurfaceMesh::SimplexID<2>> visited_edges;
    std::set<SurfaceMesh::SimplexID<3>> faces;

    std::set<SurfaceMesh::KeyType> visited_verts;
    SurfaceMesh::SimplexID<2> curr;

    for (auto eid : mesh.get_level_id<2>())
    {
        if (visited_edges.find(eid) == visited_edges.end())
        {
            // Does this surface have a boundary edgee?
            bool hasboundary = false;
            // Does this surface have consistent normals?
            bool orientable = true;
            // Is this surface have weird 3+ face intersections?
            bool pseudo_manifold = true;
            frontier_edges.push_back(eid);

            while (!frontier_edges.empty())
            {
                curr = frontier_edges.front();
                if (visited_edges.find(curr) == visited_edges.end())
                {
                    std::array<SurfaceMesh::KeyType, 2> n = curr.indices();
                    if (visited_verts.find(n[0]) == visited_verts.end()) {
                        // Have never visited n[0]
                        visited_verts.insert(n[0]);
                        --connected_components;

                        if (visited_verts.find(n[1]) == visited_verts.end()) {
                            // Both vertices unvisited
                            visited_verts.insert(n[1]);
                        }
                    }
                    else{
                        if (visited_verts.find(n[1]) == visited_verts.end()) {
                            visited_verts.insert(n[1]);
                            --connected_components;
                        } else {
                            //found a closed cycle
                            ++holes;
                        }
                    }
                    visited_edges.insert(curr);
                    // If on a boundary stop otherwise add neighboring edges to
                    // the queue
                    auto w = mesh.get_cover(curr);
                    if (w.size() == 1) {hasboundary = true;
                    }
                    else if (w.size() == 2)
                    {
                        neighbors_up(mesh, curr, std::back_inserter(frontier_edges));
                    }
                    else {
                        pseudo_manifold = false;
                        valid = false;
                    }
                    neighbors_up(mesh, curr, std::back_inserter(frontier_edges));
                    for(auto n : w ){
                        auto fid = curr.get_simplex_up(n);
                        faces.insert(fid);
                    }
                }
                frontier_edges.pop_front();
            }

            int sz = faces.size();
            if (pseudo_manifold){
                if (orientable && !hasboundary){
                    // Local space should be a sphere
                    ++voids;
                    --sz;
                }
                holes -= sz;
            } else{
                // TODO: Implement better detection of faces which complete 2-cycles s.t. this can be complete.
            }
            faces.clear();
        }
    }

    visited_edges.clear();
    visited_verts.clear();

    return std::make_tuple(valid, connected_components, holes, voids);
}
} // end namespace gamer
