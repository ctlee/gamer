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

#include <libraries/casc/casc>
#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>

#include "SurfaceMesh.h"
#include "Vertex.h"

void print(const SurfaceMesh &mesh)
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

void generateHistogram(const SurfaceMesh &mesh)
{
    // compute angle distribution
    std::array<double, 18> histogram;
    histogram.fill(0);
    for (auto face : mesh.get_level_id<3>())
    {
        auto   vertexIDs = mesh.get_name(face);
        // Unpack the ID's for convenience
        Vertex a = *mesh.get_simplex_up<1>({vertexIDs[0]});
        Vertex b = *mesh.get_simplex_up<1>({vertexIDs[1]});
        Vertex c = *mesh.get_simplex_up<1>({vertexIDs[2]});

        auto   binAngle = [&](double angle) -> int{
                return std::floor(angle/10);
            };
        histogram[binAngle(angle(a, b, c))]++;
        histogram[binAngle(angle(b, a, c))]++;
        histogram[binAngle(angle(c, a, b))]++;
    }

    int factor = mesh.size<3>()*3;
    std::for_each(histogram.begin(), histogram.end(), [&factor](double &n){
        n = 100.0*n/factor;
    });

    std::cout << "Angle Distribution:" << std::endl;
    for (int x = 0; x < 18; x++)
        std::cout << x*10 << "-" << (x+1)*10 << ": " << std::setprecision(2)
                  << std::fixed << histogram[x] << std::endl;
    std::cout << std::endl << std::endl;

    // compute the edge length distribution
    std::cout << "Edge Length Distribution:" << std::endl;
    std::vector<double> lengths;
    for (auto edge : mesh.get_level_id<2>())
    {
        auto   vertexIDs = mesh.down(edge);
        auto   t1  = *vertexIDs.cbegin();
        auto   t2  = *(++vertexIDs.cbegin());
        auto   v1  = *t1;
        auto   v2  = *t2;
        double len = magnitude(v2-v1);
        lengths.push_back(len);
    }
    std::sort(lengths.begin(), lengths.end());

    std::array<double, 20> histogramLength;
    histogramLength.fill(0);
    double                 interval = (lengths.back() - lengths.front())/20;
    double                 low = lengths.front();

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
        std::for_each(histogramLength.begin(), histogramLength.end(), [&factor](double &n){
            n = 100.0*n/factor;
        });

        for (int x = 0; x < 20; x++)
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

std::tuple<double, double, int, int> getMinMaxAngles(
    const SurfaceMesh &mesh,
    int maxMinAngle,
    int minMaxAngle)
{
    double minAngle = 360;
    double maxAngle = 0;
    int small = 0;
    int large = 0;

    // for each triangle
    for(auto nid : mesh.get_level_id<3>()){
        auto name = mesh.get_name(nid);
        auto a = *mesh.get_simplex_up({name[0]});
        auto b = *mesh.get_simplex_up({name[1]});
        auto c = *mesh.get_simplex_up({name[2]});
        std::array<double, 3> angles;
        try{
            angles[0] = angle(a,b,c);
            angles[1] = angle(b,a,c);
            angles[2] = angle(a,c,b);
        }
        catch (std::runtime_error& e){
            throw std::runtime_error("ERROR(getMinMaxAngles): Cannot compute angles of face with zero area.");
        }

        for(double angle : angles){
            if (angle < minAngle){
                minAngle = angle;
            }
            if (angle > maxAngle){
                maxAngle = angle;
            }
            if (angle < maxMinAngle){
                ++small;
            }
            if (angle > minMaxAngle){
                ++large;
            }
        }
    }
    return std::make_tuple(minAngle, maxAngle, small, large);
}

double getArea(const SurfaceMesh &mesh)
{
    double area;
    for (auto faceID : mesh.get_level_id<3>())
        area += getArea(mesh, faceID);
    return area;
}

double getArea(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID)
{
    auto name = mesh.get_name(faceID);
    auto a = *mesh.get_simplex_up({name[0]});
    auto b = *mesh.get_simplex_up({name[1]});
    auto c = *mesh.get_simplex_up({name[2]});
    auto wedge = (b-a)^(b-c);
    return std::sqrt(wedge|wedge)/2;
}

/**
 * http://research.microsoft.com/en-us/um/people/chazhang/publications/icip01_ChaZhang.pdf
 * http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
 */
double getVolume(const SurfaceMesh &mesh)
{
    bool orientError = false;
    double volume = 0;
    for (auto faceID : mesh.get_level_id<3>())
    {
        auto   name = mesh.get_name(faceID);
        auto   a = (*mesh.get_simplex_up({name[0]})).position;
        auto   b = (*mesh.get_simplex_up({name[1]})).position;
        auto   c = (*mesh.get_simplex_up({name[2]})).position;

        Vector norm;
        double tmp;
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
    if(orientError){
        std::cerr << "ERROR getVolume(): Orientation undefined for one or more "
                  << "simplices. Did you call compute_orientation()?" << std::endl;
    }
    return volume/6;
}

void translate(SurfaceMesh &mesh, Vector v)
{
    for (auto &vertex : mesh.get_level<1>())
        vertex += v;
}

void translate(SurfaceMesh &mesh, double dx, double dy, double dz)
{
    Vector v = Vector({dx, dy, dz});
    translate(mesh, v);
}

void scale(SurfaceMesh &mesh, Vector v)
{
    for (auto &vertex : mesh.get_level<1>())
    {
        vertex[0] *= v[0];
        vertex[1] *= v[1];
        vertex[2] *= v[2];
    }
}

void scale(SurfaceMesh &mesh, double sx, double sy, double sz)
{
    Vector v = Vector();
    v[0] = sx; v[1] = sy; v[2] = sz;
    scale(mesh, v);
}

void scale(SurfaceMesh &mesh, double s)
{
    Vector v = Vector();
    v[0] = s; v[1] = s; v[2] = s;
    scale(mesh, v);
}

std::pair<Vector, double> getCenterRadius(SurfaceMesh &mesh){
    Vector center = Vector();
    double radius = 0;
    if (mesh.size<1>() > 0){
        for(auto vertexID : mesh.get_level_id<1>()){
            center += *vertexID;
        }
        center /= mesh.size<1>();

        for(auto vertexID : mesh.get_level_id<1>()){
            Vector tmp = *vertexID - center;
            double dist = std::sqrt(tmp|tmp);
            if(dist > radius) radius = dist;
        }
    }
    return std::make_pair(center, radius);
}

void centeralize(SurfaceMesh &mesh){
    Vector center;
    double radius;
    std::tie(center, radius) = getCenterRadius(mesh);
    translate(mesh, -center);
}

void normalSmooth(SurfaceMesh &mesh){
    for(auto nid : mesh.get_level_id<1>()){
        surfacemesh_detail::normalSmoothH(mesh, nid);
    }
    double min, max;
    int nSmall, nLarge;
    std::tie(min,max, nSmall, nLarge) = getMinMaxAngles(mesh, 15, 150);
    std::cout << "  Min Angle: " << min << ", Max Angle: " << max
              << ", # Small Angles: " << nSmall
              << ", # Large Angles: " << nLarge << std::endl;
}

int getValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID)
{
    return mesh.get_cover(vertexID).size();
}

tensor<double, 3, 2> getTangent(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    return surfacemesh_detail::getTangentH(mesh, (*vertexID).position, vertexID);
}

tensor<double, 3, 2> getTangent(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID)
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

Vector getNormal(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    Vector norm;
    auto   faces = mesh.up(mesh.up(vertexID));
    for (auto faceID : faces)
    {
        norm += getNormal(mesh, faceID);
    }
    norm /= faces.size();
    return norm;
}

Vector getNormal(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID)
{
    Vector norm;
    auto   name = mesh.get_name(faceID);

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
        // std::cerr << "ERROR(getNormal): Orientation undefined, cannot compute "
                  // << "normal. Did you call compute_orientation()?" << std::endl;
        throw std::runtime_error("ERROR(getNormal): Orientation undefined, cannot compute normal. Did you call compute_orientation()?");
    }
    return norm;
}


// TODO: (0) Allow for smoothing of specific regions.
bool smoothMesh(SurfaceMesh &mesh, int maxMinAngle, int minMaxAngle, int maxIter, bool preserveRidges){
    bool smoothed;
    double minAngle, maxAngle;
    int nSmall, nLarge;
    int nIter = 1;

    // Check if we are smoothed already
    smoothed = minAngle > maxMinAngle && maxAngle < minMaxAngle;

    std::tie(minAngle, maxAngle, nSmall, nLarge) = getMinMaxAngles(mesh, maxMinAngle, minMaxAngle);
    std::cout << "Initial Quality: Min Angle = " << minAngle << ", "
              << "Max Angle = " << maxAngle << ", "
              << "# smaller-than-" << maxMinAngle << " = " << nSmall << ", "
              << "# larger-than-" << minMaxAngle << " = " << nLarge << std::endl;

    // while not smoothed and not at maxiter perform one round of
    // weightedVertexSmooth + edgeflipping
    while (!smoothed && nIter <= maxIter){
        for(auto vertex : mesh.get_level_id<1>()){
            surfacemesh_detail::weightedVertexSmooth(mesh, vertex, RINGS);
            //barycenterVertexSmooth(mesh, vertex);
        }

        std::vector<SurfaceMesh::SimplexID<2> > edgesToFlip;
        // Get set of good, non-interfering edges to flip according to the
        // Angle based criteria.
        surfacemesh_detail::selectFlipEdges(mesh, preserveRidges, surfacemesh_detail::checkFlipAngle,
                            std::back_inserter(edgesToFlip));
        for(auto edgeID : edgesToFlip){
            surfacemesh_detail::edgeFlip(mesh, edgeID);
        }
        init_orientation(mesh);
        check_orientation(mesh);

        // Mark for flipping by edge valence.
        // edgesToFlip.clear();
        // selectFlipEdges(mesh, preserveRidges, checkFlipValence,
        //                    std::back_inserter(edgesToFlip));
        // for(auto edgeID : edgesToFlip){
        //     edgeFlip(mesh, edgeID);
        // }
        // init_orientation(mesh);
        // check_orientation(mesh);

        std::tie(minAngle, maxAngle, nSmall, nLarge) = getMinMaxAngles(mesh, maxMinAngle, minMaxAngle);
        std::cout << "Iteration " << nIter << ":" << std::endl;
        std::cout << "Min Angle = " << minAngle << ", "
                  << "Max Angle = " << maxAngle << ", "
                  << "# smaller-than-" << maxMinAngle << " = " << nSmall << ", "
                  << "# larger-than-" << minMaxAngle << " = " << nLarge << std::endl;

        smoothed = minAngle > maxMinAngle && maxAngle < minMaxAngle;
        ++nIter;
    }
    return smoothed;
}

/**
 * @brief      Refine the mesh by quadrisection.
 *
 * Note that this function will delete all stored data on edges and faces. But
 * this can be easily fixed.
 *
 * @param      mesh  The mesh
 */
std::unique_ptr<SurfaceMesh> refineMesh(const SurfaceMesh &mesh){
    std::unique_ptr<SurfaceMesh> refinedMesh(new SurfaceMesh);

    // Copy over vertices to refinedMesh
    for (auto vertex : mesh.get_level_id<1>()){
        auto key = mesh.get_name(vertex);
        refinedMesh->insert(key, *vertex);
    }

    // Split edges and generate a map of names before to after
    std::map<std::array<int,2>, int> edgeMap;

    for(auto edge : mesh.get_level_id<2>()){
        auto edgeName = mesh.get_name(edge);
        auto v1 = *mesh.get_simplex_up({edgeName[0]});
        auto v2 = *mesh.get_simplex_up({edgeName[1]});

        auto newVertex = refinedMesh->add_vertex(Vertex(0.5*(v1+v2)));
        edgeMap.emplace(std::make_pair(edgeName, newVertex));
    }

    // Connect edges and face and copy data
    for(auto face : mesh.get_level_id<3>()){
        auto name = mesh.get_name(face);
        int a, b, c;

        // Skip checking if found
        auto it = edgeMap.find({name[0],name[1]});
        a = it->second;

        it = edgeMap.find({name[1],name[2]});
        b = it->second;

        it = edgeMap.find({name[0],name[2]});
        c = it->second;

        refinedMesh->insert({name[0], a});
        refinedMesh->insert({a, name[1]});

        refinedMesh->insert({name[1], b});
        refinedMesh->insert({b, name[2]});

        refinedMesh->insert({name[0], c});
        refinedMesh->insert({c, name[2]});

        refinedMesh->insert({a,b,c});
        refinedMesh->insert({name[0],a,c}, *face);
        refinedMesh->insert({name[1],a,b}, *face);
        refinedMesh->insert({name[2],b,c}, *face);
    }
    return refinedMesh;
}


void coarse(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight){
    // TODO: Check if all polygons are closed (0)
    std::cout << "Before coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;

    // Compute the average edge length
    double avgLen = 0;
    if (denseWeight > 0){
        avgLen = surfacemesh_detail::getMeanEdgeLength(mesh);
    }

    double sparsenessRatio = 1;
    double flatnessRatio = 1;

    // Construct list of vertices to decimate
    std::vector<SurfaceMesh::SimplexID<1> > toRemove;
    for(auto vertexID : mesh.get_level_id<1>()){
        // Sparseness as coarsening criteria
        if(denseWeight > 0){
            // Get max length of edges.
            auto edges = mesh.up(vertexID);
            double maxLen = 0;
            for(auto edgeID : edges){
                auto name =  mesh.get_name(edgeID);
                auto v = *mesh.get_simplex_down(edgeID, name[0])
                         - *mesh.get_simplex_down(edgeID, name[1]);
                double tmpLen = std::sqrt(v|v);
                if(tmpLen > maxLen) maxLen = tmpLen;
            }
            sparsenessRatio = std::pow(maxLen/avgLen, denseWeight);
        }

        // Curvature as coarsening criteria
        if(flatRate > 0){
            auto lst = surfacemesh_detail::computeLocalStructureTensor(mesh, vertexID, RINGS);
            auto eigenvalues = surfacemesh_detail::getEigenvalues(lst).eigenvalues();
            // The closer this ratio is to 0 the flatter the local region.
            flatnessRatio = std::pow(eigenvalues[1]/eigenvalues[2], flatRate);
        }

        // Add vertex to delete list
        if(sparsenessRatio * flatnessRatio < coarseRate){
            toRemove.push_back(vertexID);
        }
    }

    std::cout << toRemove.size() << " vertices are marked to be removed." << std::endl;

    for(auto vertexID : toRemove){
        surfacemesh_detail::decimateVertex(mesh, vertexID);
    }

    std::cout << "After coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;
}


void coarseIT(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight){
    std::cout << "Before coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;

    // Compute the average edge length
    double avgLen = 0;
    if (denseWeight > 0){
        avgLen = surfacemesh_detail::getMeanEdgeLength(mesh);
    }

    double sparsenessRatio = 1;
    double flatnessRatio = 1;

    // Backup vertices so we can modify in place
    auto range = mesh.get_level_id<1>();
    const std::vector<SurfaceMesh::SimplexID<1> > Vertices(range.begin(), range.end());
    for(auto vertexID : Vertices) {
        // Sparseness as coarsening criteria
        if(denseWeight > 0){
            // Get max length of edges.
            auto edges = mesh.up(vertexID);
            double maxLen = 0;
            for(auto edgeID : edges){
                auto name =  mesh.get_name(edgeID);
                auto v = *mesh.get_simplex_down(edgeID, name[0])
                         - *mesh.get_simplex_down(edgeID, name[1]);
                double tmpLen = std::sqrt(v|v);
                if(tmpLen > maxLen) maxLen = tmpLen;
            }
            sparsenessRatio = std::pow(maxLen/avgLen, denseWeight);
        }

        // Curvature as coarsening criteria
        if(flatRate > 0){
            auto lst = surfacemesh_detail::computeLocalStructureTensor(mesh, vertexID, RINGS);
            auto eigenvalues = surfacemesh_detail::getEigenvalues(lst).eigenvalues();
            // The closer this ratio is to 0 the flatter the local region.
            flatnessRatio = std::pow(eigenvalues[1]/eigenvalues[2], flatRate);
        }

        // Check if we should delete
        if(sparsenessRatio * flatnessRatio < coarseRate){
            surfacemesh_detail::decimateVertex(mesh, vertexID);
        }
    }

    std::cout << "After coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;
}

void fillHoles(SurfaceMesh &mesh){
    std::vector<std::vector<SurfaceMesh::SimplexID<2>>> holeList;
    surfacemesh_detail::findHoles(mesh, holeList);

    for(auto& holeEdges : holeList){
        std::vector<SurfaceMesh::SimplexID<1>> sortedVertices;
        surfacemesh_detail::edgeRingToVertices(mesh, holeEdges, std::back_inserter(sortedVertices));

        surfacemesh_detail::triangulateHole(mesh, sortedVertices, Face(), holeEdges);
    }
}

void flipNormals(SurfaceMesh &mesh){
    for(auto& fdata : mesh.get_level<3>()){
        fdata.orientation *= -1;
    }
}

bool hasHole(const SurfaceMesh &mesh){
    for(auto eID : mesh.get_level_id<2>()){
        auto cover = mesh.get_cover(eID);
        if(cover.size() != 2){
            return true;
        }
    }
    return false;
}

std::unique_ptr<SurfaceMesh> sphere(int order){
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    mesh->insert({0}, Vertex(0,0,-1));
    mesh->insert({1}, Vertex(-1,-1,0));
    mesh->insert({2}, Vertex(-1,1,0));
    mesh->insert({3}, Vertex(1,1,0));
    mesh->insert({4}, Vertex(1,-1,0));
    mesh->insert({5}, Vertex(0,0,1));

    mesh->insert({0,1,2});
    mesh->insert({0,2,3});
    mesh->insert({0,3,4});
    mesh->insert({0,1,4});
    mesh->insert({5,1,2});
    mesh->insert({5,2,3});
    mesh->insert({5,3,4});
    mesh->insert({5,1,4});
    for (int i = 1; i < order; ++i){
        mesh = refineMesh(*mesh);
    }

    compute_orientation(*mesh);
    if(getVolume(*mesh) < 0){
        for(auto &data : mesh->get_level<3>())
            data.orientation *= -1;
    }
    return mesh;
}

std::unique_ptr<SurfaceMesh> cube(int order){
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    mesh->insert({0}, Vertex(-1,-1,-1));
    mesh->insert({1}, Vertex(-1,1,-1));
    mesh->insert({2}, Vertex(1,1,-1));
    mesh->insert({3}, Vertex(1,-1,-1));

    mesh->insert({4}, Vertex(-1,-1,1));
    mesh->insert({5}, Vertex(-1,1,1));
    mesh->insert({6}, Vertex(1,1,1));
    mesh->insert({7}, Vertex(1,-1,1));

    mesh->insert({0,1,5});
    mesh->insert({0,4,5});

    mesh->insert({0,3,4});
    mesh->insert({3,4,7});

    mesh->insert({2,3,7});
    mesh->insert({2,6,7});

    mesh->insert({1,6,2});
    mesh->insert({1,5,6});

    mesh->insert({4,5,6});
    mesh->insert({4,6,7});

    mesh->insert({0,1,3});
    mesh->insert({1,2,3});

    for (int i = 1; i < order; ++i){
        mesh = refineMesh(*mesh);
    }

    compute_orientation(*mesh);
    if(getVolume(*mesh) < 0){
        for(auto &data : mesh->get_level<3>())
            data.orientation *= -1;
    }
    return mesh;
}

