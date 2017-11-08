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

#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <vector>

#include <libraries/casc/include/CASCTraversals.h>
#include <libraries/casc/include/stringutil.h>
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
        angles[0] = angle(a,b,c);
        angles[1] = angle(b,a,c);
        angles[2] = angle(a,c,b);

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
 */
double getVolume(const SurfaceMesh &mesh)
{
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
            std::cerr << "Orientation undefined..." << std::endl;
            // TODO: Add in helpful error
        }

        // Far less efficient but cooler way...
        // norm = getNormalFromTangent(getTangent(mesh, faceID));
        // auto wedge = a^b^c;
        // tmp = std::sqrt(wedge|wedge);

        // auto sgn = dot((a+b+c)/3, norm);
        // if(sgn <= 0) {
        //     tmp = -1*tmp;
        // }

        volume += tmp/6.0;
    }
    return std::abs(volume);
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

void edgeFlip(SurfaceMesh &mesh, SurfaceMesh::SimplexID<2> edgeID)
{
    // Assuming that the mesh is manifold and edge has been vetted for flipping
    auto name = mesh.get_name(edgeID);
    auto up   = mesh.get_cover(edgeID);
    mesh.remove<2>({name[0], name[1]});
    mesh.insert<3>({name[0], up[0], up[1]});
    mesh.insert<3>({name[1], up[0], up[1]});
}

std::vector<SurfaceMesh::SimplexID<2> > selectFlipEdges(
    const SurfaceMesh &mesh,
    bool preserveRidges,
    std::function<bool(const SurfaceMesh &, SurfaceMesh::SimplexID<2> &)> &checkFlip)
{

    std::vector<SurfaceMesh::SimplexID<2> >   edgesToFlip;
    casc::NodeSet<SurfaceMesh::SimplexID<2> > ignoredEdges;

    for (auto edgeID : mesh.get_level_id<2>())
    {
        if (!ignoredEdges.count(edgeID))
        {
            auto name = mesh.get_name(edgeID);
            std::pair<Vertex, Vertex> shared;
            shared.first  = *mesh.get_simplex_up({name[0]});
            shared.second = *mesh.get_simplex_up({name[1]});

            std::pair<Vertex, Vertex> notShared;
            auto up = mesh.get_cover(edgeID);

            if (up.size() > 2)
            {
                //std::cerr << "This edge participates in more than 2 faces.
                // Returning..." << std::endl;
                continue;
            }
            else if (up.size() < 2)
            {
                //std::cerr << "This edge participates in fewer than 2 faces.
                // Returning..." << std::endl;
                continue;
            }
            notShared.first  = *mesh.get_simplex_up({up[0]});
            notShared.second = *mesh.get_simplex_up({up[1]});

            // Add check to see if notShared.first and second are connected.
            if (mesh.exists<2>({up[0], up[1]}))
            {
                //std::cerr << "Found a tetrahedron cannot edge flip." <<
                // std::endl;
                continue;
            }

            // Check if we're on a ridge
            if (preserveRidges)
            {
                // auto t1 = getTangent(mesh, mesh.get_simplex_up(edgeID,
                // up[0]));
                // auto a = getNormalFromTangent(t1);
                // auto t2 = getTangent(mesh, mesh.get_simplex_up(edgeID,
                // up[1]));
                // auto b = getNormalFromTangent(t2);
                auto a   = getNormal(mesh, mesh.get_simplex_up(edgeID, up[0]));
                auto b   = getNormal(mesh, mesh.get_simplex_up(edgeID, up[1]));
                auto val = angle(a, b);
                if (val > 60)
                {
                    continue;
                }
            }

            // Check if flipping creates a fold
            auto f1   = (shared.first - notShared.first)^(shared.second - notShared.first);
            auto f2   = (shared.first - notShared.second)^(shared.second - notShared.second);
            auto f3   = (notShared.first - shared.first)^(notShared.second - shared.first);
            auto f4   = (notShared.first - shared.second)^(notShared.second - shared.second);
            auto area = std::pow(std::sqrt(f1|f1) + std::sqrt(f2|f2), 2);
            auto areaFlip = std::pow(std::sqrt(f3|f3) + std::sqrt(f4|f4), 2);
            if (areaFlip/area > 1.01) // TODO: this is an arbitrary area
                                      // ratio...
            {                         //std::cerr << "Suspect flipping will
                                      // create fold."  <<
                // std::endl;
                continue;
            }

            if (checkFlip(mesh, edgeID))
            {
                edgesToFlip.push_back(edgeID);
                std::vector<SurfaceMesh::SimplexID<2> > neighbors;
                std::vector<SurfaceMesh::SimplexID<2> > neighborsAway;
                neighbors_up(mesh, edgeID, std::back_inserter(neighbors));
                for (auto neighbor : neighbors)
                {
                    ignoredEdges.insert(neighbor);
                    neighbors_up(mesh, neighbor, std::back_inserter(neighborsAway));
                }
                for (auto neighbor : neighborsAway)
                {
                    ignoredEdges.insert(neighbor);
                }
            }
        }
    }

    return edgesToFlip;
}

bool checkFlipAngle(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID)
{
    auto getMinAngle = [](const Vertex &a, const Vertex &b, const Vertex &c){
            double              minAngle = 999; // dummy for now
            double              tmp;
            std::vector<Vertex> triangle = {a, b, c};
            for (int i = 0; i < 3; i++)
            {
                std::rotate(triangle.begin(), triangle.begin()+i, triangle.end());
                auto it = triangle.begin();
                tmp = angle(*it, *(it+1), *(it+2));
                if (tmp < minAngle)
                    minAngle = tmp;
            }
            return minAngle;
        };

    auto name = mesh.get_name(edgeID);
    std::pair<Vertex, Vertex> shared;
    shared.first  = *mesh.get_simplex_up({name[0]});
    shared.second = *mesh.get_simplex_up({name[1]});
    std::pair<Vertex, Vertex> notShared;
    auto up = mesh.get_cover(edgeID);
    notShared.first  = *mesh.get_simplex_up({up[0]});
    notShared.second = *mesh.get_simplex_up({up[1]});

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

    if (minAngleFlip > minAngle)
        return true;
    return false;
}

bool checkFlipValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<2> &edgeID)
{
    auto name = mesh.get_name(edgeID);
    std::pair<SurfaceMesh::SimplexID<1>, SurfaceMesh::SimplexID<1> > shared;
    shared.first  = mesh.get_simplex_up({name[0]});
    shared.second = mesh.get_simplex_up({name[1]});
    std::pair<SurfaceMesh::SimplexID<1>, SurfaceMesh::SimplexID<1> > notShared;
    auto up = mesh.get_cover(edgeID);
    notShared.first  = mesh.get_simplex_up({up[0]});
    notShared.second = mesh.get_simplex_up({up[1]});
    std::array<double, 20> valence;
    // assuming there are no boundaries
    // TODO check if it's a boundary...
    valence[0] = getValence(mesh, shared.first)-6;
    valence[1] = getValence(mesh, shared.second)-6;
    valence[2] = getValence(mesh, notShared.first)-6;
    valence[3] = getValence(mesh, notShared.second)-6;

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
    int flipExcess = 0;

    for (int i = 0; i < 4; i++)
    {
        flipExcess += std::pow(valence[i], 2);
    }
    if (flipExcess < excess)
        return true;
    return false;
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
    Vector parallel;
    tensor<double, 3, 2> identity{{
                                      1, 0, 0, 0, 1, 0, 0, 0, 1
                                  }};
    auto llproj = identity-perpProj; // perpendicular projector
    Eigen::Map<Eigen::Matrix3d> llproj_e(llproj.data());
    Eigen::Map<Eigen::Vector3d> parallel_e(parallel.data());
    parallel_e = llproj_e*disp_e;

    (*vertexID).position = (*vertexID).position + parallel;
}


void normalSmooth(SurfaceMesh &mesh){
    for(auto nid : mesh.get_level_id<1>()){
        normalSmoothH(mesh, nid);
    }
    double min, max;
    int nSmall, nLarge;
    std::tie(min,max, nSmall, nLarge) = getMinMaxAngles(mesh, 15, 150);
    std::cout << "  Min Angle: " << min << ", Max Angle: " << max 
              << ", # Small Angles: " << nSmall 
              << ", # Large Angles: " << nLarge << std::endl; 
}

/**
 * @brief      Perona-Malik normal based smoothing algorithm
 *
 * @param      mesh      The mesh
 * @param[in]  vertexID  The vertex id
 */
void normalSmoothH(SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    auto            name = mesh.get_name(vertexID)[0];
    double          areaSum = 0;

    auto            p = (*vertexID).position;
    Eigen::Vector4d pos_e;
    pos_e << p[0], p[1], p[2], 1;
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

        for (auto face : faces)
        {
            auto inorm = getNormal(mesh, face);
            avgNorm += inorm;
        }
        avgNorm /= 3;                          // each triangle has 3 incident
                                               // faces
        avgNorm /= std::sqrt(avgNorm|avgNorm); // get unit normal

        // Compute the edge (axis) to rotate about.
        auto edge = mesh.get_simplex_down(faceID, name);
        auto edgeName = mesh.get_name(edge);
        auto a  = *mesh.get_simplex_up({edgeName[0]});
        auto b  = *mesh.get_simplex_up({edgeName[1]});
        auto ab = a-b;
        ab /= std::sqrt(ab|ab); // Eigen AngleAxis requires unit vector


        // Angle between normals in radians. This is the angle to rotate the
        // normal by.
        double angle = std::copysign(std::acos(normal|avgNorm), dot(cross(normal, avgNorm), ab));
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

int getValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> SimplexID)
{
    std::vector<SurfaceMesh::SimplexID<1> > vertices;
    neighbors_up(mesh, SimplexID, std::back_inserter(vertices));
    return vertices.size();
}

tensor<double, 3, 2> getTangent(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<1> vertexID)
{
    return getTangentH(mesh, (*vertexID).position, vertexID);
}

tensor<double, 3, 2> getTangent(const SurfaceMesh &mesh, SurfaceMesh::SimplexID<3> faceID)
{
    auto cover = mesh.get_name(faceID);
    auto vertexID = mesh.get_simplex_up({cover[0]});
    std::set<SurfaceMesh::KeyType> next(std::begin(cover)+1, std::end(cover));
    return getTangentF(mesh, (*vertexID).position, vertexID, next);
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
        std::cerr << "Orientation undefined..." << std::endl;
    }
    return norm;
}

Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double, 3, 2> mat)
{
    Eigen::Map<Eigen::Matrix3d> emat(mat.data());
    // TODO: How much optimization can we get from having a persistent
    // eigensolver?
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(emat);
    if (eigensolver.info() != Eigen::Success)
        abort();
    return eigensolver;
}
