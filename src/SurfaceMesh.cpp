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

#include <libraries/casc/include/CASCFunctions.h>
#include <libraries/casc/include/CASCTraversals.h>
#include <libraries/casc/include/stringutil.h>
#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>

#include <libraries/casc/include/typetraits.h>
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
        if ((*faceID).orientation == -1)
        {
            // a->b->c
            norm = cross(b, c);
            tmp  = dot(a, norm);

        }
        else if ((*faceID).orientation == 1)
        {
            // c->b->a
            norm = cross(b, a);
            tmp  = dot(c, norm);
        }
        else
        {
            orientError = true;
        }

        // Far less efficient but cooler way...
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

void edgeFlip(SurfaceMesh &mesh, SurfaceMesh::SimplexID<2> edgeID)
{
    // Assuming that the mesh is manifold and edge has been vetted for flipping
    auto name = mesh.get_name(edgeID);
    auto up   = mesh.get_cover(edgeID);
    mesh.remove<2>({name[0], name[1]});
    mesh.insert<3>({name[0], up[0], up[1]});
    mesh.insert<3>({name[1], up[0], up[1]});
    // TODO: this gets rid of boundary markings...
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
    // TODO: (5) check if it's a boundary...
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

        for (auto face : faces)
        {
            auto inorm = getNormal(mesh, face);
            inorm /= std::sqrt(inorm|inorm);
            avgNorm += inorm;
        }
        avgNorm /= 3;

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

int getValence(const SurfaceMesh &mesh, const SurfaceMesh::SimplexID<1> vertexID)
{
    return mesh.get_cover(vertexID).size();
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

    if ((*faceID).orientation == -1)
    {
        norm = cross(c-b, a-b);
    }
    else if ((*faceID).orientation == 1)
    {
        norm = cross(a-b, c-b);
    }
    else
    {
        std::cerr << "ERROR(getNormal): Orientation undefined, cannot compute "
                  << "normal. Did you call compute_orientation()?" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return norm;
}

Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> getEigenvalues(tensor<double, 3, 2> mat)
{
    Eigen::Map<Eigen::Matrix3d> emat(mat.data());
    // TODO: (99) How much optimization can we get from having a persistent eigensolver?
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(emat);
    if (eigensolver.info() != Eigen::Success)
        abort();
    return eigensolver;
}

void weightedVertexSmooth(SurfaceMesh &mesh,
    SurfaceMesh::SimplexID<1> vertexID,
    int rings)
{
    // TODO: (2) Fix problems when manipulating boundary vertices.
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
        auto edgeName = mesh.get_name(edge);
        Vertex shared   = *mesh.get_simplex_up({(edgeName[0] == centerName) ? edgeName[1] : edgeName[0]});

        // Get the vertices connected to adjacent edge by face
        auto up   = mesh.get_cover(edge);

        Vertex prev = *mesh.get_simplex_up({up[0]});
        Vertex next = *mesh.get_simplex_up({up[1]});

        Vector pS = prev - shared;
        pS  /= std::sqrt(pS|pS);
        Vector nS = next - shared;
        nS /= std::sqrt(nS|nS);

        // Bisector of the 'rhombus'
        Vector bisector = (pS + nS)/2;
        // TODO (0): Create a better method for catching divide by zeros

        try {
            normalize(bisector);
        }
        catch (std::runtime_error& e){
            continue;
        }

        // Get a reference vecter to shared which lies on the plane of interest.
        Vector disp = center - shared;
        Eigen::Map<Eigen::Vector3d> disp_e(disp.data());

        // Normal of tangent plane
        //auto tanNorm = getNormalFromTangent(pS^nS);
        Vector tanNorm = cross(pS, nS);

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
        perp_e = perpProj_e*disp_e; // matrix (3x3) * vector = vector

        auto alpha = (pS|nS)+1; // keep the dot product positive
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
    auto lst = computeLocalStructureTensor(mesh, vertexID, rings);

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
             ).matrix(); // vector 3x1
    newPos_e = eigen_result.eigenvectors()*w; // matrix 3x3 * vector = vector
    center.position += newPos;
}

tensor<double,3,2> computeLocalStructureTensor(const SurfaceMesh &mesh,
        const SurfaceMesh::SimplexID<1> vertexID,
        const int rings){
    // Set of neighbors
    std::set<SurfaceMesh::SimplexID<1> > nbors;
    // Get list of neighbors
    casc::kneighbors_up(mesh, vertexID, rings, nbors);
    // local structure tensor
    tensor<double, 3, 2> lst = tensor<double, 3, 2>();
    for (SurfaceMesh::SimplexID<1> nid : nbors)
    {
        auto norm = getNormal(mesh, nid);           // Get Vector normal
        norm /= std::sqrt(norm|norm);               // normalize
        lst += norm*norm;                           // tensor product
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
            weightedVertexSmooth(mesh, vertex, RINGS);
            //barycenterVertexSmooth(mesh, vertex);
        }

        std::vector<SurfaceMesh::SimplexID<2> > edgesToFlip;
        // Get set of good, non-interfering edges to flip according to the
        // Angle based criteria.
        selectFlipEdges(mesh, preserveRidges, checkFlipAngle,
                            std::back_inserter(edgesToFlip));
        for(auto edgeID : edgesToFlip){
            edgeFlip(mesh, edgeID);
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
        for (auto edgeID : mesh.get_level_id<2>()){
            auto name =  mesh.get_name(edgeID);
            auto v = *mesh.get_simplex_down(edgeID, name[0])
                     - *mesh.get_simplex_down(edgeID, name[1]);
            avgLen += std::sqrt(v|v);
        }
        avgLen /= mesh.size<2>();
    }

    double sparsenessRatio = 1;
    double flatnessRatio = 1;

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
            // TODO: how many rings to consider? (0)
            auto lst = computeLocalStructureTensor(mesh, vertexID, RINGS);
            auto eigenvalues = getEigenvalues(lst).eigenvalues();
            // std::cout << "E[0]: " << eigenvalues[0] << std::endl;
            // std::cout << "E[1]: " << eigenvalues[1] << std::endl;
            // std::cout << "E[2]: " << eigenvalues[2] << std::endl;
            // The closer this ratio is to 0 the flatter the local region.
            flatnessRatio = std::pow(eigenvalues[1]/eigenvalues[2], flatRate);
        }

        //std::cout << "Coarse value: " << sparsenessRatio*flatnessRatio << std::endl;
        // Add vertex to delete list
        if(sparsenessRatio * flatnessRatio < coarseRate){
            toRemove.push_back(vertexID);
        }
    }

    std::cout << toRemove.size() << " vertices are marked to be removed." << std::endl;

    for(auto vertexID : toRemove){
        auto fdata = **mesh.up(std::move(mesh.up(vertexID))).begin();
        fdata.orientation = 0; // Reset the orientation accordingly

        std::set<SurfaceMesh::SimplexID<1> > boundary;
        casc::neighbors_up(mesh, vertexID, std::inserter(boundary, boundary.end()));
        std::set<SurfaceMesh::SimplexID<1>> backupBoundary(boundary);

        // Remove the vertex
        mesh.remove(vertexID);

        // Sort vertices into ring order
        std::vector<SurfaceMesh::SimplexID<1>> sortedVerts;
        std::set<int> bNames; // boundary names
        std::vector<SurfaceMesh::SimplexID<2>> edgeList;

        auto it = boundary.begin();
        int firstName = mesh.get_name(*it)[0];
        int n1;

        while(boundary.size() > 0)
        {
            std::vector<SurfaceMesh::SimplexID<1>> nbors;
            auto currID = *it;

            n1 = mesh.get_name(currID)[0];
            bNames.insert(n1);

            bool success = false;
            std::move(it, std::next(it), std::back_inserter(sortedVerts));
            boundary.erase(it);

            if(boundary.size() == 0){
                if(mesh.exists({n1,firstName})){
                    edgeList.push_back(mesh.get_simplex_up(currID, firstName));
                    break; // Break out of while loop
                }
            }
            else{
                // Get neighbors and search for next vertex
                casc::neighbors_up(mesh, currID, std::back_inserter(nbors));
                for(auto nbor : nbors){
                    auto result = boundary.find(nbor);
                    if(result != boundary.end()){
                        // Check that the edge is a boundary
                        auto tmp = mesh.get_simplex_up(*result, n1);
                        if(mesh.get_cover(tmp).size() == 1){
                            it = result;
                            edgeList.push_back(tmp);
                            success = true;
                            break;
                        }
                    }
                }
            }
            // The ring isn't really a ring
            if(!success){
                std::cerr << "ERROR(coarse): Hole ring is not closed." << std::endl;
                abort(); // Something bad happened...
                //return;
            }
        }

        triangulateHole(mesh, sortedVerts, fdata, std::back_inserter(edgeList));

        // Set the orientation of each edge
        // TODO: (0) Make orientation automatic
        orientHoleHelper<std::integral_constant<std::size_t,1>>::apply(mesh, std::move(bNames), backupBoundary.begin(), backupBoundary.end());

        // Compute the face orientations
        if(!computeHoleOrientation(mesh, std::move(edgeList))){
            std::cerr << "ERROR(coarse): Mesh became non-orientable." << std::endl;
            abort();
            //return;
        }

        // Smooth vertices around the filled hole
        for(auto v : backupBoundary){
            weightedVertexSmooth(mesh, v, RINGS);
        }
    }

    std::cout << "After coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;
}


void coarseIT(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight){
    // TODO: Check if all polygons are closed (0)
    std::cout << "Before coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;

    // Compute the average edge length
    double avgLen = 0;
    if (denseWeight > 0){
        for (auto edgeID : mesh.get_level_id<2>()){
            auto name =  mesh.get_name(edgeID);
            auto v = *mesh.get_simplex_down(edgeID, name[0])
                     - *mesh.get_simplex_down(edgeID, name[1]);
            avgLen += std::sqrt(v|v);
        }
        avgLen /= mesh.size<2>();
    }

    double sparsenessRatio = 1;
    double flatnessRatio = 1;
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
            // TODO: how many rings to consider? (0)
            auto lst = computeLocalStructureTensor(mesh, vertexID, RINGS);
            auto eigenvalues = getEigenvalues(lst).eigenvalues();
            // std::cout << "E[0]: " << eigenvalues[0] << std::endl;
            // std::cout << "E[1]: " << eigenvalues[1] << std::endl;
            // std::cout << "E[2]: " << eigenvalues[2] << std::endl;
            // The closer this ratio is to 0 the flatter the local region.
            flatnessRatio = std::pow(eigenvalues[1]/eigenvalues[2], flatRate);
        }

        //std::cout << "Coarse value: " << sparsenessRatio*flatnessRatio << std::endl;
        // Add vertex to delete list
        if(sparsenessRatio * flatnessRatio < coarseRate){
            auto fdata = **mesh.up(std::move(mesh.up(vertexID))).begin();
            fdata.orientation = 0; // Reset the orientation accordingly

            std::set<SurfaceMesh::SimplexID<1> > boundary;
            casc::neighbors_up(mesh, vertexID, std::inserter(boundary, boundary.end()));
            std::set<SurfaceMesh::SimplexID<1>> backupBoundary(boundary);

            // Remove the vertex
            mesh.remove(vertexID);

            // Sort vertices into ring order
            std::vector<SurfaceMesh::SimplexID<1>> sortedVerts;
            std::set<int> bNames; // boundary names
            std::vector<SurfaceMesh::SimplexID<2>> edgeList;

            auto it = boundary.begin();
            int firstName = mesh.get_name(*it)[0];
            int n1;
            while(boundary.size() > 0)
            {
                std::vector<SurfaceMesh::SimplexID<1>> nbors;
                auto currID = *it;

                n1 = mesh.get_name(currID)[0];
                bNames.insert(n1);

                bool success = false;
                std::move(it, std::next(it), std::back_inserter(sortedVerts));
                boundary.erase(it);

                if(boundary.size() == 0){
                    if(mesh.exists({n1,firstName})){
                        edgeList.push_back(mesh.get_simplex_up(currID, firstName));
                        break; // Break out of while loop
                    }
                }
                else{
                    // Get neighbors and search for next vertex
                    casc::neighbors_up(mesh, currID, std::back_inserter(nbors));
                    for(auto nbor : nbors){
                        auto result = boundary.find(nbor);
                        if(result != boundary.end()){
                            // Check that the edge is a boundary
                            auto tmp = mesh.get_simplex_up(*result, n1);
                            if(mesh.get_cover(tmp).size() == 1){
                                // std::cout << "Sorted: " << tmp << std::endl;
                                it = result;
                                edgeList.push_back(tmp);
                                success = true;
                                break;
                            }
                        }
                    }
                }
                // The ring isn't really a ring
                if(!success){
                    std::cerr << "ERROR(coarse): Hole ring is not closed." << std::endl;
                    abort(); // Something bad happened...
                    //return;
                }
            }

            triangulateHole(mesh, sortedVerts, fdata, std::back_inserter(edgeList));

            // Set the orientation of each edge
            // TODO: (0) Make orientation automatic
            orientHoleHelper<std::integral_constant<std::size_t,1>>::apply(mesh, std::move(bNames), backupBoundary.begin(), backupBoundary.end());

            // TODO: COMPUTE THE FACE ORIENTATIONS!.....
            if(!computeHoleOrientation(mesh, std::move(edgeList))){
                std::cerr << "ERROR(coarse): Mesh became non-orientable." << std::endl;
                return;
            }

            // Smooth vertices around the filled hole
            for(auto v : backupBoundary){
                weightedVertexSmooth(mesh, v, RINGS);
            }
        }
    }

    std::cout << "After coarsening: " << mesh.size<1>() << " " << mesh.size<2>() << " " << mesh.size<3>() << std::endl;
}

bool computeHoleOrientation(SurfaceMesh &mesh, const std::vector<SurfaceMesh::SimplexID<2> > &&edgeList){
    std::deque<SurfaceMesh::SimplexID<2> > frontier;
    std::set<SurfaceMesh::SimplexID<2> > visited;
    bool orientable = true;

    for(auto outer : edgeList)
    {
        if(visited.find(outer) == visited.end()){
            frontier.push_back(outer);

            while(!frontier.empty()){
                auto curr = frontier.front();

                if(visited.find(curr) == visited.end())
                {
                    visited.insert(curr);
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
                                }
                            }
                        }
                        std::vector<SurfaceMesh::SimplexID<2> > tmp;
                        neighbors_up(mesh, curr, std::back_inserter(tmp));
                        for(auto e : tmp){
                            if(std::find(edgeList.begin(), edgeList.end(), e) != edgeList.end())
                                frontier.push_back(e);
                        }
                    }
                    else{
                        std::cerr << "ERROR(computeHoleOrientation): Found an edge"
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

