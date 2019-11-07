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

#include "gamer/EigenDiagonalization.h"
#include "gamer/OsculatingJets.h"
#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{
std::tuple<REAL*, REAL*, REAL*, REAL*, std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType>>
curvatureViaMDSB(const SurfaceMesh& mesh){
    std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType> sigma;

    REAL *Amix = new REAL[mesh.size<1>()];
    REAL *kg = new REAL[mesh.size<1>()];
    REAL *kh = new REAL[mesh.size<1>()];
    REAL *k1 = new REAL[mesh.size<1>()];
    REAL *k2 = new REAL[mesh.size<1>()];
    Vector *Kh = new Vector[mesh.size<1>()];
    Vector *normals = new Vector[mesh.size<1>()];

    // Map VertexIDs to indices
    std::size_t i = 0;
    for (const auto vertexID : mesh.get_level_id<1>()) {
        Amix[i] = 0;
        kg[i] = 0;
        kh[i] = 0;
        k1[i] = 0;
        k2[i] = 0;
        normals[i] = getNormal(mesh, vertexID);
        sigma[vertexID.indices()[0]] = i++;
    }

    for (const auto faceID : mesh.get_level_id<3>()) {
        auto indices = faceID.indices(); // Face vertex indices
        std::array<Vertex, 3> vertices;  // Vertex data for indices
        std::array<typename SurfaceMesh::KeyType, 2> keysDown; // Tmp to store
                                                               // keys

        // Fill vertices in order corresponding to indices
        for (std::size_t skip = 0; skip < 3; ++skip) {
            std::size_t j = 0;
            std::size_t k = 0;
            while (k < 2) {
                if (j == skip) ++j;
                else{
                    keysDown[k++] = indices[j++];
                }
            }
            vertices[skip] = *mesh.get_simplex_down(faceID, keysDown);
        }

        // TODO: (15) This section computes the same distances a bunch of time
        std::array<REAL, 3> dist;
        dist[0] = distance(vertices[0], vertices[1]);
        dist[1] = distance(vertices[1], vertices[2]);
        dist[2] = distance(vertices[0], vertices[2]);
        std::sort(dist.begin(), dist.end());

        // Check if the triangle is obtuse...
        bool obtuse = dist[0]*dist[0] + dist[1]*dist[1] < dist[2]*dist[2];

        REAL t_area = 0; // Area of the face
        // Populate t_area if obtuse
        if (obtuse) t_area = getArea(vertices[0], vertices[1], vertices[2]);

        // List of indices to rotate
        std::array<std::size_t, 3> idxmap = {0, 1, 2};
        for (std::size_t i = 0; i < 3; ++i) {
            // idxmap[0] is the current vertex
            REAL ang = angle(vertices[idxmap[2]], vertices[idxmap[0]], vertices[idxmap[1]]);

            std::size_t i0 = sigma[indices[idxmap[0]]];
            std::size_t i1 = sigma[indices[idxmap[1]]];
            std::size_t i2 = sigma[indices[idxmap[2]]];

            // Add angle to Gaussian Curvature
            kg[i0] += ang;

            // Vectors of other edges
            Vector v1 = vertices[idxmap[1]] - vertices[idxmap[2]];
            Vector v2 = vertices[idxmap[2]] - vertices[idxmap[1]];

            REAL cot = 1.0/tan(ang);

            if (obtuse) {
                if (ang > M_PI/2.0) {
                    Amix[i0] += t_area/2.0;
                }
                else{
                    Amix[i0] += t_area/4.0;
                }
            }
            else{
                REAL tmp = cot/8.0;
                REAL lenSq = length(v1);
                lenSq *= lenSq;
                Amix[i1] += tmp*lenSq;
                Amix[i2] += tmp*lenSq;
            }

            // Add value to Mean Curvature
            Kh[i1] += cot*v1;
            Kh[i2] += cot*v2;

            std::rotate(idxmap.begin(), idxmap.begin() + 1, idxmap.end());
        }
    }

    for (std::size_t i = 0; i < mesh.size<1>(); ++i) {
        Kh[i] = Kh[i]/(2.0*Amix[i]);
        kh[i] = std::copysign(length(Kh[i])/2.0, -dot(Kh[i], normals[i]));
        kg[i] = (2.0*M_PI - kg[i])/Amix[i];

        REAL kh2 = kh[i]*kh[i];
        REAL tmp = kh2 < kg[i] ? 0 : std::sqrt(kh2-kg[i]);

        k1[i] = kh[i] + tmp;
        k2[i] = kh[i] - tmp;
    }

    delete[] Amix;
    delete[] Kh;
    delete[] normals;
    return std::make_tuple(kh, kg, k1, k2, sigma);
}


std::tuple<REAL*, REAL*, REAL*, REAL*, std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType>>
curvatureViaJets(const SurfaceMesh& mesh, std::size_t dJet, std::size_t dPrime){
    std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType> sigma;

    REAL *kg = new REAL[mesh.size<1>()];
    REAL *kh = new REAL[mesh.size<1>()];
    REAL *k1 = new REAL[mesh.size<1>()];
    REAL *k2 = new REAL[mesh.size<1>()];

    int min_nb_points = (dJet + 1) * (dJet + 2) / 2;

    // Map VertexIDs to indices
    std::size_t i = 0;
    for (const auto vertexID : mesh.get_level_id<1>()) {
        std::vector<SurfaceMesh::SimplexID<1>> nbors;
        nbors.push_back(vertexID);
        surfacemesh_detail::vertexGrabber(mesh, min_nb_points-1, nbors, vertexID);

        if (nbors.size() < min_nb_points) {
            std::cerr << "Not enough pts (have: " << nbors.size() << ", need: "      << min_nb_points << ") for fitting this vertex: "
                      << vertexID << std::endl;
            continue;
        }

        Monge_via_jet_fitting monge_fit;
        auto mongeForm = monge_fit(nbors.begin(), nbors.end(), dJet, dPrime);
        mongeForm.comply_wrt_given_normal(getNormal(mesh, vertexID));

        REAL tk1 = mongeForm.principal_curvatures(0);
        REAL tk2 = mongeForm.principal_curvatures(1);

        k1[i] = tk1;
        k2[i] = tk2;

        kg[i] = tk1*tk2;
        kh[i] = (tk1+tk2)/2.;
        sigma[vertexID.indices()[0]] = i++;
    }
    return std::make_tuple(kh, kg, k1, k2, sigma);
}
}