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

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <ostream>
#include <regex>
#include <cmath>
#include <vector>
#include <limits>

#include "gamer/SurfaceMesh.h"
#include "gamer/MarchingCube.h"
#include "gamer/PDBReader.h"
#include "gamer/Vertex.h"

/// Namespace for all things gamer
namespace gamer
{

std::unique_ptr<SurfaceMesh> readPDB_distgrid(const std::string &filename, const float radius)
{
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::vector<Atom>            atoms;
    // If readPDB errors return nullptr
    if (!readPDB(filename, std::back_inserter(atoms)))
    {
        mesh.reset();
        return mesh;
    }
    std::cout << "Atoms: " << atoms.size() << std::endl;
    Vector3f min, max;
    getMinMax(atoms.cbegin(), atoms.cend(), min, max, [&radius](const float atomRadius) -> float {
            return DIM_SCALE*(atomRadius + radius);
        });

    float min_dimension = std::min((max[0] - min[0]), std::min((max[1] - min[1]), (max[2] - min[2])));
    std::cout << "Min Dimension: " << min_dimension << std::endl;

    Vector3i dim;
    Vector3f maxMin = max-min;

    dim = static_cast<Vector3i>((maxMin) + Vector3f({1, 1, 1})) * DIM_SCALE;

    std::cout << "Dimension: " << dim << std::endl;
    std::cout << "Min:" << min << std::endl;
    std::cout << "Max:" << max << std::endl;

    Vector3f span = (maxMin).ElementwiseDivision(static_cast<Vector3f>(dim) - Vector3f({1, 1, 1}));
    std::cout << "Delta: " << span << std::endl;

    // Move dataset to positive octant and scale
    for (auto &atom : atoms)
    {
        atom.pos = (atom.pos-min).ElementwiseDivision(span);
        atom.radius = (atom.radius + radius)/((span[0] + span[1] + span[2]) / 3.0);
    }

    float* dataset = new float[dim[0]*dim[1]*dim[2]];
    for (int i = 0; i < dim[0]*dim[1]*dim[2]; ++i)
    {
        dataset[i] = -5.0f;
    }
    gridSAS(atoms.cbegin(), atoms.cend(), dim, dataset);

    std::vector<Vertex>          holelist;
    std::unique_ptr<SurfaceMesh> SASmesh = std::move(marchingCubes(dataset, 5.0f, dim, span, 0.0f, std::back_inserter(holelist)));

    for (auto curr = atoms.cbegin(); curr != atoms.cend(); ++curr)
    {
        Vector3f pos = curr->pos;
        // compute the dataset coordinates of the atom's center
        Vector3i c;
        std::transform(pos.begin(), pos.end(), c.begin(), [](float v) -> int {
                return round(v);
            });
    }

    // Reset dataset
    for (int i = 0; i < dim[0]*dim[1]*dim[2]; ++i)
    {
        dataset[i] = -5.0f;
    }

    auto SASverts = SASmesh->get_level<1>();
    gridSES(SASverts.begin(), SASverts.end(), dim, dataset, radius);
    // for(int i = 0; i < dim[0]*dim[1]*dim[2]; ++i){
    //     std::cout << dataset[i] << std::endl;
    // }

    mesh = std::move(marchingCubes(dataset, 5.0f, dim, span, 0.0f, std::back_inserter(holelist)));
    delete[] dataset;

    // Translate back to the original position from the positive octant
    for (auto &v : mesh->get_level<1>())
    {
        v += min;
    }
    return mesh;
}

/**
 * @brief      Reads a pdb gauss.
 *
 * @param[in]  filename    The filename
 * @param[in]  blobbyness  The blobbyness
 * @param[in]  isovalue    The isovalue
 *
 * @return     { description_of_the_return_value }
 */
std::unique_ptr<SurfaceMesh> readPDB_gauss(const std::string &filename,
                                           const float        blobbyness,
                                           float              isovalue)
{
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::vector<Atom>            atoms;
    // If readPDB errors return nullptr
    if (!readPDB(filename, std::back_inserter(atoms)))
    {
        mesh.reset();
        return mesh;
    }

    std::cout << "Atoms: " << atoms.size() << std::endl;

    Vector3f min, max;
    getMinMax(atoms.cbegin(), atoms.cend(), min, max,
              [&blobbyness](const float atomRadius) -> float{
            return atomRadius * sqrt(1.0 + log(pdbreader_detail::EPSILON) / blobbyness);
        });

    float min_dimension = std::min((max[0] - min[0]), std::min((max[1] - min[1]), (max[2] - min[2])));

    std::cout << "Min Dimension: " << min_dimension << std::endl;

    Vector3i dim;

    Vector3f maxMin = max-min;

    dim = static_cast<Vector3i>((maxMin) + Vector3f({1, 1, 1})) * DIM_SCALE;

    std::cout << "Dimension: " << dim << std::endl;
    std::cout << "Min:" << min << std::endl;
    std::cout << "Max:" << max << std::endl;

    Vector3f span = (maxMin).ElementwiseDivision(static_cast<Vector3f>(dim) - Vector3f({1, 1, 1}));
    std::cout << "Delta: " << span << std::endl;

    float* dataset = new float[dim[0]*dim[1]*dim[2]]();

    // Bring atoms to +++ quadrant
    // for(auto& atom : atoms){
    //     atom.pos = (atom.pos-min).ElementwiseDivision(span);
    // }

    std::cout << "Begin blurring coordinates" << std::endl;
    blurAtoms(atoms.cbegin(), atoms.cend(), dataset, min, maxMin, dim, blobbyness);
    std::cout << "Done blurring coords" << std::endl;

    float minval = std::numeric_limits<float>::infinity();
    float maxval = -std::numeric_limits<float>::infinity();

    for (int i = 0; i < dim[2] * dim[1] * dim[0]; ++i)
    {
        float cval = dataset[i];
        if (cval < minval){
            minval = cval;
        }
        if (cval > maxval)
        {
            maxval = cval;
        }
    }
    // std::cout << "Min Density: " << minval << ", Max Density: " << maxval << std::endl;
    float data_isoval = 0.44 * maxval; // Override the user's isovalue... is
                                       // this a good idea?
    if (data_isoval < isovalue)
    {
        isovalue = data_isoval;
    }
    std::cout << "Isovalue: " << isovalue << std::endl;

    std::vector<Vertex> holelist;
    mesh = std::move(marchingCubes(dataset, maxval, dim, span, isovalue, std::back_inserter(holelist)));
    delete[] dataset;

    // Translate back to the original position from the positive octant
    for (auto &v : mesh->get_level<1>())
    {
        v += min;
    }
    // TODO: (50) What to do with holelist...
    return mesh;
}

std::unique_ptr<SurfaceMesh> readPQR_gauss(const std::string &filename,
                                           const float        blobbyness,
                                           float              isovalue)
{
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::vector<Atom>            atoms;
    // If readPDB errors return nullptr
    if (!readPQR(filename, std::back_inserter(atoms)))
    {
        mesh.reset();
        return mesh;
    }

    std::cout << "Atoms: " << atoms.size() << std::endl;

    Vector3f min, max;
    getMinMax(atoms.cbegin(), atoms.cend(), min, max,
              [&blobbyness](const float atomRadius) -> float{
            return atomRadius * sqrt(1.0 + log(pdbreader_detail::EPSILON) / blobbyness);
        });

    float min_dimension = std::min((max[0] - min[0]), std::min((max[1] - min[1]), (max[2] - min[2])));

    std::cout << "Min Dimension: " << min_dimension << std::endl;

    Vector3i dim;

    Vector3f maxMin = max-min;

    dim = static_cast<Vector3i>((maxMin) + Vector3f({1, 1, 1})) * DIM_SCALE;

    std::cout << "Dimension: " << dim << std::endl;
    std::cout << "Min:" << min << std::endl;
    std::cout << "Max:" << max << std::endl;

    Vector3f span = (maxMin).ElementwiseDivision(static_cast<Vector3f>(dim) - Vector3f({1, 1, 1}));
    std::cout << "Delta: " << span << std::endl;

    float* dataset = new float[dim[0]*dim[1]*dim[2]]();

    // Bring atoms to +++ quadrant
    // for(auto& atom : atoms){
    //     atom.pos = (atom.pos-min).ElementwiseDivision(span);
    // }

    std::cout << "Begin blurring coordinates" << std::endl;
    blurAtoms(atoms.cbegin(), atoms.cend(), dataset, min, maxMin, dim, blobbyness);
    std::cout << "Done blurring coords" << std::endl;

    // float minval;
    float maxval;
    // minval = std::numeric_limits<float>::infinity();
    maxval = -std::numeric_limits<float>::infinity();

    for (int i = 0; i < dim[2] * dim[1] * dim[0]; ++i)
    {
        float cval = dataset[i];
        // if (cval < minval){
        //     minval = cval;
        // }
        if (cval > maxval)
        {
            maxval = cval;
        }
    }
    // std::cout << "Min Density: " << minval << ", Max Density: " << maxval <<
    // std::endl;
    float data_isoval = 0.44 * maxval; // Override the user's isovalue... is
                                       // this a good idea?
    if (data_isoval < isovalue)
    {
        isovalue = data_isoval;
    }
    std::cout << "Isovalue: " << isovalue << std::endl;

    std::vector<Vertex> holelist;
    mesh = std::move(marchingCubes(dataset, maxval, dim, span, isovalue, std::back_inserter(holelist)));
    delete[] dataset;

    // Translate back to the original position from the positive octant
    for (auto &v : mesh->get_level<1>())
    {
        v += min;
    }
    // TODO: (50) What to do with holelist...
    return mesh;
}

} // end namespace gamer
