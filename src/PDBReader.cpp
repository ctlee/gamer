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

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <ostream>
#include <regex>
#include <cmath>
#include <vector>
#include <limits>

#include "SurfaceMesh.h"
#include "MarchingCube.h"
#include "PDBReader.h"
#include "Vertex.h"


std::unique_ptr<SurfaceMesh> readPDB_distgrid(const std::string& filename, const float radius){
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::vector<AtomType> atomTypes;
    // If readPDB errors return nullptr
    if(!readPDB(filename, std::back_inserter(atomTypes)))
    {
        mesh.reset();
        return mesh;
    }
    std::cout << "Atoms: " << atomTypes.size() << std::endl;
    f3Vector min, max;    
    getMinMax(atomTypes.cbegin(), atomTypes.cend(), min, max, [&radius](const float atomRadius)->float {return DIM_SCALE*(atomRadius + radius);});

    float min_dimension = std::min((max[0] - min[0]), std::min((max[1] - min[1]), (max[2] - min[2])));
    std::cout << "Min Dimension: " << min_dimension << std::endl;

    i3Vector dim;
    f3Vector maxMin = max-min;

    dim = static_cast<i3Vector>((maxMin) + f3Vector({1,1,1})) * DIM_SCALE;

    std::cout << "Dimension: " << dim << std::endl;
    std::cout << "Min:" << min << std::endl;
    std::cout << "Max:" << max << std::endl;

    f3Vector span = (maxMin).ElementwiseDivision(static_cast<f3Vector>(dim) - f3Vector({1,1,1}));
    std::cout << "Delta: " << span << std::endl;

    // Move dataset to positive octant and scale
    for(auto& atom : atomTypes){
        atom.pos = (atom.pos-min).ElementwiseDivision(span);
        atom.radius = (atom.radius + radius)/((span[0] + span[1] + span[2]) / 3.0);
    }

    float* dataset = new float[dim[0]*dim[1]*dim[2]];
    for(int i = 0; i < dim[0]*dim[1]*dim[2]; ++i){
        dataset[i] = -5.0f;
    }
    gridSAS(atomTypes.cbegin(), atomTypes.cend(), dim, dataset);

    // for(int i = 0; i < dim[0]*dim[1]*dim[2]; ++i){
    //     std::cout << dataset[i] << std::endl;
    // }

    std::vector<Vertex> holelist;
    mesh = std::move(marchingCubes(dataset, 5.0f, dim, span, min, 0.0f, std::back_inserter(holelist)));


    delete[] dataset;
    return mesh;
}

std::unique_ptr<SurfaceMesh> readPDB_gauss(const std::string& filename, 
 		    const float blobbyness, 
			float isovalue){
	std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

	std::vector<AtomType> atomTypes;
    // If readPDB errors return nullptr
    if(!readPDB(filename, std::back_inserter(atomTypes)))
    {
        mesh.reset();
        return mesh;
    }

    std::cout << "Atoms: " << atomTypes.size() << std::endl;

    f3Vector min, max;    
    getMinMax(atomTypes.cbegin(), atomTypes.cend(), min, max, 
            [&blobbyness](const float atomRadius)->float{return atomRadius * sqrt(1.0 + log(detail::EPSILON) / blobbyness);});

    float min_dimension = std::min((max[0] - min[0]), std::min((max[1] - min[1]), (max[2] - min[2])));
   
    std::cout << "Min Dimension: " << min_dimension << std::endl;

    i3Vector dim;

    f3Vector maxMin = max-min;

    dim = static_cast<i3Vector>((maxMin) + f3Vector({1,1,1})) * DIM_SCALE;

    std::cout << "Dimension: " << dim << std::endl;
    std::cout << "Min:" << min << std::endl;
    std::cout << "Max:" << max << std::endl;

    f3Vector span = (maxMin).ElementwiseDivision(static_cast<f3Vector>(dim) - f3Vector({1,1,1}));
    std::cout << "Delta: " << span << std::endl;

    float* dataset = new float[dim[0]*dim[1]*dim[2]]();

    // Bring atoms to +++ quadrant
    // for(auto& atom : atomTypes){
    //     atom.pos = (atom.pos-min).ElementwiseDivision(span);
    // }

    std::cout << "Begin blurring coordinates" << std::endl;
    blurAtoms(atomTypes.cbegin(), atomTypes.cend(), dataset, min, maxMin, dim, blobbyness);
    std::cout << "Done blurring coords" << std::endl;

    float minval, maxval;
    minval = std::numeric_limits<float>::infinity();
    maxval = -std::numeric_limits<float>::infinity();

    for (int i = 0; i < dim[2] * dim[1] * dim[0]; ++i)
    {
        float cval = dataset[i];
        if (cval < minval){
            minval = cval;
        }
        if (cval > maxval){
            maxval = cval;
        }
    }
    std::cout << "Min Density: " << minval << ", Max Density: " << maxval << std::endl;

    float data_isoval = 0.44 * maxval;

    if (data_isoval < isovalue)
    {
        isovalue = data_isoval;
    }
    std::cout << "Isovalue: " << isovalue << std::endl;
 
    std::vector<Vertex> holelist;
    mesh = std::move(marchingCubes(dataset, maxval, dim, span, min, isovalue, std::back_inserter(holelist)));
    delete[] dataset;
    // TODO: (0) What to do with holelist...
    return mesh;
}