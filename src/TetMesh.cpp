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

//#include <libraries/casc/include/typetraits.h>
#include "TetMesh.h"





std::unique_ptr<TetMesh> tetgenToTetMesh(tetgenio &tetio){

    std::unique_ptr<TetMesh> mesh(new TetMesh);
    auto &metadata = *mesh.get_simplex_up();

    // Check for higher order cells
    const bool higher_order = tetio.numberofcorners == 10;
   
    metadata.higher_order = higher_order;

    // Copy over vertex data
    for(int i = 0; i < tetio.numberofpoints; ++i){
        double *ptr = tetio.pointlist[i*3];
        int marker = 0;
        mesh->insert({i}, Vertex(ptr[0], ptr[1], ptr[2], marker, false));
    }

    // Copy over tetrahedron data
    for(int i = 0; i < tetio.numberoftetrahedra; ++i){
        // Set material
        int material = 0;
        if(tetio.numberoftetrahedronattributes > 0) 
            material = (int) tetio.tetrahedronattributelist[i * tetio.numberoftetrahedronattributes];

        // Get vertex id's
        //auto idx = i*tetio.numberofcorners;
        int *ptr = tetio.tetrahedronlist[i*tetio.numberofcorners];

        // DETERMINE ORIENTATION TODO...

        mesh.insert<4>({ptr[0], ptr[1], ptr[2], ptr[3]}, Cell(-1, material));

        if (higher_order){
        	// Assign data to edges...
        	// ptr[4], ptr[5], ptr[6], ptr[7], ptr[8], ptr[9]
        }
    }

    // Go over faces and copy data
    for (int i = 0; i < tetio.numberoffacets; ++i){
    	// tetio.trifacelist
    	// tetio.trifacemarkerlist
    }

    // Go over edges and copy data
    // edgelist
    // edgemarkerlist

    return mesh;
}