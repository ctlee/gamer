/*
 * ***************************************************************************
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 1994-- Michael Holst and Zeyun Yu
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
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     PDB2Volume.C    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Convert a PDB into a 3D Volume
 *
 * Source:   Part of this file was adapted from the PDBParser developed in
 *           Chandrajit Bajaj's group at The University of Texas at Austin.
 * ***************************************************************************
 */


#include "PDB2Volume.h"
#include <limits>
#include <regex>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "ReadPDB.h"

#define MAX_STRING       256

void write_rawiv_float(FILE *,
                       float *,
                       int *,
                       float *,
                       float *);



/*
 * ***************************************************************************
 * Routine:  evalDensity    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the Gaussian function value using the blobbyness
 * ***************************************************************************
 */
float evalDensity(const ATOM& atom, float pnt[3], double maxRadius)
{
    double expval;

    double r = (atom.x - pnt[0]) * (atom.x - pnt[0]) +
               (atom.y - pnt[1]) * (atom.y - pnt[1]) +
               (atom.z - pnt[2]) * (atom.z - pnt[2]);
    double r0 = atom.radius; r0 *= r0;

    // expval = BLOBBYNESS*(r/r0 - 1.0);
    expval = BLOBBYNESS * (r - r0);

    // truncated gaussian
    if (sqrt(r) > maxRadius)
    {
        return 0.0;
    }

    return (float)(exp(expval));
}



/*
 * ***************************************************************************
 * Routine:  blurAtoms    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Blur all atoms with Gaussian functions
 *           Each atom spreads in a region of certain radius
 * ***************************************************************************
 */
void blurAtoms(std::vector<ATOM>::const_iterator begin,
               std::vector<ATOM>::const_iterator end,
               float min[3],
               float max[3],
               float *dataset,
               int dim[3])
{
    float  orig[3], span[3];
    int    i, j, k;
    int    m;
    int    xdim, ydim, zdim;
    int    amax[3], amin[3];
    double c[3], maxRad;


    xdim = dim[0];
    ydim = dim[1];
    zdim = dim[2];

    for (k = 0; k < zdim; k++)
    {
        for (j = 0; j < ydim; j++)
        {
            for (i = 0; i < xdim; i++)
            {
                dataset[IndexVect(i, j, k)] = 0;
            }
        }
    }

    orig[0] = min[0];
    orig[1] = min[1];
    orig[2] = min[2];
    span[0] = (max[0] - min[0]) / (float)(xdim - 1);
    span[1] = (max[1] - min[1]) / (float)(ydim - 1);
    span[2] = (max[2] - min[2]) / (float)(zdim - 1);

    m = 0;
    for (auto curr = begin; curr != end; ++curr)
    {
        maxRad = curr->radius * sqrt(1.0 + log(detail::EPSILON) / (2.0 * BLOBBYNESS));

        // compute the dataset coordinates of the atom's center
        c[0] = (curr->x - orig[0]) / span[0];
        c[0] = ((c[0] - floor(c[0])) >= 0.5) ? ceil(c[0]) : floor(c[0]);
        c[1] = (curr->y - orig[1]) / span[1];
        c[1] = ((c[1] - floor(c[1])) >= 0.5) ? ceil(c[1]) : floor(c[1]);
        c[2] = (curr->z - orig[2]) / span[2];
        c[2] = ((c[2] - floor(c[2])) >= 0.5) ? ceil(c[2]) : floor(c[2]);

        // then compute the bounding box of the atom (maxRad^3)
        for (j = 0; j < 3; j++)
        {
            int tmp;
            tmp     = (int)(c[j] - (maxRad / span[j]) - 1);
            tmp     = (tmp < 0) ? 0 : tmp;
            amin[j] = tmp;
            tmp     = (int)(c[j] + (maxRad / span[j]) + 1);
            tmp     = (tmp > (dim[j] - 1)) ? (dim[j] - 1) : tmp;
            amax[j] = tmp;
        }

        // begin blurring kernel
        for (k = amin[2]; k <= amax[2]; k++)
        {
            for (j = amin[1]; j <= amax[1]; j++)
            {
                for (i = amin[0]; i <= amax[0]; i++)
                {
                    float pnt[3], density;

                    pnt[0] = orig[0] + i * span[0];
                    pnt[1] = orig[1] + j * span[1];
                    pnt[2] = orig[2] + k * span[2];

                    density                      = evalDensity(*curr, pnt, maxRad);
                    dataset[IndexVect(i, j, k)] += density;
                }
            }
        }
/*
        if ((((m + 1) % 20) == 0) || ((m + 1) == size))
        {
            printf("%2.2f%% done (%08d)\r", 100.0 * (m + 1) / (float)size, m + 1);
            fflush(stdout);
        }
        */
    }
    printf("\n"); fflush(stdout);
}



/*
 * ***************************************************************************
 * Routine:  PDB2Volume    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Convert a PDB/PQR/XYZR formats into a 3D Volume
 *
 * Notes:    Each atom is treated as a Gaussian function defined by
 *           - center and radius: (available in the PDB/PQR files)
 *           - blobbyness: the decay rate of the Gaussian function
 * ***************************************************************************
 */
float PDB2Volume(std::string filename, float **data, int *xd, int *yd, int *zd,
                 float min[3], float max[3], std::vector<ATOM>& atom_list)
{
    int    dim[3];
    float  min_dimension;
    float *dataset;
    PDBelementInformation eInfo;

    readAtomFile(filename, std::back_inserter(atom_list));

    getMinMax(atom_list.cbegin(), atom_list.cend(), min, max);

    min_dimension = std::min((max[0] - min[0]), std::min((max[1] - min[1]), (max[2] - min[2])));

    if (min_dimension < 64.0f)
    {
        min_dimension = 64.0f / min_dimension;
        dim[0]        = (int)((max[0] - min[0]) * min_dimension) + 1;
        dim[1]        = (int)((max[1] - min[1]) * min_dimension) + 1;
        dim[2]        = (int)((max[2] - min[2]) * min_dimension) + 1;
    }
    else
    {
        dim[0] = (int)(max[0] - min[0]) + 1;
        dim[1] = (int)(max[1] - min[1]) + 1;
        dim[2] = (int)(max[2] - min[2]) + 1;
    }
    dim[0] = (int)(dim[0] * DIM_SCALE);
    dim[1] = (int)(dim[1] * DIM_SCALE);
    dim[2] = (int)(dim[2] * DIM_SCALE);

    printf("dimension: %d X %d X %d\n", dim[0], dim[1], dim[2]);

    dataset = (float *)malloc(sizeof(float) * dim[0] * dim[1] * dim[2]);
    blurAtoms(atom_list.begin(), atom_list.end(), min, max, dataset, dim);
    printf("min[3]: %f %f %f \n",  min[0], min[1], min[2]);
    printf("max[3]: %f %f %f \n",  max[0], max[1], max[2]);
    printf("span[3]: %f %f %f \n", (max[0] - min[0]) / (float)(dim[0] - 1),
           (max[1] - min[1]) / (float)(dim[1] - 1), (max[2] - min[2]) / (float)(dim[2] - 1));

    float minval, maxval;
    minval = std::numeric_limits<float>::infinity();
    maxval = -std::numeric_limits<float>::infinity();

    for (int i = 0; i < dim[2] * dim[1] * dim[0]; ++i)
    {
        float cval = dataset[i];

        if (cval < minval)
        {
            minval = cval;
        }

        if (cval > maxval)
        {
            maxval = cval;
        }
    }

    printf("min_density: %f   max_density: %f \n", minval, maxval);


    /* write the volume to disk
       if ((fp=fopen("test.rawiv", "w"))==NULL){
       printf("write error...\n");
       exit(0);
       };
       write_rawiv_float(fp,dataset,dim,min,max);
       exit(0);
     */


    *data     = dataset;
    *xd       = dim[0];
    *yd       = dim[1];
    *zd       = dim[2];
    return maxval;
}




/*
 * ***************************************************************************
 * Routine:  SurfaceMeshOld::readPDB_gauss    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Convert a PDB/PQR/XYZR formats into a SurfaceMeshOld using Gaussian
 *           scalar function
 *
 * Notes:    Each atom is treated as a Gaussian function defined by
 *           - center and radius: (available in the PDB/PQR files)
 *           - blobbyness: the decay rate of the Gaussian function
 * ***************************************************************************
 */
std::tuple<SurfaceMeshOld*,SurfaceMesh_ASC*> SurfaceMeshOld::readPDB_gauss(const char *filename, float blobbyness,
                                         float iso_value)
{
    float  max_density;
    float *dataset;
    float  data_iso_val;
    time_t t1, t2;
    int    xdim, ydim, zdim;
    float  min[3], max[3], span[3];
    SPNT *holelist;
    std::vector<ATOM> atom_list;

    printf("\nbegin blurring PDB/PQR coordinates ... \n");
    time(&t1);
    max_density = PDB2Volume(filename, &dataset, &xdim, &ydim, &zdim, min, max,
                             atom_list);
    (void)time(&t2);
    printf("time to generate volume from PDB: %d seconds. \n\n", (int)(t2 - t1));
    span[0] = (max[0] - min[0]) / (float)(xdim - 1);
    span[1] = (max[1] - min[1]) / (float)(ydim - 1);
    span[2] = (max[2] - min[2]) / (float)(zdim - 1);

    printf("begin extracting isosurfaces ... \n");
    (void)time(&t1);
    data_iso_val = 0.44 * max_density;

    if (data_iso_val < iso_value)
    {
        iso_value = data_iso_val;
    }

    printf("isovalue: %f \n", iso_value);

    auto meshes = SurfaceMeshOld::marchingCube(xdim, ydim, zdim, dataset, iso_value, &holelist);
    auto surfmesh = std::get<0>(meshes);
    auto pF = std::get<1>(meshes);

    time(&t2);
    std::cout << "vertices: " << pF->size<1>() << ", faces: " << pF->size<3>() << std::endl;
    printf("time to extract isosurface: %d seconds. \n\n", (int)(t2 - t1));

    free(dataset);

    // convert from pixel to angstrom
    for (int j = 0; j < surfmesh->num_vertices; j++)
    {
        surfmesh->vertex[j].x = surfmesh->vertex[j].x * span[0] + min[0];
        surfmesh->vertex[j].y = surfmesh->vertex[j].y * span[1] + min[1];
        surfmesh->vertex[j].z = surfmesh->vertex[j].z * span[2] + min[2];
    }

    for (auto& curr : pF->get_level<1>())
    {
        curr.position.get(0) = curr.position.get(0) * span[0] + min[0];
        curr.position.get(1) = curr.position.get(1) * span[1] + min[1];
        curr.position.get(2) = curr.position.get(2) * span[2] + min[2];
    }

    // Flip normals so they now points outwards
    surfmesh->flipNormals();

    // Return generated SurfaceMeshOld
    return meshes;
}
