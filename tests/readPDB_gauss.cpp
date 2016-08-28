
#include "gamer.h"
#include <iostream>
#include <regex>

int main(int argc, char *argv[])
{
    char* filename = argv[1];

    float  max_density;
    float *dataset;
    int    xdim, ydim, zdim;
    int    atom_num  = 0;
    ATOM  *atom_list = NULL;
    float  min[3], max[3];
    char   IsXYZR = 0;

    max_density = PDB2Volume(filename, &dataset, &xdim, &ydim, &zdim, min, max,
                             &atom_list, &atom_num, IsXYZR);


    std::cout << atom_num << std::endl;

    for(int i = 0; i < atom_num; ++i)
    {
        std::cout << atom_list[i].radius << std::endl;
    }
/*
    for(int i = 0; i < xdim; ++i)
    {
        for(int j = 0; j < ydim; ++j)
        {
            for(int k = 0; k < zdim; ++k)
            {
                std::cout << dataset[i + xdim*(j + ydim*k)] << std::endl;
            }
        }
    }
*/
    std::cout << xdim << " " << ydim << " " << zdim << std::endl;

    std::cout << max_density << std::endl;
//    SurfaceMesh* mesh = SurfaceMesh::readPDB_gauss(filename, -0.2, 2.5);

//    std::cerr << mesh->numVertices() << std::endl;

    return -1;
}
