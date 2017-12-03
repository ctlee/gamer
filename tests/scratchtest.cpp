/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "stringutil.h"
#include "SurfaceMeshOld.h"
#include "SurfaceMesh.h"
#include "Vertex.h"
#include "tensor.h"

#include <libraries/casc/include/CASCFunctions.h>
#include <libraries/casc/include/SimplexSet.h>
#include <libraries/casc/include/SimplexMap.h>
#include <libraries/casc/include/decimate.h>
#include <libraries/casc/include/typetraits.h>

int  main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    auto mesh = readOFF(argv[1]);
    if(mesh == nullptr){
        std::cout << "Something bad happened..." << std::endl;
        exit(1);
    }
    std::cout << "Done reading" << std::endl;

    compute_orientation(*mesh);

    double volume = getVolume(*mesh);
    std::cout << "Volume: " << volume << std::endl;

    //coarseIT(*mesh, 0.0160, 1, 0);
    coarseIT(*mesh, 2.5, 0, 10);
    writeOFF("test.off", *mesh);



    // auto oldmesh = SurfaceMeshOld::readOFF(argv[1]);
    // //oldmesh->coarse(0.016,1,0,-1);
    // oldmesh->coarse(2.5,0,10,-1);
    // 
    std::cout << "EOF" << std::endl;
}
