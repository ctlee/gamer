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
#include "SurfaceMesh.h"
#include "TetMesh.h"
#include "Vertex.h"
#include "tensor.h"
#include "PDBReader.h"

#include <libraries/casc/include/CASCFunctions.h>
#include <libraries/casc/include/SimplexSet.h>
#include <libraries/casc/include/SimplexMap.h>
#include <libraries/casc/include/decimate.h>
#include <libraries/casc/include/typetraits.h>

int  main(int argc, char *argv[])
{
    auto mesh = readOFF(argv[1]);

    fillHoles(*mesh);

    writeOFF("test.off", *mesh);

    auto vol = getVolume(*mesh);
    if (vol < 0){
        flipNormals(*mesh);
    }
    // vector of vectors >.<
    // std::vector<std::vector<SurfaceMesh::SimplexID<2>>> test;
    // surfacemesh_detail::findHoles(*mesh, test);

    // std::cout << "Number of holes: " << test.size()  << std::endl << std::endl;

    // int i = 0;
    // for(auto v : test){
    //     std::cout << "Begin ring: ";
    //     for(auto s : v){
    //         ++i;
    //         std::cout << s << " ";
    //     }
    //     std::cout << std::endl;

    //     std::cout << "Sorted Vertices: ";
    //     std::vector<SurfaceMesh::SimplexID<1>> foo;
    //     surfacemesh_detail::edgeRingToVertices(*mesh, v, std::back_inserter(foo));
    //     for(auto vert : foo){
    //         std::cout << vert << " ";
    //     }
    //     std::cout << std::endl << std::endl;
    // }

    // std::cout << "Number of edges: " << i << std::endl;
    std::cout << "EOF" << std::endl;
}
