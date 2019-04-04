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

    auto vol = getVolume(*mesh);
    if (vol < 0){
        flipNormals(*mesh);
    }

    auto s = *mesh->get_level_id<1>().begin();

    auto curve = getMeanCurvature(*mesh, s);

    std::cout << "Curvature = " << std::sqrt(curve|curve)/2 << std::endl;

    std::cout << "EOF" << std::endl;
}
