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
    // auto mesh = readPDB_molsurf("2jho.pdb");

    // smoothMesh(*mesh, 10, true, false);

    // writeOFF("2jho.off", *mesh);

    auto mesh = readOFF(argv[1]);

    auto vol = getVolume(*mesh);
    if (vol < 0){
        flipNormals(*mesh);
    }

    double max = 0;
    for (auto s : mesh->get_level_id<1>()){
        auto curve = getMeanCurvature(*mesh, s);
        std::cout << s << ": curvature = " << curve << std::endl;
        if (curve > max)
            max = curve;
    }
    std::cout << "Max: " << max << std::endl;

    std::cout << "EOF" << std::endl;
}
