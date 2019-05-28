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

#include "gamer/gamer"

#include <casc/casc>
#include <libraries/casc/include/typetraits.h>

int  main(int argc, char *argv[])
{
    auto v = Vertex(1,2,3);
    std::cout << type_name<decltype(&Vertex::position)>() << std::endl;

    // auto mesh = readPDB_molsurf("5kp9_bioassembly1.pdb");

    // smoothMesh(*mesh, 10, true, false);

    // writeOFF("5kp9_bioassembly1.off", *mesh);

    // auto mesh = readOFF(argv[1]);

    // auto vol = getVolume(*mesh);
    // if (vol < 0){
    //     flipNormals(*mesh);
    // }


    // double max = 0;
    // for (auto s : mesh->get_level_id<1>()){
    //     auto curve = getGaussianCurvature(*mesh, s);
    //     std::cout << s << ": curvature = " << curve << std::endl;
    //     if (curve > max)
    //         max = curve;
    // }
    // std::cout << "Max: " << max << std::endl;
    //
    // auto mesh = readDolfin(argv[1]);

    // writeDolfin("test.xml", *mesh);
    std::cout << "EOF" << std::endl;
}
