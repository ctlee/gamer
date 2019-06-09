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

int  main(int argc, char *argv[])
{
    auto mesh = readPDB_molsurf("2jho.pdb");
    compute_orientation(*mesh);
    smoothMesh(*mesh, 10, true, false);

    auto vol = getVolume(*mesh);
    if (vol < 0){
        flipNormals(*mesh);
    }

    auto lid = mesh->get_level_id<1>();
    for(auto vid = lid.begin(); vid != lid.end();){
        std::cout << *vid << std::endl;
        if((*(*vid))[1] > 1){
            auto value = *vid;
            ++vid;
            mesh->remove(value);
            continue;
        }
        ++vid;
    }

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
