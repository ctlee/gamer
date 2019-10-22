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
#include <array>
#include <memory>

#include "gamer/gamer"
// #include <casc/casc>


int  main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    auto mesh = gamer::readOFF(argv[1]);

    casc::compute_orientation(*mesh);
    if (gamer::getVolume(*mesh) < 0) {
        gamer::flipNormals(*mesh);
    }

    gamer::curvatureViaJets(*mesh, 2, 2);

    std::cout << "EOF" << std::endl;
}