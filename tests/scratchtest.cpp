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
    gamer::tensor<double, 3, 3> v0;

    std::cout << gamer::Alt(v0) << std::endl;

    std::cout << "EOF" << std::endl;
}