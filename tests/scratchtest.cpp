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
#include <Eigen/core>

// #include <casc/casc>


int  main(int argc, char *argv[])
{

    // std::array<REAL, 6> covariance = {{ 1, 2, 3, 4, 5, 6 }};
    Eigen::Matrix<double,3,3> m;
    m << 1,2,3,4,5,6,7,8,9;

    for (int i =0; i < 9; ++i){
      std::cout << m(i) << std::endl;
    }

    std::cout << m(1,2) << std::endl;

    std::cout << "EOF" << std::endl;
}