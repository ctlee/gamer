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


template <typename Complex>
struct Callback
{
    using SimplexSet = typename casc::SimplexSet<Complex>;
    using KeyType = typename Complex::KeyType;

    Face operator()(Complex& F,
            const std::array<KeyType, 3>& new_name,
            const SimplexSet& merged){
      std::cout << merged << " -> " << casc::to_string(new_name) << std::endl;
      return **merged.template cbegin<3>();
    }

    Edge operator()(Complex& F,
            const std::array<KeyType, 2>& new_name,
            const SimplexSet& merged){
      std::cout << merged << " -> " << casc::to_string(new_name) << std::endl;
      return **merged.template cbegin<2>();
    }

    Vertex operator()(Complex& F,
            const std::array<KeyType, 1>& new_name,
            const SimplexSet& merged){
      std::cout << merged << " -> " << casc::to_string(new_name) << std::endl;
      return **merged.template cbegin<1>();
    }

    // template <std::size_t k>
    // Complex::NodeDataTypes<k> operator()(Complex& F,
    //         const std::array<KeyType, k>& new_name,
    //         const SimplexSet& merged){
    //     // std::cout << merged << " -> " << new_name << std::endl;
    //     // return 0;
    // }
};

int  main(int argc, char *argv[])
{
    auto mesh = readOFF(argv[1]);
    // auto mesh = readPQR_gauss(argv[1],-0.2, 2.5);
    // writeOFF("test.off", *mesh);

    std::cout << "Mesh read in" << std::endl;

    auto edge = *mesh->get_level_id<2>().begin();

    casc::decimate(*mesh, edge, Callback<SurfaceMesh>());

    std::cout << "EOF" << std::endl;
}
