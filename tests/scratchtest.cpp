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

    template <std::size_t k>
    int operator()(Complex& F,
                   const std::array<KeyType, k>& new_name,
                   const SimplexSet& merged){
        // std::cout << merged << " -> " << new_name << std::endl;
        return 0;
    }
};

int  main(int argc, char *argv[])
{
    auto mesh = readOFF(argv[1]);
    // auto mesh = readPQR_gauss(argv[1],-0.2, 2.5);
    // writeOFF("test.off", *mesh);

    std::cout << "Mesh read in" << std::endl;

    std::vector<SurfaceMesh*> vec;
    vec.push_back(mesh.get());
    auto tetmesh = makeTetMesh(vec, "test");

    auto edge = *tetmesh->get_level_id<2>().begin();

    edgeCollapse(*tetmesh, edge, 0, Callback<TetMesh>());

    std::cout << "EOF" << std::endl;
}
