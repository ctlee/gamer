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
#include "SurfaceMeshOld.h"
#include "SurfaceMesh.h"
#include "TetMesh.h"
#include "Vertex.h"
#include "tensor.h"

#include <libraries/casc/include/CASCFunctions.h>
#include <libraries/casc/include/SimplexSet.h>
#include <libraries/casc/include/SimplexMap.h>
#include <libraries/casc/include/decimate.h>
#include <libraries/casc/include/typetraits.h>

int  main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    auto mesh = readOFF(argv[1]);
    if(mesh == nullptr){
        std::cout << "Something bad happened..." << std::endl;
        exit(1);
    }
    std::cout << "Done reading" << std::endl;

    compute_orientation(*mesh);
    if(getVolume(*mesh) < 0){
        for(auto &face : mesh->get_level<3>())
            face.orientation *= -1;
    }

    // coarseIT(*mesh, 0.016, 1, 0);
    // coarseIT(*mesh, 2.5, 0, 10);

    // smoothMesh(*mesh, 15, 165, 5, true);
    writeOFF("tmp.off", *mesh);

    // auto &global = *mesh->get_simplex_up();
    // global.closed = true;
    // global.ishole = true;

    // auto box = sphere(5);    

    // auto &global2 = *box->get_simplex_up();
    // global2.closed = true;
    // global2.ishole = false;

    // Vector center;
    // double radius;
    // std::tie(center, radius) = getCenterRadius(*mesh);

    // scale(*box, radius*10);
    // writeOFF("box.off", *box);

    // std::vector<std::unique_ptr<SurfaceMesh>> meshes;

    // std::cout << "Volume of mesh: " << getVolume(*mesh) << std::endl;
    // std::cout << "Volume of box: " << getVolume(*box) << std::endl;

    // meshes.push_back(std::move(mesh));
    // meshes.push_back(std::move(box));

    // std::cout << "\n\n\nBEGINNING TETRAHEDRALIZATION..." << std::endl;
    // auto temesh = makeTetMesh(meshes, "pq1.4zYAV");

    // writeOFF("test.off", *mesh);
    // auto oldmesh = SurfaceMeshOld::readOFF(argv[1]);
    // //oldmesh->coarse(0.016,1,0,-1);
    // oldmesh->coarse(2.5,0,10,-1);
    // 
    std::cout << "EOF" << std::endl;
}
