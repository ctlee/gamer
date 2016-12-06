/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>
#include "SurfaceMesh.h"
#include "SimplicialComplexVisitors.h"
#include "Vertex.h"

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    std::cout << "Begin reading Mesh..." << std::endl;
    auto result = readOFF(argv[1]);

    if(result.second == false){
        std::cout << "Something bad happened...";
        exit(1);
    }
    auto mesh = result.first;
   

    compute_orientation(*mesh);
    std::cout << "Surface Area: " << getArea(*mesh) << std::endl;
    std::cout << "Volume: " << getVolume(*mesh) << std::endl;

    // for(auto i = 0; i < 10; ++i){
    //     for(auto nid : mesh->get_level_id<1>()){
    //         weightedVertexSmooth<1>(*mesh, nid);
    //     }
    // }

    // writeOFF("../data/vsmooth.off", *mesh);
    for(auto i = 0; i < 1; ++i){
        for(auto nid : mesh->get_level_id<1>())
        {
            normalSmooth(*mesh, nid);
        }
    }

    std::cout << "Surface Area: " << getArea(*mesh) << std::endl;
    std::cout << "Volume: " << getVolume(*mesh) << std::endl;

    // compute_orientation(*mesh);

    // generateHistogram(*mesh);

    // for(auto i = 0; i < 10; ++i){
    //     for(auto nid : mesh->get_level_id<1>())
    //     {
    //         weightedVertexSmooth<3>(*mesh, nid);
    //     }
    // }
    
    // generateHistogram(*mesh);
    // compute_orientation(*mesh);

    writeOFF("../data/smoothed.off", *mesh);
    std::cout << "EOF" << std::endl;
}
