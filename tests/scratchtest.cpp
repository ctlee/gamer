/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include "SurfaceMeshOld.h"
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


    int trials = 10;
    std::chrono::duration<double> elapsed_seconds;

    for(int trial=0; trial < trials; ++trial){
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        // CODE GOES HERE

        auto mesh = SurfaceMeshOld::readOFF(argv[1]);
        for(int i=0; i < 10; ++i){
            mesh->normalSmooth();
        }

        // CODE ENDS
        end = std::chrono::system_clock::now();
        free(mesh);
        elapsed_seconds += end-start;
    }
    elapsed_seconds /= trials;
    std::cout << "Average time old: " << elapsed_seconds.count() << "s\n";


    
    // New code benchmark 
    for(int trial=0; trial < trials; ++trial){
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        // CODE GOES HERE


        auto result = readOFF(argv[1]);

        if(result.second == false){
            std::cout << "Something bad happened...";
            exit(1);
        }
        auto mesh = result.first;
        compute_orientation(*mesh);

        for(auto i = 0; i < 1; ++i){
            for(auto nid : mesh->get_level_id<1>())
            {
                normalSmooth(*mesh, nid);
            }
        }

        // CODE ENDS
        end = std::chrono::system_clock::now();
        free(mesh);
        elapsed_seconds += end-start;
    }
    elapsed_seconds /= trials;
    std::cout << "Average time new: " << elapsed_seconds.count() << "s\n";

    // std::cout << "Surface Area: " << getArea(*mesh) << std::endl;
    // std::cout << "Volume: " << getVolume(*mesh) << std::endl;

    // // for(auto i = 0; i < 10; ++i){
    // //     for(auto nid : mesh->get_level_id<1>()){
    // //         weightedVertexSmooth<1>(*mesh, nid);
    // //     }
    // // }

    // // writeOFF("../data/vsmooth.off", *mesh);
    // for(auto i = 0; i < 1; ++i){
    //     for(auto nid : mesh->get_level_id<1>())
    //     {
    //         normalSmooth(*mesh, nid);
    //     }
    // }

    // std::cout << "Surface Area: " << getArea(*mesh) << std::endl;
    // std::cout << "Volume: " << getVolume(*mesh) << std::endl;

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

    // writeOFF("../data/smoothed.off", *mesh);
    std::cout << "EOF" << std::endl;
}
