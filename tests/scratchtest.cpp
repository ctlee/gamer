/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include "ASCFunctions.h"
#include "stringutil.h"
#include "SurfaceMeshOld.h"
#include "SurfaceMesh.h"
#include "SimplicialComplexVisitors.h"
#include "Vertex.h"

int  main(int argc, char *argv[])
{

    auto mesh = AbstractSimplicialComplex<int,int,int,int,int,int>();

    mesh.insert<3>({0,1,2});

    auto v = mesh.get_simplex_up({0,1,2});
    if(v != nullptr){
        for(auto c : mesh.get_name(v)){
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    else{
        std::cout << "Simplex v doesn't exist" << std::endl;
    }
    // if(argc != 2)
    // {
    //     std::cerr << "Wrong arguments passed" << std::endl;
    //     return -1;
    // }
    // auto mesh = readOBJ(argv[1]);
    // if(mesh == nullptr){
    //     std::cout << "Something bad happened...";
    //     exit(1);
    // }
    // std::cout << "Done reading" << std::endl;
    // compute_orientation(*mesh);
    // auto volume = getVolume(*mesh);
    
    // writeOFF("test.off", *mesh);

    // std::cout << "Volume: " << volume << std::endl;
    // //mesh->genGraph("test.dot");
    // writeDOT("test.dot", *mesh);

    // int trials = 10;
    // std::chrono::duration<double> elapsed_seconds;

    // for(int trial=0; trial < trials; ++trial){
    //     std::chrono::time_point<std::chrono::system_clock> start, end;
    //     start = std::chrono::system_clock::now();
    //     // CODE GOES HERE

    //     auto mesh = SurfaceMeshOld::readOFF(argv[1]);
    //     for(int i=0; i < 10; ++i){
    //         mesh->normalSmooth();
    //     }

    //     // CODE ENDS
    //     end = std::chrono::system_clock::now();
    //     free(mesh);
    //     elapsed_seconds += end-start;
    // }
    // elapsed_seconds /= trials;
    // std::cout << "Average time old: " << elapsed_seconds.count() << "s\n";
    
    // // New code benchmark 
    // for(int trial=0; trial < trials; ++trial){
    //     std::chrono::time_point<std::chrono::system_clock> start, end;
    //     start = std::chrono::system_clock::now();
    //     // CODE GOES HERE

    //     auto mesh = readOFF(argv[1]);
    //     if(mesh == nullptr){
    //         std::cout << "Something bad happened...";
    //         exit(1);
    //     }

    //     compute_orientation(*mesh);

    //     for(auto i = 0; i < 1; ++i){
    //         for(auto nid : mesh->get_level_id<1>())
    //         {
    //             normalSmooth(*mesh, nid);
    //         }
    //     }

    //     // CODE ENDS
    //     end = std::chrono::system_clock::now();
    //     elapsed_seconds += end-start;
    // }
    // elapsed_seconds /= trials;
    // std::cout << "Average time new: " << elapsed_seconds.count() << "s\n";


    std::cout << "EOF" << std::endl;
}
