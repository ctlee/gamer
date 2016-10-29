/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include "SurfaceMesh.h"
#include "Vertex.h"
#include <vector>
#include <iostream>
#include <string>


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

    init_orientation(*mesh);
    clear_orientation(*mesh);
    compute_orientation(*mesh);
    
    for(auto v : mesh->get_level_id<1>())
    {
        auto T = getTangent(*mesh, v);
        std::cout << T << std::endl;
        /*
        auto edges = mesh->up(v);
        auto faces = mesh->up(edges);

        std::cout << v << std::endl;
        for(auto f : faces)
        {
            std::cout << "    " << f << std::endl;
        }
        */
    }
}
