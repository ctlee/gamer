/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include "SurfaceMesh.h"
#include "Vertex.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

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
    mesh->remove({0});

    mesh->insert({4}, Vertex(0,0,2));
    mesh->insert({1,4,3});
    mesh->insert({2,4,1});
    mesh->insert({3,4,2});


    mesh->renumber();
    writeOFF("../data/test.off", *mesh);
    std::cout << "EOF" << std::endl;
}
