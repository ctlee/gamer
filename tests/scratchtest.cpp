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
    auto nid = *(++++mesh->get_level_id<1>().begin());
    std::cout << *nid << std::endl << std::endl;
    
    std::set<SurfaceMesh::NodeID<1>> nodes;

    neighbors_up(*mesh, nid, std::inserter(nodes, nodes.begin()));
    auto next = nodes;
    for(auto node : next){
        neighbors_up(*mesh, node, std::inserter(nodes, nodes.begin()));
    }
    for(auto node : nodes){
        std::cout << *node << std::endl;
    }
    
    std::cout << std::endl << std::endl;
    auto nodes2 = neighbors_up(*mesh, nid, 3);

    for(auto node : nodes2){
        std::cout << *node << std::endl;
    }
    // writeOFF("../data/test.off", *mesh);
    std::cout << "EOF" << std::endl;
}
