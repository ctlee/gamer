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
    

    init_orientation(*mesh);
    clear_orientation(*mesh);
    auto orient = compute_orientation(*mesh);

    std::cout << "Connected Components: " << std::get<0>(orient) << std::endl;
    std::cout << "Orientable: " << std::get<1>(orient) << std::endl;
    std::cout << "Psuedo-manifold: " << std::get<2>(orient) << std::endl;

    //print(*mesh);
    mesh->genGraph("test.dot");

    std::cout << "Generating Histogram..." << std::endl;
    generateHistogram(*mesh);

    auto edges = selectFlipEdgesByAngle(*mesh, true);

    for (auto edge : edges){
        edgeFlip(*mesh, edge);
    }

    // Fishy code to draw the tangent planes...
    // auto range = mesh->get_level_id<1>();
    // auto nd = --range.end();
    // int name = mesh->size<1>()-1;
    // for (auto node : range){
    //     auto T = getTangent(*mesh, node);
    //     T = T/std::sqrt(T|T);
    //     std::cout << *node << std::endl;
    //     std::cout << T << std::endl;
    //     auto v = *node;
    //     for(int i = 0; i < 3; i++){
    //         auto q = v + Vertex(T.get(0,i), T.get(1,i), T.get(2,i));
    //         std::cout << q << std::endl;
    //         mesh->insert<1>({name+i+1}, q);
    //     }
    //     mesh->insert<3>({name+1, name+2, name+3});
    //     name += 3;
    //     if (node == *nd) break;
    // }
    writeOFF("test.off", *mesh);
    std::cout << "EOF" << std::endl;
}
