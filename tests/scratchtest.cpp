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
    SurfaceMesh x = SurfaceMesh();
    x.insert<1>({1}, Vertex(1,2,3));
    x.insert<1>({2}, Vertex(3,4,5));
    x.insert<1>({3}, Vertex(5,6,7));
    x.insert<1>({6}, Vertex(7,8,9));
    x.insert<1>({5}, Vertex(9,10,11));
    x.insert<1>({4}, Vertex(11,12,13));
    
    x.insert<3>({1,2,3});
    x.insert<1>({7}, Vertex(-1,-1,-1));
    x.insert<3>({2,3,4});
    print_vertices(x);

    Vertex v = x.get<1>({7});
    std::cout << "Node<7>=" << v << std::endl;
    print_nodes<0>(x);
    print_nodes<1>(x);
    print_nodes<2>(x);
    print_nodes<3>(x);
    std::cout<<std::endl;

    writeOFF("test.off", x);
    
    std::cout << "Try to get neighbors" << std::endl;
    
    for(auto nodeID : x.get_level_id<1>()){
        std::cout << *nodeID << std::endl;
        auto w = x.get_cover<1>(nodeID);
        for(int i=0; i < w.size(); i++)
            std::cout << w[i] << std::endl;
        std::vector<SurfaceMesh::NodeID<1>> tmp;
        neighbors(x, nodeID, std::back_inserter(tmp));
    }

    x.remove<1>({1});

    v = x.get<1>({7});
    std::cout << "Node<7>=" << v << std::endl;

    print_nodes<0>(x);
    print_nodes<1>(x);
    print_nodes<2>(x);
    print_nodes<3>(x);

    std::cout << "EOF" << std::endl;
}
