#include "SurfaceMesh.h"
#include "Vertex.h"
#include <iostream>

int main(int argc, char *argv[])
{
    SurfaceMesh x = SurfaceMesh();
    x.insert<1>({1}, Vertex(1,2,3));
    x.insert<1>({2}, Vertex(3,4,5));
    x.insert<1>({3}, Vertex(5,6,7));
    x.insert<1>({4}, Vertex(7,8,9));
    x.insert<1>({5}, Vertex(9,10,11));
    x.insert<1>({6}, Vertex(11,12,13));
    
    x.insert<3>({1,2,3});
    x.insert<1>({7}, Vertex(-1,-1,-1));
    x.insert<3>({2,3,4});
    x.print();

    Vertex v = x.get<1>({7});
    std::cout << "Node<7>=" << v << std::endl;
    x.print_id<0>();
    x.print_id<1>();
    x.print_id<2>();
    x.print_id<3>();
    std::cout<<std::endl;
    
    x.remove<1>({1});

    v = x.get<1>({7});
    std::cout << "Node<7>=" << v << std::endl;
    x.print_id<0>();
    x.print_id<1>();
    x.print_id<2>();
    x.print_id<3>();
    std::cout << "EOF" << std::endl;
}
