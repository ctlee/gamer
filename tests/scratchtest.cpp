#include "SurfaceMesh.h"
#include "Vertex.h"
#include <iostream>

int main(int argc, char *argv[])
{
    SurfaceMesh x = SurfaceMesh();
    x.insert<1>({1}, Vertex(1,2,3));
    x.insert<1>({2}, Vertex(3,4,5));
    x.insert<1>({3}, Vertex(5,6,7));
    x.print();
    size_t result = x.remove<1>({1});
    std::cout << result << std::endl;
    x.print();
}
