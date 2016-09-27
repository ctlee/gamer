#include <iostream>
#include <map>
#include <cmath>
#include "Vertex.h"


int main(int argc, char *argv[])
{
    Vertex v0 = Vertex(0,0,0);
    Vertex v1 = Vertex(0,0,0);
    Vertex v2 = Vertex(1,1,1);
    Vertex v3 = Vertex(2,3,4);

    assert(v0 == v1);
    std::cout << "v0 == v1" << std::endl;
    assert(v0 != v2);
    std::cout << "v0 != v2" << std:: endl;
    assert(v2 + v3 == Vertex(3,4,5));
    std::cout << "v2+v3" << std::endl;

    std::cout << "Test Complete" << std::endl;
	return -1;
}
