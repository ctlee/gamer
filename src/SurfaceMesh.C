#include "SurfaceMesh.h"

void print_vertices(const SurfaceMesh& mesh){
    for(auto x: mesh.get_level<1>()) {
        std::cout << x << ", ";
    }
    std::cout << std::endl;
}