/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <libraries/Eigen/Dense>
#include <libraries/Eigen/Eigenvalues>
#include "SurfaceMesh.h"
#include "SimplicialComplexVisitors.h"
#include "Vertex.h"



template <typename T, std::size_t k>
std::ostream& operator<<(std::ostream& out, const std::array<T,k>& A)
{
    out << "[";
    for(int i = 0; i + 1 < k; ++i)
    {
        out << A[i] << " ";
    }
    if(k > 0)
    {
        out << A[k-1];
    }
    out << "]";
    return out;
}

template <typename Complex>
struct PrintVisitor
{
    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        std::cout << F.get_name(s) << std::endl;
        return true;
    }
};

template <typename Complex>
auto make_print_visitor(const Complex& F)
{
    return PrintVisitor<Complex>();
}

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
    compute_orientation(*mesh);

    for(auto nid: mesh->get_level_id<1>())
    {
        weightedVertexSmooth(*mesh, nid);
    }
    // for(auto nid : mesh->get_level_id<3>()){
    //     //auto v = LocalStructureTensorVisitor();
    //     // auto v = make_print_visitor(*mesh);
    //     // visit_neighbors_down<1>(v, *mesh, nid);
    //     visit_neighbors_down<1>(make_print_visitor(*mesh), *mesh, nid);
    //     //std::cout << *nid << " " << v.lst << std::endl;
    // }

    // for(auto nid : mesh->get_level_id<1>()){
    //     auto v = LocalStructureTensorVisitor();
    //     visit_neighbors_up<1>(v, *mesh, nid);
    //     std::cout << *nid << " " << v.lst << std::endl;
    //     getEigenvalues(v.lst);
    // }
    
    writeOFF("../data/test.off", *mesh);
    std::cout << "EOF" << std::endl;
}
