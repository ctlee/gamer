/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include "SurfaceMesh.h"
#include "Vertex.h"
#include "SimplicialComplexVisitors.h"
#include <vector>
#include <iostream>
#include <string>
#include <type_traits>
#include <list>


/*

template <std::size_t dimension>
auto getTangentH(SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::NodeID<SurfaceMesh::topLevel> curr)
{
    return 1.0;
}

template <std::size_t level, std::size_t dimension>
auto getTangentH(SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::NodeID<level> curr)
{
    tensor<double, dimension, SurfaceMesh::topLevel - level> rval;

    for(auto alpha : mesh.get_cover(curr))
    {
        auto edge = *mesh.get_edge_up(curr, alpha);
        const auto& v = (*mesh.get_node_up({alpha})).position;
        auto next = mesh.get_node_up(curr,alpha);
        rval += edge.orientation * (v-origin) * getTangentH(mesh, origin, next);
    }

    return rval;
}

auto getTangent_jbm(SurfaceMesh& mesh, SurfaceMesh::NodeID<1> vertexID)
{
    return getTangentH(mesh, (*vertexID).position, vertexID);
}

template <class Complex>
void SimplexCollapse(const Complex& F)
{

}
*/

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

template <typename T>
void print_name(const T& name)
{
    for(auto x : name)
    {
        std::cout << x << " ";
    }
}

template <typename NodeID>
void grab(const typename NodeID::complex& F, NodeID s)
{
    std::cout << F.get_name(s) << std::endl;
}

template <typename Complex>
struct PrintEdgeVisitor
{
    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template EdgeID<level> s)
    {
        auto down = s.down();
        std::cout << F.get_name(down);
        std::cout << " -> ";
        std::cout << F.get_name(s.up()) << std::endl;
        return true;
    }
};
template <typename Complex>
PrintEdgeVisitor<Complex> make_print_edge_visitor(const Complex& F)
{
    return PrintEdgeVisitor<Complex>();
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
PrintVisitor<Complex> make_print_visitor(const Complex& F)
{
    return PrintVisitor<Complex>();
}


template <typename Complex>
struct TestAVisitor
{
    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        std::cout << F.get_name(s) << std::endl;
        return true;
    }
};

template <typename Complex>
PrintVisitor<Complex> make_testA_visitor(const Complex& F)
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

    if(result.second == false)
    {
        std::cout << "Something bad happened...";
        exit(1);
    }
    auto mesh = result.first;

    compute_orientation(*mesh);
    
    {
        auto s = *(++++mesh->get_level_id<3>().begin());

        visit_node_down(make_print_visitor(*mesh), *mesh, s);
//        edge_up(make_print_edge_visitor(*mesh), *mesh, mesh->get_edge_up(mesh->get_node_up(),1));
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
