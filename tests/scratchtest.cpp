/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <vector>
#include <iostream>
#include <string>
#include <cmath>

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

template <typename Visitor, typename Complex, typename K, typename R>
struct Neighbors_Up_Node {};


template <typename Visitor, typename Complex, std::size_t k, std::size_t ring>
struct Neighbors_Up_Node<Visitor, Complex, std::integral_constant<std::size_t, k>, std::integral_constant<std::size_t, ring>>
{
    static constexpr auto level = k;
    using NodeID = typename Complex::template NodeID<level>;

    using Neighbors_Up_Node_Next = Neighbors_Up_Node<Visitor,Complex,
            std::integral_constant<std::size_t,level>,std::integral_constant<std::size_t,ring-1>>;

    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, std::set<NodeID>& nodes, Iterator begin, Iterator end)
    {
        std::set<NodeID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            if(v.visit(F, *curr))
            {
                for(auto a : F.get_cover(*curr))
                {
                    auto id = F.get_node_up(*curr,a);
                    for(auto b : F.get_name(id)) 
                    {
                        auto nbor = F.get_node_down(id,b);
                        if(nodes.insert(nbor).second)
                        {
                            next.insert(nbor);
                        }
                    }
                }
            }
        }
        Neighbors_Up_Node_Next::apply(v, F, nodes, next.begin(), next.end());
    }
};

template <typename Visitor, typename Complex, std::size_t k>
struct Neighbors_Up_Node<Visitor, Complex, std::integral_constant<std::size_t, k>, std::integral_constant<std::size_t, 0>>
{
    static constexpr auto level = k;
    using NodeID = typename Complex::template NodeID<level>;
    
    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, std::set<NodeID>& nodes, Iterator begin, Iterator end)
    {
        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);
        }
    }
};

template <std::size_t rings, typename Visitor, typename NodeID>
void visit_neighbors_up(Visitor v, const typename NodeID::complex& F, NodeID s)
{
    std::set<NodeID> nodes{s};
    Neighbors_Up_Node<Visitor, typename NodeID::complex, 
            std::integral_constant<std::size_t,NodeID::level>, 
            std::integral_constant<std::size_t, rings>>::apply(v,F,nodes,&s,&s+1);
}


template <typename Complex>
struct PrintVisitor
{
    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        std::cout << *s << std::endl;
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


    auto nid = *(++++mesh->get_level_id<1>().begin());
    std::cout << *nid << std::endl << std::endl;

    visit_neighbors_up<2>(make_print_visitor(*mesh), *mesh, nid);
    // std::set<SurfaceMesh::NodeID<1>> nodes;

    // neighbors_up(*mesh, nid, std::inserter(nodes, nodes.begin()));
    // auto next = nodes;
    // for(auto node : next){
    //     neighbors_up(*mesh, node, std::inserter(nodes, nodes.begin()));
    // }
    // for(auto node : nodes){
    //     std::cout << *node << std::endl;
    // }
    
    // std::cout << std::endl << std::endl;
    // auto nodes2 = neighbors_up(*mesh, nid, 3);

    // for(auto node : nodes2){
    //     std::cout << *node << std::endl;
    // }
    // writeOFF("../data/test.off", *mesh);
    std::cout << "EOF" << std::endl;
}
