#pragma once

#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include "SimplicialComplex.h"

template <typename Visitor, typename Traits, typename Complex, typename K>
struct BFS_Up_Node {};

template <typename Visitor, typename Traits, typename Complex, std::size_t k>
struct BFS_Up_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t,k>>
{
    static constexpr auto level = k;
    using CurrSimplexID = typename Complex::template SimplexID<level>;
    using NextSimplexID = typename Complex::template SimplexID<level+1>;
    template <typename T> using Container = typename Traits::template Container<T>;

    using BFS_Up_Node_Next = BFS_Up_Node<Visitor,Traits,Complex,std::integral_constant<std::size_t,level+1>>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, Iterator begin, Iterator end)
    {
        Container<NextSimplexID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            if(v.visit(F, *curr))
            {
                F.get_cover(*curr, [&](typename Complex::KeyType a)
                {
                    auto id = F.get_simplex_up(*curr,a);
                    next.insert(id);
                });
            }
        }

        BFS_Up_Node_Next::apply(std::forward<Visitor>(v), F, next.begin(), next.end());
    }
};

template <typename Visitor, typename Traits, typename Complex>
struct BFS_Up_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t, Complex::topLevel>>
{
    static constexpr auto level = Complex::topLevel;
    using CurrSimplexID = typename Complex::template SimplexID<level>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, Iterator begin, Iterator end)
    {
        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);
        }
    }
};



template <typename Visitor, typename Traits, typename Complex, typename K>
struct BFS_Down_Node {};

template <typename Visitor, typename Traits, typename Complex, std::size_t k>
struct BFS_Down_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t,k>>
{
    static constexpr auto level = k;
    using CurrSimplexID = typename Complex::template SimplexID<level>;
    using NextSimplexID = typename Complex::template SimplexID<level-1>;
    template <typename T> using Container = typename Traits::template Container<T>;

    using BFS_Down_Node_Next = BFS_Down_Node<Visitor,Traits,Complex,std::integral_constant<std::size_t,level-1>>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, Iterator begin, Iterator end)
    {
        Container<NextSimplexID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);

            F.get_name(*curr, [&](typename Complex::KeyType a)
            {
                auto id = F.get_simplex_down(*curr,a);
                next.insert(id);
            });
        }

        BFS_Down_Node_Next::apply(std::forward<Visitor>(v), F, next.begin(), next.end());
    }
};

template <typename Visitor, typename Traits, typename Complex>
struct BFS_Down_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t, 1>>
{
    static constexpr auto level = Complex::topLevel;
    using CurrSimplexID = typename Complex::template SimplexID<level>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, Iterator begin, Iterator end)
    {
        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);
        }
    }
};




template <typename Visitor, typename Traits, typename Complex, typename K>
struct BFS_Edge {};

template <typename Visitor, typename Traits, typename Complex, std::size_t k>
struct BFS_Edge<Visitor, Traits, Complex, std::integral_constant<std::size_t,k>>
{
    static constexpr auto level = k;
    using CurrEdgeID = typename Complex::template EdgeID<level>;
    using NextEdgeID = typename Complex::template EdgeID<level+1>;
    using CurrSimplexID = typename Complex::template SimplexID<level>;
    template <typename T> using Container = typename Traits::template Container<T>;
    using BFS_Edge_Next = BFS_Edge<Visitor,Traits,Complex,std::integral_constant<std::size_t,level+1>>;


    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, Iterator begin, Iterator end)
    {
        Container<NextEdgeID> next;
        std::vector<typename Complex::KeyType> cover;

        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);

            CurrSimplexID n = curr->up();
            F.get_cover(n, std::back_inserter(cover));
            for(auto a : cover)
            {
                NextEdgeID id = F.get_edge_up(n,a);
                next.insert(next.end(), id);
            }
            cover.clear();
        }

        BFS_Edge_Next::apply(std::forward<Visitor>(v), F, next.begin(), next.end());
    }
};

template <typename Visitor, typename Traits, typename Complex>
struct BFS_Edge<Visitor, Traits, Complex, std::integral_constant<std::size_t, Complex::topLevel>>
{
    static constexpr auto level = Complex::topLevel;
    using CurrEdgeID = typename Complex::template EdgeID<level>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, Iterator begin, Iterator end)
    {
        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);
        }
    }
};

template <typename Visitor, typename Complex, std::size_t k, std::size_t ring>
struct Neighbors_Up_Node
{
    static constexpr auto level = k;
    using SimplexID = typename Complex::template SimplexID<level>;

    using Neighbors_Up_Node_Next = Neighbors_Up_Node<Visitor,Complex,level,ring-1>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, NodeSet<SimplexID>& nodes, Iterator begin, Iterator end)
    {
        NodeSet<SimplexID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            if(v.visit(F, *curr))
            {
                for(auto a : F.get_cover(*curr))
                {
                    auto id = F.get_simplex_up(*curr,a);
                    for(auto b : F.get_name(id)) 
                    {
                        auto nbor = F.get_simplex_down(id,b);
                        if(nodes.insert(nbor).second)
                        {
                            next.insert(nbor);
                        }
                    }
                }
            }
        }

        Neighbors_Up_Node_Next::apply(std::forward<Visitor>(v), F, nodes, next.begin(), next.end());
    }
};

template <typename Visitor, typename Complex, std::size_t k>
struct Neighbors_Up_Node<Visitor, Complex, k, 0>
{
    static constexpr auto level = k;
    using SimplexID = typename Complex::template SimplexID<level>;
    
    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, NodeSet<SimplexID>& nodes, Iterator begin, Iterator end)
    {
        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);
        }
    }
};



template <typename Visitor, typename Complex, std::size_t k, std::size_t ring>
struct Neighbors_Down_Node
{
    static constexpr auto level = k;
    using SimplexID = typename Complex::template SimplexID<level>;

    using Neighbors_Down_Node_Next = Neighbors_Down_Node<Visitor,Complex,
            level,ring-1>;

    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, NodeSet<SimplexID>& nodes, Iterator begin, Iterator end)
    {
        NodeSet<SimplexID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            if(v.visit(F, *curr))
            {
                for(auto a : F.get_name(*curr))
                {
                    auto id = F.get_simplex_down(*curr,a);
                    for(auto b : F.get_cover(id)) 
                    {
                        auto nbor = F.get_simplex_up(id,b);
                        if(nodes.insert(nbor).second)
                        {
                            next.insert(nbor);
                        }
                    }
                }
            }
        }

        Neighbors_Down_Node_Next::apply(std::forward<Visitor>(v), F, nodes, next.begin(), next.end());
    }
};

template <typename Visitor, typename Complex, std::size_t k>
struct Neighbors_Down_Node<Visitor, Complex, k, 0>
{
    static constexpr auto level = k;
    using SimplexID = typename Complex::template SimplexID<level>;
    
    template <typename Iterator>
    static void apply(Visitor&& v, const Complex& F, NodeSet<SimplexID>& nodes, Iterator begin, Iterator end)
    {
        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);
        }
    }
};



template <typename T> using AllowRepeat = std::vector<T>;

struct BFS_NoRepeat_Node_Traits
{
    template <typename T> using Container = NodeSet<T>;
};

struct BFS_NoRepeat_Edge_Traits
{
    template <typename T> using Container = NodeSet<T>;
//    template <typename Complex, typename SimplexID> auto node_next(Complex F, SimplexID s);
};

template <typename Visitor, typename SimplexID>
void visit_node_up(Visitor&& v, const typename SimplexID::complex& F, SimplexID s)
{
    BFS_Up_Node<Visitor, BFS_NoRepeat_Node_Traits, typename SimplexID::complex, std::integral_constant<std::size_t,SimplexID::level>>::apply(std::forward<Visitor>(v),F,&s,&s+1);
}

template <typename Visitor, typename SimplexID>
void visit_node_down(Visitor&& v, const typename SimplexID::complex& F, SimplexID s)
{
    BFS_Down_Node<Visitor, BFS_NoRepeat_Node_Traits, typename SimplexID::complex, std::integral_constant<std::size_t,SimplexID::level>>::apply(std::forward<Visitor>(v),F,&s,&s+1);
}

template <typename Visitor, typename EdgeID>
void edge_up(Visitor&& v, const typename EdgeID::complex& F, EdgeID s)
{
    BFS_Edge<Visitor, BFS_NoRepeat_Edge_Traits, typename EdgeID::complex, std::integral_constant<std::size_t,EdgeID::level>>::apply(std::forward<Visitor>(v),F,&s,&s+1);
}

template <std::size_t rings, typename Visitor, typename SimplexID>
void visit_neighbors_up(Visitor&& v, const typename SimplexID::complex& F, SimplexID s)
{
    NodeSet<SimplexID> nodes{s};
    Neighbors_Up_Node<Visitor,typename SimplexID::complex,SimplexID::level,rings>::apply(std::forward<Visitor>(v),F,nodes,&s,&s+1);
}

template <std::size_t rings, typename Visitor, typename SimplexID>
void visit_neighbors_down(Visitor&& v, const typename SimplexID::complex& F, SimplexID s)
{
    NodeSet<SimplexID> nodes{s};
    Neighbors_Down_Node<Visitor,typename SimplexID::complex,SimplexID::level,rings>::apply(std::forward<Visitor>(v),F,nodes,&s,&s+1);
}