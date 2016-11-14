#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <type_traits>


template <typename Visitor, typename Traits, typename Complex, typename K>
struct BFS_Up_Node {};

template <typename Visitor, typename Traits, typename Complex, std::size_t k>
struct BFS_Up_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t,k>>
{
    static constexpr auto level = k;
    using CurrNodeID = typename Complex::template NodeID<level>;
    using NextNodeID = typename Complex::template NodeID<level+1>;
    template <typename T> using Container = typename Traits::template Container<T>;

    using BFS_Up_Node_Next = BFS_Up_Node<Visitor,Traits,Complex,std::integral_constant<std::size_t,level+1>>;

    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, Iterator begin, Iterator end)
    {
        Container<NextNodeID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            if(v.visit(F, *curr))
            {
                for(auto a : F.get_cover(*curr))
                {
                    auto id = F.get_node_up(*curr,a);
                    next.insert(id);
                }
            }
        }

        BFS_Up_Node_Next::apply(v, F, next.begin(), next.end());
    }
};

template <typename Visitor, typename Traits, typename Complex>
struct BFS_Up_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t, Complex::topLevel>>
{
    static constexpr auto level = Complex::topLevel;
    using CurrNodeID = typename Complex::template NodeID<level>;

    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, Iterator begin, Iterator end)
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
    using CurrNodeID = typename Complex::template NodeID<level>;
    using NextNodeID = typename Complex::template NodeID<level-1>;
    template <typename T> using Container = typename Traits::template Container<T>;

    using BFS_Down_Node_Next = BFS_Down_Node<Visitor,Traits,Complex,std::integral_constant<std::size_t,level-1>>;

    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, Iterator begin, Iterator end)
    {
        Container<NextNodeID> next;

        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);

            for(auto a : F.get_name(*curr))
            {
                auto id = F.get_node_down(*curr,a);
                next.insert(id);
            }
        }

        BFS_Down_Node_Next::apply(v, F, next.begin(), next.end());
    }
};

template <typename Visitor, typename Traits, typename Complex>
struct BFS_Down_Node<Visitor, Traits, Complex, std::integral_constant<std::size_t, 1>>
{
    static constexpr auto level = Complex::topLevel;
    using CurrNodeID = typename Complex::template NodeID<level>;

    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, Iterator begin, Iterator end)
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
    using CurrNodeID = typename Complex::template NodeID<level>;
    template <typename T> using Container = typename Traits::template Container<T>;
    using BFS_Edge_Next = BFS_Edge<Visitor,Traits,Complex,std::integral_constant<std::size_t,level+1>>;


    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, Iterator begin, Iterator end)
    {
        Container<NextEdgeID> next;
        std::vector<typename Complex::KeyType> cover;

        for(auto curr = begin; curr != end; ++curr)
        {
            v.visit(F, *curr);

            CurrNodeID n = curr->up();
            F.get_cover(n, std::back_inserter(cover));
            for(auto a : cover)
            {
                NextEdgeID id = F.get_edge_up(n,a);
                next.insert(next.end(), id);
            }
            cover.clear();
        }

        BFS_Edge_Next::apply(v, F, next.begin(), next.end());
    }
};

template <typename Visitor, typename Traits, typename Complex>
struct BFS_Edge<Visitor, Traits, Complex, std::integral_constant<std::size_t, Complex::topLevel>>
{
    static constexpr auto level = Complex::topLevel;
    using CurrEdgeID = typename Complex::template EdgeID<level>;

    template <typename Iterator>
    static void apply(Visitor& v, const Complex& F, Iterator begin, Iterator end)
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
    template <typename T> using Container = std::set<T>;
};

struct BFS_NoRepeat_Edge_Traits
{
    template <typename T> using Container = std::set<T>;
//    template <typename Complex, typename NodeID> auto node_next(Complex F, NodeID s);
};

template <typename Visitor, typename NodeID>
void visit_node_up(Visitor v, const typename NodeID::complex& F, NodeID s)
{
    BFS_Up_Node<Visitor, BFS_NoRepeat_Node_Traits, typename NodeID::complex, std::integral_constant<std::size_t,NodeID::level>>::apply(v,F,&s,&s+1);
}

template <typename Visitor, typename NodeID>
void visit_node_down(Visitor v, const typename NodeID::complex& F, NodeID s)
{
    BFS_Down_Node<Visitor, BFS_NoRepeat_Node_Traits, typename NodeID::complex, std::integral_constant<std::size_t,NodeID::level>>::apply(v,F,&s,&s+1);
}

template <typename Visitor, typename EdgeID>
void edge_up(Visitor v, const typename EdgeID::complex& F, EdgeID s)
{
    BFS_Edge<Visitor, BFS_NoRepeat_Edge_Traits, typename EdgeID::complex, std::integral_constant<std::size_t,EdgeID::level>>::apply(v,F,&s,&s+1);
}
