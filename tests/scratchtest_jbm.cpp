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
std::ostream& operator<<(std::ostream& out, const std::set<T>& A)
{
    out << "[";
    for(auto a : A)
    {
        std::cout << a << " ";
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
auto make_print_visitor(const Complex& F)
{
    return PrintVisitor<Complex>();
}


namespace jbm
{
    template <class T> using set = std::set<T>;
    template <typename Complex>
    struct SimplexSet
    {
        template <std::size_t j>
        using Simplex = typename Complex::template NodeID<j>;
        using LevelIndex = typename std::make_index_sequence<Complex::numLevels>;
        using NodeIDLevel = typename util::int_type_map<std::size_t, std::tuple, LevelIndex, Simplex>::type;
        using type = typename util::type_map<NodeIDLevel, jbm::set>::type;
    };
}

template <typename Complex>
struct TestAVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;

    TestAVisitor(SimplexSet* p) : pLevels(p) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        if(std::get<level>(*pLevels).find(s) == std::get<level>(*pLevels).end())
        {
            std::get<level>(*pLevels).insert(s);
//            std::cout << F.get_name(s) << std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }

private:
    SimplexSet* pLevels;
};



template <typename Complex>
struct TestBVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;

    TestBVisitor(SimplexSet* p) : pLevels(p) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        return true;
    }

    bool visit(const Complex& F, typename Complex::template NodeID<1> s)
    {
        visit_node_up(TestAVisitor<Complex>(pLevels), F, s);
        return false;
    }

private:
    SimplexSet* pLevels;
};


template <typename Complex>
struct GrabVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;

    GrabVisitor(SimplexSet* p) : pLevels(p) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        if(std::get<level>(*pLevels).find(s) != std::get<level>(*pLevels).end())
        {
            std::get<level>(*pLevels).erase(s);
//            std::cout << F.get_name(s) << std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }

private:
    SimplexSet* pLevels;
};



template <typename Complex, std::size_t level>
struct InnerVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;
    using Simplex = typename Complex::template NodeID<level>;
    using KeyType = typename Complex::KeyType;

    InnerVisitor(SimplexSet* p, Simplex s) : pLevels(p), simplex(s), new_point(-1) {}

    template <std::size_t k>
    bool visit(const Complex& F, typename Complex::template NodeID<k> s)
    {
        if(std::get<k>(*pLevels).find(s) != std::get<k>(*pLevels).end())
        {
            std::set<KeyType> new_name;
            for(auto a : F.get_name(s))
            {
                new_name.insert(a);
            }
            for(auto a : F.get_name(simplex))
            {
                new_name.erase(a);
            }
            new_name.insert(new_point);
            // This print out gives me a list of all the changes which need to occur.
            // And the subsequent visitor gets all the simplices which will merge into new_name.
            // However, these changes cannot be applied in the order they are generated here.
            // These changes must be stored and executed later.
            std::cout << "Inner: " << F.get_name(s) << " --> " << new_name << std::endl;
            visit_node_down(GrabVisitor<Complex>(pLevels), F, s);
        }
        return true;
    }

private:
    SimplexSet* pLevels;
    Simplex     simplex;
    KeyType     new_point;
};


template <typename Complex>
struct MainVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;

    MainVisitor(SimplexSet* p) : pLevels(p) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
//        std::cout << "MAIN: " << F.get_name(s) << std::endl;
        visit_node_up(InnerVisitor<Complex,level>(pLevels,s), F, s);
        return true;
    }

private:
    SimplexSet* pLevels;
};



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
    
    for(auto s : mesh->get_level_id<2>())
    {
        std::cout << mesh->get_name(s) << std::endl;
        typename jbm::SimplexSet<SurfaceMesh>::type levels;
        visit_node_down(TestBVisitor<SurfaceMesh>(&levels), *mesh, s);
        visit_node_down(MainVisitor<SurfaceMesh>(&levels), *mesh, s);
//        visit_node_down(make_testB_visitor(mesh->get_node_up<1>({1})), *mesh, s);
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
