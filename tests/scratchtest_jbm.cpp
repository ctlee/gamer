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
    template <class T> using vector = std::vector<T>;

    template <typename Complex>
    struct SimplexSet
    {
        template <std::size_t j>
        using Simplex = typename Complex::template NodeID<j>;
        using LevelIndex = typename std::make_index_sequence<Complex::numLevels>;
        using NodeIDLevel = typename util::int_type_map<std::size_t, std::tuple, LevelIndex, Simplex>::type;
        using type = typename util::type_map<NodeIDLevel, jbm::set>::type;
    };

    template <typename Complex>
    struct SimplexDataSet
    {
        using KeyType = typename Complex::KeyType;

        template <std::size_t k, typename T>
        struct DataType
        {
            using type = std::pair<std::array<KeyType,k>, T>;
        };

        template <std::size_t k>
        struct DataType<k, void>
        {
            using type = std::array<KeyType,k>;
        };

        template <std::size_t j>
        using DataSet = typename DataType<j, typename Complex::template NodeData<j>>::type;
        using LevelIndex = typename std::make_index_sequence<Complex::numLevels>;
        using NodeIDLevel = typename util::int_type_map<std::size_t, std::tuple, LevelIndex, DataSet>::type;
        using type = typename util::type_map<NodeIDLevel, jbm::vector>::type;
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

    GrabVisitor(SimplexSet* p, SimplexSet* grab) : pLevels(p), pGrab(grab) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
        if(std::get<level>(*pLevels).find(s) != std::get<level>(*pLevels).end())
        {
            std::get<level>(*pLevels).erase(s);
            std::get<level>(*pGrab).insert(s);
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
    SimplexSet* pGrab;
};




template <typename Complex, std::size_t BaseLevel, template <typename> class Callback>
struct InnerVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;
    using Simplex = typename Complex::template NodeID<BaseLevel>;
    using KeyType = typename Complex::KeyType;

    InnerVisitor(SimplexSet* p, Simplex s, KeyType np, Callback<Complex>* c) : pLevels(p), simplex(s), new_point(np), callback(c) {}

    template <std::size_t OldLevel>
    bool visit(const Complex& F, typename Complex::template NodeID<OldLevel> s)
    {
        constexpr std::size_t NewLevel = OldLevel - BaseLevel + 1;

        if(std::get<OldLevel>(*pLevels).find(s) != std::get<OldLevel>(*pLevels).end())
        {
            auto old_name = F.get_name(s);
            auto base_name = F.get_name(simplex);
            std::array<KeyType,NewLevel> new_name;

            std::size_t i = 0; // new_name
            std::size_t j = 0; // old_name
            std::size_t k = 0; // base_name

            new_name[i++] = new_point;

            while(i < NewLevel)
            {
                if(base_name[k] == old_name[j])
                {
                    ++j; ++k;
                }
                else
                {
                    new_name[i++] = old_name[j++];
                }
            }

            SimplexSet grab;
            visit_node_down(GrabVisitor<Complex>(pLevels,&grab), F, s);
            (*callback)(F, old_name, new_name, grab);
        }
        return true;
    }

private:
    SimplexSet* pLevels;
    Simplex     simplex;
    KeyType     new_point;
    Callback<Complex>* callback;
};


template <typename Complex, template <typename> class Callback>
struct MainVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;
    using KeyType = typename Complex::KeyType;

    MainVisitor(SimplexSet* p, Callback<Complex>* c, KeyType np) : pLevels(p), callback(c), new_point(np) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template NodeID<level> s)
    {
//        std::cout << "MAIN: " << F.get_name(s) << std::endl;
        visit_node_up(InnerVisitor<Complex,level,Callback>(pLevels,s,new_point,callback), F, s);
        return true;
    }

private:
    SimplexSet* pLevels;
    Callback<Complex>* callback;
    KeyType     new_point;
};


template <typename Complex>
struct Callback
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;
    template <std::size_t level> using Type = typename Complex::template NodeData<level>;
    using KeyType = typename Complex::KeyType;

    template <std::size_t OldLevel, std::size_t NewLevel>
    void operator()(const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,NewLevel>& new_name,
                    const SimplexSet& merged)
    {
        std::cout << "Inner: " << old_name << " --> " << new_name << std::endl;
    }

    template <std::size_t OldLevel>
    void operator()(const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,1>& new_name,
                    const SimplexSet& merged)
    {
        Vertex center;
        std::size_t cnt = 0;
        for(auto v : std::get<1>(merged))
        {
            center = center + (*v);
            ++cnt;
        }
        center = center / (double)(cnt);
        std::cout << "1nner: " << old_name << " --> " << new_name << " : " << center << std::endl;

        std::get<1>(data).push_back(std::make_pair(new_name, center));
    }

    template <std::size_t OldLevel>
    void operator()(const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,3>& new_name,
                    const SimplexSet& merged)
    {
        std::cout << "3nner: " << old_name << " --> " << new_name << " : " << std::get<3>(merged).size() << std::endl;
        std::get<3>(data).push_back(std::make_pair(new_name,Face()));
    }

    typename jbm::SimplexDataSet<Complex>::type data;
};


template <typename Complex, typename NodeDataType, typename T>
struct PerformInsertion {};

template <typename Complex, typename NodeDataType, std::size_t level>
struct PerformInsertion<Complex, NodeDataType, std::integral_constant<std::size_t,level>>
{
    static void apply(Complex& F, typename jbm::SimplexDataSet<Complex>::type& data)
    {
        for(auto curr : std::get<level>(data))
        {
//            std::cout << curr.first << " " << curr.second << std::endl;
            F.template insert<level>(curr.first, curr.second);
        }
        PerformInsertion<Complex,typename Complex::template NodeData<level+1>,std::integral_constant<std::size_t,level+1>>::apply(F,data);
    }
};

template <typename Complex, std::size_t level>
struct PerformInsertion<Complex, void, std::integral_constant<std::size_t,level>>
{
    static void apply(Complex& F, typename jbm::SimplexDataSet<Complex>::type& data)
    {
        std::cout << level << std::endl;
        for(auto curr : std::get<level>(data))
        {
//            std::cout << curr << std::endl;
            F.template insert<level>(curr);
        }
        PerformInsertion<Complex,typename Complex::template NodeData<level+1>,std::integral_constant<std::size_t,level+1>>::apply(F,data);
    }
};

template <typename Complex, typename NodeDataType>
struct PerformInsertion<Complex, NodeDataType, std::integral_constant<std::size_t,Complex::topLevel>>
{
    static void apply(Complex& F, typename jbm::SimplexDataSet<Complex>::type& data)
    {
        std::cout << Complex::topLevel << std::endl;
        for(auto curr : std::get<Complex::topLevel>(data))
        {
//            std::cout << curr.first << std::endl;
            F.template insert<Complex::topLevel>(curr.first, curr.second);
        }
    }
};

template <typename Complex, std::size_t k>
struct PerformRemoval
{
    static void apply(Complex& F, typename jbm::SimplexSet<Complex>::type& data)
    {
        for(auto curr : std::get<k>(data))
        {
            std::array<typename Complex::KeyType,k> name(F.get_name(curr));
            std::cout << "---------|| " << name << " ||---------" << std::endl;
            F.template remove<k>(curr);
        }
        PerformRemoval<Complex,k-1>::apply(F,data);
    }
};

template <typename Complex>
struct PerformRemoval<Complex,0>
{
    static void apply(Complex& F, typename jbm::SimplexSet<Complex>::type& data) {}
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
    
    std::vector<Callback<SurfaceMesh>> callbacks;
    int np = 100;
    int i = 0;
    typename jbm::SimplexSet<SurfaceMesh>::type doomed;

    for(auto s : mesh->get_level_id<3>())
    {
        if(i++ == 15)
        {
            Callback<SurfaceMesh> clbk;
            std::cout << mesh->get_name(s) << std::endl;
            typename jbm::SimplexSet<SurfaceMesh>::type levels;
            visit_node_down(TestBVisitor<SurfaceMesh>(&levels), *mesh, s);
            doomed = levels;
            visit_node_down(MainVisitor<SurfaceMesh,Callback>(&levels,&clbk,np), *mesh, s);
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
            callbacks.push_back(clbk);
            ++np;
        }
    }

    PerformRemoval<SurfaceMesh,3>::apply(*mesh, doomed);

    for(auto curr : callbacks)
    {
        PerformInsertion<SurfaceMesh,Vertex,std::integral_constant<std::size_t,1>>::apply(*mesh, curr.data);
    }

    for(auto s : mesh->get_level_id<3>())
    {
        std::cout << mesh->get_name(s) << std::endl;
    }

//    mesh->renumber();


    for(auto s : mesh->get_level_id<3>())
    {
        std::cout << mesh->get_name(s) << std::endl;
    }

    writeOFF("awesome.off", *mesh);




    std::map<int,int> fn;
    fn[-1] = 2;
    fn[-2] = 3;
    fn[-3] = 5;

    for(auto curr : fn)
    {
        std::cout << curr.first << std::endl;
    }
}
