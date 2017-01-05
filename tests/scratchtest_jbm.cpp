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
auto getTangentH(SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::SimplexID<SurfaceMesh::topLevel> curr)
{
    return 1.0;
}

template <std::size_t level, std::size_t dimension>
auto getTangentH(SurfaceMesh& mesh, const tensor<double, dimension, 1>& origin, SurfaceMesh::SimplexID<level> curr)
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

auto getTangent_jbm(SurfaceMesh& mesh, SurfaceMesh::SimplexID<1> vertexID)
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

template <typename SimplexID>
void grab(const typename SimplexID::complex& F, SimplexID s)
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
    bool visit(const Complex& F, typename Complex::template SimplexID<level> s)
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
        using Simplex = typename Complex::template SimplexID<j>;
        using LevelIndex = typename std::make_index_sequence<Complex::numLevels>;
        using SimplexIDLevel = typename util::int_type_map<std::size_t, std::tuple, LevelIndex, Simplex>::type;
        using type = typename util::type_map<SimplexIDLevel, NodeSet>::type;
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
        using SimplexIDLevel = typename util::int_type_map<std::size_t, std::tuple, LevelIndex, DataSet>::type;
        using type = typename util::type_map<SimplexIDLevel, jbm::vector>::type;
    };
}

template <typename Complex>
struct GetCompleteNeighborhoodHelper
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;

    GetCompleteNeighborhoodHelper(SimplexSet* p) : pLevels(p) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template SimplexID<level> s)
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
struct GetCompleteNeighborhood
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;

    GetCompleteNeighborhood(SimplexSet* p) : pLevels(p) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template SimplexID<level> s)
    {
        return true;
    }

    bool visit(const Complex& F, typename Complex::template SimplexID<1> s)
    {
        visit_node_up(GetCompleteNeighborhoodHelper<Complex>(pLevels), F, s);
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
    bool visit(const Complex& F, typename Complex::template SimplexID<level> s)
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
    using Simplex = typename Complex::template SimplexID<BaseLevel>;
    using KeyType = typename Complex::KeyType;
    using ReturnValues = typename jbm::SimplexDataSet<Complex>::type;

    template <typename ReturnType, std::size_t OldLevel>
    struct PerformCallback
    {
        constexpr static std::size_t NewLevel = OldLevel - BaseLevel + 1;
        static void apply(Callback<Complex>* callback,
                    ReturnValues* data,
                    const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,NewLevel>& new_name,
                    const SimplexSet& merged)
        {
            ReturnType rval = (*callback)(F,old_name,new_name,merged);
            std::get<NewLevel>(*data).push_back(std::make_pair(new_name,rval));
        }
    };

    template <std::size_t OldLevel>
    struct PerformCallback<void,OldLevel>
    {
        constexpr static std::size_t NewLevel = OldLevel - BaseLevel + 1;
        static void apply(Callback<Complex>* callback,
                    ReturnValues* data,
                    const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,NewLevel>& new_name,
                    const SimplexSet& merged)
        {
            (*callback)(F,old_name,new_name,merged);
            std::get<NewLevel>(*data).push_back(new_name);
        }
    };

    InnerVisitor(SimplexSet* p, Simplex s, KeyType np, ReturnValues* rv, Callback<Complex>* c)
        : pLevels(p), simplex(s), new_point(np), callback(c), data(rv) {}

    template <std::size_t OldLevel>
    bool visit(const Complex& F, typename Complex::template SimplexID<OldLevel> s)
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
            PerformCallback<typename Complex::template NodeData<NewLevel>, OldLevel>::apply(callback,data,F, old_name, new_name, grab);
        }
        return true;
    }

private:
    SimplexSet* pLevels;
    Simplex     simplex;
    KeyType     new_point;
    Callback<Complex>* callback;
    ReturnValues* data;
};


template <typename Complex, template <typename> class Callback>
struct MainVisitor
{
    using SimplexSet = typename jbm::SimplexSet<Complex>::type;
    using KeyType = typename Complex::KeyType;
    using ReturnValues = typename jbm::SimplexDataSet<Complex>::type;

    MainVisitor(SimplexSet* p, Callback<Complex>* c, KeyType np, ReturnValues* rv)
        : pLevels(p), callback(c), new_point(np), data(rv) {}

    template <std::size_t level>
    bool visit(const Complex& F, typename Complex::template SimplexID<level> s)
    {
//        std::cout << "MAIN: " << F.get_name(s) << std::endl;
        visit_node_up(InnerVisitor<Complex,level,Callback>(pLevels,s,new_point,data,callback), F, s);
        return true;
    }

private:
    SimplexSet* pLevels;
    Callback<Complex>* callback;
    KeyType     new_point;
    ReturnValues* data;
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
    Vertex operator()(const Complex& F,
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

//        std::get<1>(data).push_back(std::make_pair(new_name, center));
        return center;
    }

    template <std::size_t OldLevel>
    Face operator()(const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,3>& new_name,
                    const SimplexSet& merged)
    {
        std::cout << "3nner: " << old_name << " --> " << new_name << " : " << std::get<3>(merged).size() << std::endl;
//        std::get<3>(data).push_back(std::make_pair(new_name,Face()));
        return Face();
    }
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
        for(auto curr : std::get<level>(data))
        {
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
        for(auto curr : std::get<Complex::topLevel>(data))
        {
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
//            std::cout << "---------|| " << F.get_name(curr) << " ||---------" << std::endl;
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


template <typename Complex, typename Simplex, template<typename> class Callback>
void decimate(Complex& F, Simplex s, Callback<Complex>& clbk)
{
    static int np = -1;
    int i = 0;
    typename jbm::SimplexDataSet<SurfaceMesh>::type rv;
    typename jbm::SimplexSet<SurfaceMesh>::type levels;

    visit_node_down(GetCompleteNeighborhood<SurfaceMesh>(&levels), F, s);
    typename jbm::SimplexSet<SurfaceMesh>::type doomed = levels;
    visit_node_down(MainVisitor<SurfaceMesh,Callback>(&levels,&clbk,np,&rv), F, s);

    PerformRemoval<SurfaceMesh,3>::apply(F, doomed);
    PerformInsertion<SurfaceMesh,Vertex,std::integral_constant<std::size_t,1>>::apply(F, rv);
    --np;
}

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    std::cout << "Begin reading Mesh..." << std::endl;
    auto mesh = readOFF(argv[1]);
    if(mesh == nullptr){
        std::cout << "Something bad happened...";
        exit(1);
    }
    
    compute_orientation(*mesh);
    
    int np = 100;
    typename jbm::SimplexSet<SurfaceMesh>::type doomed;
    std::vector<typename jbm::SimplexDataSet<SurfaceMesh>::type> rvals;

    Callback<SurfaceMesh> clbk;
    for(int i = 0; i < 45000; ++i)
    {
        auto s = *(mesh->get_level_id<2>().begin());
        decimate(*mesh, s, clbk);
    }

    writeOFF("awesome.off", *mesh);
}
