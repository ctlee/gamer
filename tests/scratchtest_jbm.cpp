/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include "SurfaceMesh.h"
#include "Vertex.h"
//#include "SimplicialComplexVisitors.h"
#include <vector>
#include <iostream>
#include <string>
#include <type_traits>
#include <list>
#include <libraries/casc/include/decimate.h>


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


template <typename Complex>
struct Callback
{
    using SimplexSet = typename casc::SimplexSet<Complex>::type;
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

        return center;
    }

    template <std::size_t OldLevel>
    Face operator()(const Complex& F,
                    const std::array<KeyType,OldLevel>& old_name,
                    const std::array<KeyType,3>& new_name,
                    const SimplexSet& merged)
    {
        std::cout << "3nner: " << old_name << " --> " << new_name << " : " << std::get<3>(merged).size() << std::endl;
        return Face();
    }
};


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
    
    DecimateExample::Callback<SurfaceMesh> clbk;
    for(int i = 0; i < 1; ++i)
    {
        auto s = *(mesh->get_level_id<2>().begin());
        decimate(*mesh, s, clbk);
    }

    compute_orientation(*mesh);

    writeOFF("awesome.off", *mesh);
}
