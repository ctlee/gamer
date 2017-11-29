/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "stringutil.h"
#include "SurfaceMeshOld.h"
#include "SurfaceMesh.h"
#include "Vertex.h"
#include "tensor.h"

#include <libraries/casc/include/CASCFunctions.h>
#include <libraries/casc/include/SimplexSet.h>
#include <libraries/casc/include/SimplexMap.h>
#include <libraries/casc/include/decimate.h>
#include <libraries/casc/include/typetraits.h>

void CardanosMethod(double A[3][3]){
    double c0, c1, c2;
    double x1, x2, x3;
    double a, b, Q;
    double theta, p;
    double B[6];

    c0 = A[0][0] * A[1][1] * A[2][2] + 2 * A[0][1] * A[0][2] * A[1][2] - A[0][0] * A[1][2] * A[1][2]
         - A[1][1] * A[0][2] * A[0][2] - A[2][2] * A[0][1] * A[0][1];
    c1 = A[0][0] * A[1][1] - A[0][1] * A[0][1] + A[0][0] * A[2][2] -
         A[0][2] * A[0][2] + A[1][1] * A[2][2] - A[1][2] * A[1][2];
    c2 = A[0][0] + A[1][1] + A[2][2];

    a = (3.0 * c1 - c2 * c2) / 3.0;
    b = (-2.0 * c2 * c2 * c2 + 9.0 * c1 * c2 - 27.0 * c0) / 27.0;
    Q = b * b / 4.0 + a * a * a / 27.0;

    theta = atan2(sqrt(-Q), -0.5 * b);
    p     = sqrt(0.25 * b * b - Q);

    x1 = c2 / 3.0 + 2.0 * pow(p, 1.0 / 3.0) * cos(theta / 3.0);
    x2 = c2 / 3.0 - pow(p, 1.0 / 3.0) * (cos(theta / 3.0) + sqrt(3.0) * sin(theta / 3.0));
    x3 = c2 / 3.0 - pow(p, 1.0 / 3.0) * (cos(theta / 3.0) - sqrt(3.0) * sin(theta / 3.0));


    double tx, ty, tz;

    // TODO: use a vector to represent these... and call SORT!
    tx = std::max(x1, std::max(x2, x3));

    if (tx == x1)
    {
        if (x2 >= x3)
        {
            ty = x2;
            tz = x3;
        }
        else
        {
            ty = x3;
            tz = x2;
        }
    }
    else if (tx == x2)
    {
        if (x1 >= x3)
        {
            ty = x1;
            tz = x3;
        }
        else
        {
            ty = x3;
            tz = x1;
        }
    }
    else // if (tx == x3)
    {
        if (x1 >= x2)
        {
            ty = x1;
            tz = x2;
        }
        else
        {
            ty = x2;
            tz = x1;
        }
    }
    double evx, evy, evz;

    x1             = tx;
    x2             = ty;
    x3             = tz;
    evx = tx;
    evy = ty;
    evz = tz;

    std::cout << "Old Eigenvalues: " << evx << " " << evy << " " << evz << std::endl;
}


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

template <typename Complex>
struct Callback
{
    using SimplexSet = typename casc::SimplexSet<Complex>;
    using KeyType = typename Complex::KeyType;

    template <std::size_t k>
    void operator()(Complex& F,
            const std::array<KeyType, k>& new_name,
            const SimplexSet& merged){
        std::cout << merged << " -> " << new_name << std::endl;
    }

    Vertex operator()(Complex& F,
            const std::array<KeyType, 1>& new_name,
            const SimplexSet& merged){
        std::cout << merged << " -> " << new_name << std::endl;

        Vertex center;
        std::size_t cnt = 0;
        for(auto v : merged.template get<1>())
        {
            center = center + (*v);
            ++cnt;
        }
        center = center / (double)(cnt);
        return center;
    }

    Face operator()(Complex& F,
            const std::array<KeyType, 3>& new_name,
            const SimplexSet& merged){
        std::cout << merged << " -> " << new_name << std::endl;
        return Face();
    }
};


int  main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    auto mesh = readOFF(argv[1]);
    if(mesh == nullptr){
        std::cout << "Something bad happened..." << std::endl;
        exit(1);
    }
    std::cout << "Done reading" << std::endl;

    compute_orientation(*mesh);

    double volume = getVolume(*mesh);
    std::cout << "Volume: " << volume << std::endl;

    //coarseIT(*mesh, 0.0160, 1, 0);
    coarseIT(*mesh, 1, 0, 10);
    writeOFF("test.off", *mesh);

    // auto oldmesh = SurfaceMeshOld::readOFF(argv[1]);
    // oldmesh->coarse(0.016,1,0,-1);

    // compute_orientation(*mesh);
   

    // auto it = mesh->get_level_id<1>().begin();
    // auto nid = *it;

    // std::cout << "GetTangent of: " << nid << std::endl;

    // auto tangent = getTangent(*mesh, nid);
    // std::cout << "Tangent: " << tangent << std::endl;

    // auto o = (*mesh->get_simplex_up({0})).position;
    // auto a = (*mesh->get_simplex_up({1})).position - o;
    // auto b = (*mesh->get_simplex_up({2})).position - o;
    // auto c = (*mesh->get_simplex_up({3})).position - o;

    // std::cout << ((c^b) + (b^a) + (a^c))/6 << std::endl;
    // std::cout << ((b^c) + (a^b) + (c^a))/6 << std::endl;


    // auto res = (-1*b*((c - a)/2) - c*((-1*b+a)/2) -a*((b-c)/2))/3;
    // //auto res = ((-1*b*c + b*a)/2 + (c*b - c*a)/2 -(a*b - a*c)/2)/3;
    // //auto res = (-1*b*c/2 + b*a/2 + c*b/2 - c*a/2 -a*b/2 + a*c/2)/3;
    // //auto res = (-1*b*c/2 + b*a/2 + c*b/2 - c*a/2 -a*b/2 + a*c/2)/3;
    // std::cout << res << std::endl << std::endl;


    // auto a = tensor<double, 3,1>({1,0,0});
    // auto b = tensor<double, 3,1>({0,1,0});

    // auto res = a^b;

    // std::cout << (a^b) << std::endl;
    // std::cout << (a*b) - (b*a) << std::endl;
    // std::cout << std::sqrt(res|res) << std::endl;

    // using SimplexSet = typename casc::SimplexSet<SurfaceMesh>;
    // SimplexSet S;

    // auto it = mesh->get_level_id<2>().begin();

    // auto nid = *it;
    // S.insert(nid);
    // nid = *(++it);
    // S.insert(nid);

    // SimplexSet dest;
    // std::cout << "Star of: " << S << std::endl;
    // getStar(*mesh, S, dest);
    // std::cout << ">> " << dest << std::endl << std::endl;

    // dest.clear();
    // std::cout << "Closure of: " << S << std::endl;
    // getClosure(*mesh, S, dest);
    // std::cout << ">> " << dest << std::endl << std::endl;

    // dest.clear();
    // S.clear(); 
    // S.insert(*(mesh->get_level_id<1>().begin()));
    // std::cout << "Link of: " << S << std::endl;
    // getLink(*mesh, S, dest);
    // std::cout << ">> " << dest << std::endl << std::endl; 


    /*
     * TIMING CODE STARTS HERE
     */

    // int trials = 10;
    // std::chrono::duration<double> elapsed_seconds;

    // for(int trial=0; trial < trials; ++trial){
    //     std::chrono::time_point<std::chrono::system_clock> start, end;
    //     start = std::chrono::system_clock::now();
    //     // CODE GOES HERE

    //     auto mesh = SurfaceMeshOld::readOFF(argv[1]);
    //     for(int i=0; i < 10; ++i){
    //         mesh->normalSmooth();
    //     }

    //     // CODE ENDS
    //     end = std::chrono::system_clock::now();
    //     free(mesh);
    //     elapsed_seconds += end-start;
    // }
    // elapsed_seconds /= trials;
    // std::cout << "Average time old: " << elapsed_seconds.count() << "s\n";
    
    // // New code benchmark 
    // for(int trial=0; trial < trials; ++trial){
    //     std::chrono::time_point<std::chrono::system_clock> start, end;
    //     start = std::chrono::system_clock::now();
    //     // CODE GOES HERE

    //     auto mesh = readOFF(argv[1]);
    //     if(mesh == nullptr){
    //         std::cout << "Something bad happened...";
    //         exit(1);
    //     }

    //     compute_orientation(*mesh);

    //     for(auto i = 0; i < 1; ++i){
    //         for(auto nid : mesh->get_level_id<1>())
    //         {
    //             normalSmooth(*mesh, nid);
    //         }
    //     }

    //     // CODE ENDS
    //     end = std::chrono::system_clock::now();
    //     elapsed_seconds += end-start;
    // }
    // elapsed_seconds /= trials;
    // std::cout << "Average time new: " << elapsed_seconds.count() << "s\n";


    std::cout << "EOF" << std::endl;
}
