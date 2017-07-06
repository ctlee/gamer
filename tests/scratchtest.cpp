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

#include <libraries/casc/include/SimplexSets.h>

namespace test{
    template <typename Integer, typename IntegerSequence, typename Fn>
    struct int_for_each_helper {};

    template <class Integer, template <class, Integer...> class InHolder, Integer I, typename Fn>
    struct int_for_each_helper<Integer, InHolder<Integer, I>, Fn>
    {
        static void apply(Fn&& f){
            static constexpr auto k = I;
            f.template apply<k>();
            std::cout << I << std::endl;
        }
    };

    template <class Integer, template <class, Integer...> class InHolder, Integer I, Integer... Is, typename Fn>
    struct int_for_each_helper<Integer, InHolder<Integer, I, Is...>, Fn> 
    {
        static void apply(Fn&& f){
            f.template apply<I>();
            std::cout << I << std::endl;
            int_for_each_helper<Integer, InHolder<Integer, Is...>, Fn>::apply(std::forward<Fn>(f));
        }
    };



    template <class Integer, typename IntegerSequence, typename Fn>
    void int_for_each(Fn&& f){
            int_for_each_helper<Integer, IntegerSequence, Fn>::apply(std::forward<Fn>(f));
    }

    struct Foo
    {
        template <std::size_t k>
        void apply(){
            std::cout << "foo" <<  k << std::endl;
        }
    };
}

int  main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }
    auto mesh = readOFF(argv[1]);
    if(mesh == nullptr){
        std::cout << "Something bad happened...";
        exit(1);
    }
    std::cout << "Done reading" << std::endl;

    std::cout << "Sizes: " << mesh->size<1>() << " " << mesh->size<2>() << " " << mesh->size<3>() << std::endl;

    compute_orientation(*mesh);
    double volume = getVolume(*mesh);
   
    std::cout << "Volume: " << volume << std::endl;

    using SimplexSet = typename casc::SimplexSet<SurfaceMesh>;
    SimplexSet A;
    SimplexSet B;

    for(auto nid : mesh->get_level_id<1>())
    {
        std::cout << nid << std::endl;
        A.insert(nid);
    }

    casc::printSS<SurfaceMesh>(A);




    //SimplexSet C = casc::set_union(A, B);
    //A.print();

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