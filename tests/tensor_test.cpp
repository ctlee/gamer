
#include <iostream>
#include <map>
#include <cmath>
#include "tensor.h"
#include "SurfaceMesh.h"
#include "kdtree.h"

int main(int argc, char *argv[])
{
    tensor<double, 4, 1> v0;
    tensor<double, 4, 1> v1;
    tensor<double, 4, 1> v2;
    tensor<double, 4, 1> v3;
/*
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;

    for(auto& curr : v0)
    {
        curr = distribution(generator);
    }
    for(auto& curr : v1)
    {
        curr = distribution(generator);
    }
*/
    v0[{0}] = 1; v0[{1}] = 1; v0[{2}] = 0; v0[{3}] = 0;
    v1[{0}] = 1; v1[{1}] = 1; v1[{2}] = 0; v1[{3}] = 0;
    v2[{0}] = 0; v2[{1}] = 0; v2[{2}] = 1; v2[{3}] = 0;
    v3[{0}] = 1; v3[{1}] = 0; v3[{2}] = 0; v3[{3}] = 2;

    auto uu = v0 * v1;
    auto w  = v0 ^ v1 ^ v2 ^ v3;

    std::cout << uu.get(0,0) << std::endl;

    int m = 0;
    for(auto curr = uu.index_begin(); curr != uu.index_end() && m < 100; ++curr, ++m)
    {
        for(auto x : (*curr))
        {
            std::cout << x << " ";
        }
        std::cout << ": " << uu.get(*curr) << std::endl;
    }

    for(auto curr = v0.index_begin(); curr != v0.index_end(); ++curr)
    {
        std::cout << ": " << uu.get(*curr,*curr) << std::endl;
    }

    std::cout << std::sqrt(dot(uu,uu)) << std::endl;
    std::cout << "Test Complete" << std::endl;

	std::cout << "Test complete" << std::endl;
	return -1;
}
