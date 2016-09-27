
#include <iostream>
#include <map>
#include <cmath>
#include "tensor.h"

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
    v0.get(0) = 1; v0.get(1) = 1; v0.get(2) = 0; v0.get(3) = 0;
    v1.get(0) = 1; v1.get(1) = 1; v1.get(2) = 0; v1.get(3) = 0;
    v2.get(0) = 0; v2.get(1) = 0; v2.get(2) = 1; v2.get(3) = 0;
    v3.get(0) = 1; v3.get(1) = 0; v3.get(2) = 0; v3.get(3) = 2;

    auto uu = v0 * v0;
    auto w  = v0 ^ v1 ^ v2 ^ v3;

    for(auto x : w)
    {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    for(auto curr = uu.index_begin(); curr != uu.index_end(); ++curr)
    {
        for(auto x : (*curr))
        {
            std::cout << x << " ";
        }
        std::cout << ":: " << uu[*curr] << std::endl;
    }

    for(auto curr = v0.index_begin(); curr != v0.index_end(); ++curr)
    {
        std::cout << ": " << uu.get(*curr,*curr) << std::endl;
    }

    std::cout << std::sqrt(dot(w,w)) << std::endl;
    std::cout << "Test Complete" << std::endl;
	return -1;
}
