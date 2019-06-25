
#include <iostream>
#include <map>
#include <cmath>
#include "gamer/tensor.h"
#include "gtest/gtest.h"

/// Namespace for all things gamer
namespace gamer
{


TEST(TensorTest, MathOperations){
    tensor<double, 3, 1> v0;
    tensor<double, 3, 1> v1;
    tensor<double, 3, 1> v2;

    v0[0] = 1; v0[1] = 0; v0[2] = 0;
    v1[0] = 0; v1[1] = 1; v1[2] = 0;
    v2[0] = 0; v2[1] = 0; v2[2] = 1;

    auto v0p1 = v0 + v1;
    auto v1p0 = v1 + v0;
    EXPECT_EQ(v0p1,v1p0);

    auto v01 = v0 ^ v1;
    auto v10 = v1 ^ v0;
    EXPECT_EQ(v01, -v10);

    tensor<double,3,2> v3({1,2,3,4,5,6,7,8,9});
    tensor<double,3,2> v4({1,3,5,3,5,7,5,7,9});
    EXPECT_EQ(Sym(v3), v4);

    // TODO: (0) add test for elementwise division etc...
}

TEST(TensorTest, TensorIndexIterator){

    tensor<short, 3, 4> t34;

    std::array<std::size_t, 4> e34i = {3,0,0,0};
    std::array<std::size_t, 4> b34i = {0,0,0,0};
    EXPECT_EQ(e34i, *(t34.index_end()));
    EXPECT_EQ(b34i, *(t34.index_begin()));

    tensor<short, 5, 3> t53;
    std::array<std::size_t, 3> e53i = {5,0,0};
    std::array<std::size_t, 3> b53i = {0,0,0};
    EXPECT_EQ(e53i, *(t53.index_end()));
    EXPECT_EQ(b53i, *(t53.index_begin()));

    for (auto it = t34.index_begin(); it != t34.index_end(); ++it){
        for(std::size_t k = 0; k < 4; ++k){
            ASSERT_LT((*it)[k], 3);
        }
    }

    for (auto it = --t34.index_end(); it != t34.index_begin(); --it){
        for(std::size_t k = 0; k < 4; ++k){
            ASSERT_LT((*it)[k], 3);
        }
    }

    std::size_t tot = 0;
    std::size_t tdim = tensor<short, 3, 4>::total_dimension;
    for (auto it = t34.index_begin(); it != t34.index_end(); ++it)
        ++tot;
    EXPECT_EQ(tot, tdim);
    tot = 0;
    for (auto it = t34.index_end(); it != t34.index_begin(); it--)
        ++tot;
    EXPECT_EQ(tot, tdim);

    tot = 0;
    tdim = tensor<short, 5, 3>::total_dimension;
    for (auto it = t53.index_begin(); it != t53.index_end(); it++)
        ++tot;
    EXPECT_EQ(tot, tdim);
    tot = 0;
    for (auto it = t53.index_end(); it != t53.index_begin(); it--)
        ++tot;
    EXPECT_EQ(tot, tdim);
}

TEST(TensorTest, Constructor){
    tensor<std::size_t, 3, 1> v0;

    for(std::size_t i = 0; i < 3; ++i){
        ASSERT_EQ(v0[i], 0);
    }

    tensor<std::size_t, 3, 2> v1({0,0,0, 0,1,2, 0,2,4});
    for(std::size_t k = 0; k < 3; ++k){
        for(std::size_t i = 0; i < 3; ++i){
            ASSERT_EQ(v1.get(i,k), i*k);
        }
    }

    tensor<double, 3, 1> v2({0.1,1.2,2.3});
    v0 = v2; // Implicit cast truncates double to size_t
    for(std::size_t i = 0; i < 3; ++i){
        ASSERT_EQ(v0[i], i);
    }

    tensor<std::size_t, 3, 3> v3(5);
    for (auto it = v3.index_begin(); it != v3.index_end(); it++){
        // Test three methods of index accession...
        ASSERT_EQ(v3.get(*it), 5);
        ASSERT_EQ(v3[*it], 5);
        ASSERT_EQ(v3.get((*it)[0],(*it)[1],(*it)[2]), 5);
    }
}
} // end namespace gamer