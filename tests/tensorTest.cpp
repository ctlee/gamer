
#include <iostream>
#include <map>
#include <cmath>
#include "tensor.h"
#include "gtest/gtest.h"


TEST(TensorTest,MathOperations){
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
}