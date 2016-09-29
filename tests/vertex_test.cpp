#include <iostream>
#include <map>
#include <cmath>
#include "Vertex.h"
#include "gtest/gtest.h"

TEST(VertexTests,Construction)
{
    Vertex v0 = Vertex(0,0,0);
    Vertex v1 = Vertex(0,0,0,0,false);
    Vertex v2 = Vertex(0,0,0,1,false);
    Vertex v3 = Vertex(0.1,0.2,0.3,4,true);

    EXPECT_EQ(0, v0[0]);
    EXPECT_EQ(0, v0[1]);
    EXPECT_EQ(0, v0[2]);
    EXPECT_EQ(0, v0.marker);

    EXPECT_EQ(v1, v0);
    EXPECT_NE(v2, v1);
}

TEST(VertexTests,MathOperations)
{
    Vertex v0 = Vertex(1.0,1.0,1.0);
    Vertex v1 = Vertex(2.1,3.2,4.3);
    
    EXPECT_EQ(Vertex(3.1,4.2,5.3), v0+v1); // check addition
    EXPECT_EQ(v1+v0, v0+v1); // check + commutative
    
    EXPECT_EQ(Vertex(-1.1,-2.2,-3.3), v0-v1);
    EXPECT_EQ(Vertex(1.1,2.2,3.3), v1-v0);
}
