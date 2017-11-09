#include <iostream>
#include <map>
#include <cmath>
#include <random>
#include "Vertex.h"
#include "gtest/gtest.h"


// Tests of the constructor
TEST(VertexTest,TestConstructor)
{
    Vertex v0 = Vertex(0,0,0);
    Vertex v1 = Vertex(0,0,0,0,false);
    Vertex v2 = Vertex(0,0,0,1,false);
    Vertex v3 = Vertex(0,0,0,0,true); 
    Vertex v4 = Vertex();

    EXPECT_EQ(0, v0[0]);
    EXPECT_EQ(0, v0[1]);
    EXPECT_EQ(0, v0[2]);
    EXPECT_EQ(0, v0.marker);
    EXPECT_TRUE(false == v0.selected);
    EXPECT_TRUE(true == v3.selected);
    
    EXPECT_EQ(v1, v0);
    
    EXPECT_NE(v1, v2);
    EXPECT_NE(v2, v1);
    EXPECT_NE(v1, v3);
    EXPECT_NE(v3, v1);
    
    EXPECT_EQ(v0, v4);
   
    Vertex v5 = Vertex(v2);
    EXPECT_EQ(v2, v5);

    v3.selected = false;
    EXPECT_EQ(v0,v3);
}

TEST(VertexTest, MathOperations)
{
    Vertex v0 = Vertex(1.0,1.0,1.0);
    Vector vec = (Vector) Vertex(2.1,3.2,4.3); // TODO: change this so the cast isn't weird
    
    EXPECT_EQ(Vertex(3.1,4.2,5.3), v0+vec); // check addition

    EXPECT_EQ(Vertex(-1.1,-2.2,-3.3), v0-vec);
    EXPECT_EQ(Vertex(1,1,1), 1*v0);
    EXPECT_EQ(Vertex(3.4,3.4,3.4), 3.4*v0);
    EXPECT_EQ(Vertex(1.0,1.0,1.0), v0/1.0);
}

TEST(VertexTest, TypeCastOperations)
{
    // TODO: implement some tests to ensure that return types are correct after each operation
}

TEST(VertexTest, FunctionOperations){
    // TODO: check the accuracy of different functions
}

#if GTEST_HAS_PARAM_TEST
using ::testing::TestWithParam;
using ::testing::Range;

class VertexTest : public TestWithParam<int> {
public:
    VertexTest() {}
    ~VertexTest() {}
    virtual void SetUp() {
		std::random_device rd;
    	std::mt19937 gen(rd());
    	std::uniform_real_distribution<double> dis(-10000, 10000);
		for(int i=0; i<3; i++) val0[i] = dis(gen); 	
		for(int i=0; i<3; i++) val1[i] = dis(gen); 	
		v0 = Vertex(val0[0], val0[1], val0[2]);
		v1 = Vector();
        v1[0] = val1[0]; v1[1] = val1[1]; v1[2] = val1[2];
	}
    virtual void TearDown() {}

protected:	
	Vertex v0;
    Vector v1;
	double val0[3], val1[3];
};


TEST_P(VertexTest,RandomDoubleAddition)
{
	EXPECT_EQ(Vertex(val0[0] + val1[0], 
				val0[1] + val1[1],
				val0[2] + val1[2]), v0+v1);
    EXPECT_EQ(Vertex(val0[0] - val1[0], 
                val0[1] - val1[1],
                val0[2] - val1[2]), v0-v1);
}

// Run the randomdoubles 100 times
INSTANTIATE_TEST_CASE_P(RandomDoubles,
                        VertexTest,
                        Range(1,10));

#else
TEST(DummyTest, ValueParameterizedTestsAreNotSupportedOnThisPlatform) {}
#endif // GTEST_HAS_PARAM_TEST
