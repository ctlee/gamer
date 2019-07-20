
#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <array>
#include <memory>
#include "gamer/SurfaceMesh.h"
#include "gtest/gtest.h"

/// Namespace for all things gamer
namespace gamer
{


class SurfaceMeshTest : public testing::Test {
protected:
	SurfaceMeshTest() {}
	~SurfaceMeshTest() {}
	virtual void SetUp() {
		mesh = std::make_unique<SurfaceMesh>();
		// Idealized icosahedron!
		std::vector<Vertex> vectors;
		vectors.push_back(Vertex(0.0, 0.0, 2.0));
		vectors.push_back(Vertex(1.788854, 0.000000, 0.894427));
  		vectors.push_back(Vertex(1.552786, 1.701302, 0.894427));
 		vectors.push_back(Vertex(-1.447214, 1.051462, 0.894427));
 		vectors.push_back(Vertex(-1.447214, -1.051462, 0.894427));
  		vectors.push_back(Vertex(0.552786, -1.701302, 0.894427));
  		vectors.push_back(Vertex(1.447214, 1.051462, -0.894427));
 		vectors.push_back(Vertex(-0.552786, 1.701302, -0.894427));
 		vectors.push_back(Vertex(-1.788855, 0.000000, -0.894427));
 		vectors.push_back(Vertex(-0.552786, -1.701302, -0.894427));
  		vectors.push_back(Vertex(1.447214, -1.051462, -0.894427));
  		vectors.push_back(Vertex(0.0, 0.0, -2.0));

  		int i = 0;
  		for(auto vector : vectors){
  			mesh->insert({++i}, vector);
  		}

  		mesh->insert({0,1,2});
  		mesh->insert({0,2,3});
  		mesh->insert({4,0,3});
		mesh->insert({5,0,4});
		mesh->insert({1,0,5});
		mesh->insert({2,1,6});
		mesh->insert({7,2,6});
		mesh->insert({3,2,7});
		mesh->insert({8,3,7});
		mesh->insert({4,3,8});
		mesh->insert({9,4,8});
		mesh->insert({5,4,9});
		mesh->insert({10,5,9});
		mesh->insert({6,1,10});
		mesh->insert({1,5,10});
		mesh->insert({6,11,7});
		mesh->insert({7,11,8});
		mesh->insert({8,11,9});
		mesh->insert({9,11,10});
		mesh->insert({10,11,6});
	}
	virtual void TearDown() {}

    std::unique_ptr<SurfaceMesh> mesh;
};

TEST_F(SurfaceMeshTest, Refinement){
    int fbefore = mesh->size<3>();
    mesh = refineMesh(*mesh);
    int fafter = mesh->size<3>();
    EXPECT_EQ(fbefore*4, fafter);
}

TEST_F(SurfaceMeshTest, FillHoles){
    int vbefore = mesh->size<1>();
    int ebefore = mesh->size<2>();
    int fbefore = mesh->size<3>();

    EXPECT_EQ(vbefore, 42);
    EXPECT_EQ(ebefore, 120);
    EXPECT_EQ(fbefore, 80);

    mesh->remove({0});
    mesh->remove({8});
    fillHoles(*mesh);
    int vafter = mesh->size<1>();
    int eafter = mesh->size<2>();
    int fafter = mesh->size<3>();

    EXPECT_EQ(vbefore-2, vafter);
    EXPECT_EQ(ebefore-6, eafter);
    EXPECT_EQ(fbefore-4, fafter);
}

TEST_F(SurfaceMeshTest, Smooth){
    int vbefore = mesh->size<1>();
    int ebefore = mesh->size<2>();
    int fbefore = mesh->size<3>();
    smoothMesh(*mesh, 3, true, false);

    EXPECT_EQ(vbefore, 42);
    EXPECT_EQ(ebefore, 120);
    EXPECT_EQ(fbefore, 80);
}

} // end namespace gamer