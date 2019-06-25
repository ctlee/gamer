
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
        mesh = sphere(0);
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