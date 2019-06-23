
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
        // Idealized icosphere!
        std::vector<SMVertex> vectors;
        vectors.push_back(SMVertex(0, 0, -4.70864));
        vectors.push_back(SMVertex(3.40721, -2.47545, -2.1058));
        vectors.push_back(SMVertex(-1.30141, -4.0054, -2.1058));
        vectors.push_back(SMVertex(-4.21153, 0, -2.10578));
        vectors.push_back(SMVertex(-1.30141, 4.0054, -2.1058));
        vectors.push_back(SMVertex(3.40721, 2.47545, -2.1058));
        vectors.push_back(SMVertex(1.30141, -4.0054, 2.1058));
        vectors.push_back(SMVertex(-3.40721, -2.47545, 2.1058));
        vectors.push_back(SMVertex(-3.40721, 2.47545, 2.1058));
        vectors.push_back(SMVertex(1.30141, 4.0054, 2.1058));
        vectors.push_back(SMVertex(4.21153, 0, 2.10578));
        vectors.push_back(SMVertex(0, 0, 4.70864));
        vectors.push_back(SMVertex(-0.764945, -2.3543, -4.00542));
        vectors.push_back(SMVertex(2.00269, -1.45502, -4.00542));
        vectors.push_back(SMVertex(1.23775, -3.80934, -2.47551));
        vectors.push_back(SMVertex(4.00539, 0, -2.4755));
        vectors.push_back(SMVertex(2.00269, 1.45502, -4.00542));
        vectors.push_back(SMVertex(-2.47547, 0, -4.00541));
        vectors.push_back(SMVertex(-3.24044, -2.35431, -2.4755));
        vectors.push_back(SMVertex(-0.764945, 2.3543, -4.00542));
        vectors.push_back(SMVertex(-3.24044, 2.35431, -2.4755));
        vectors.push_back(SMVertex(1.23775, 3.80934, -2.47551));
        vectors.push_back(SMVertex(4.47819, -1.45503, 0));
        vectors.push_back(SMVertex(4.47819, 1.45503, 0));
        vectors.push_back(SMVertex(0, -4.70864, 0));
        vectors.push_back(SMVertex(2.76767, -3.80937, 0));
        vectors.push_back(SMVertex(-4.47819, -1.45503, 0));
        vectors.push_back(SMVertex(-2.76767, -3.80937, 0));
        vectors.push_back(SMVertex(-2.76767, 3.80937, 0));
        vectors.push_back(SMVertex(-4.47819, 1.45503, 0));
        vectors.push_back(SMVertex(2.76767, 3.80937, 0));
        vectors.push_back(SMVertex(0, 4.70864, 0));
        vectors.push_back(SMVertex(3.24044, -2.35431, 2.4755));
        vectors.push_back(SMVertex(-1.23775, -3.80934, 2.47551));
        vectors.push_back(SMVertex(-4.00539, 0, 2.4755));
        vectors.push_back(SMVertex(-1.23775, 3.80934, 2.47551));
        vectors.push_back(SMVertex(3.24044, 2.35431, 2.4755));
        vectors.push_back(SMVertex(0.764945, -2.3543, 4.00542));
        vectors.push_back(SMVertex(2.47547, 0, 4.00541));
        vectors.push_back(SMVertex(-2.00269, -1.45502, 4.00542));
        vectors.push_back(SMVertex(-2.00269, 1.45502, 4.00542));
        vectors.push_back(SMVertex(0.764945, 2.3543, 4.00542));
        int i = 0;
        for(auto vector : vectors){
            mesh->insert({i++}, vector);
        }

        mesh->insert({0,12,13});
        mesh->insert({1,13,15});
        mesh->insert({0,12,17});
        mesh->insert({0,17,19});
        mesh->insert({0,16,19});
        mesh->insert({1,15,22});
        mesh->insert({2,14,24});
        mesh->insert({3,18,26});
        mesh->insert({4,20,28});
        mesh->insert({5,21,30});
        mesh->insert({1,22,25});
        mesh->insert({2,24,27});
        mesh->insert({3,26,29});
        mesh->insert({4,28,31});
        mesh->insert({5,23,30});
        mesh->insert({6,32,37});
        mesh->insert({7,33,39});
        mesh->insert({8,34,40});
        mesh->insert({9,35,41});
        mesh->insert({10,36,38});
        mesh->insert({11,38,41});
        mesh->insert({36,38,41});
        mesh->insert({9,36,41});
        mesh->insert({11,40,41});
        mesh->insert({35,40,41});
        mesh->insert({8,35,40});
        mesh->insert({11,39,40});
        mesh->insert({34,39,40});
        mesh->insert({7,34,39});
        mesh->insert({11,37,39});
        mesh->insert({33,37,39});
        mesh->insert({6,33,37});
        mesh->insert({11,37,38});
        mesh->insert({32,37,38});
        mesh->insert({10,32,38});
        mesh->insert({10,23,36});
        mesh->insert({23,30,36});
        mesh->insert({9,30,36});
        mesh->insert({9,31,35});
        mesh->insert({28,31,35});
        mesh->insert({8,28,35});
        mesh->insert({8,29,34});
        mesh->insert({26,29,34});
        mesh->insert({7,26,34});
        mesh->insert({7,27,33});
        mesh->insert({24,27,33});
        mesh->insert({6,24,33});
        mesh->insert({6,25,32});
        mesh->insert({22,25,32});
        mesh->insert({10,22,32});
        mesh->insert({9,30,31});
        mesh->insert({21,30,31});
        mesh->insert({4,21,31});
        mesh->insert({8,28,29});
        mesh->insert({20,28,29});
        mesh->insert({3,20,29});
        mesh->insert({7,26,27});
        mesh->insert({18,26,27});
        mesh->insert({2,18,27});
        mesh->insert({6,24,25});
        mesh->insert({14,24,25});
        mesh->insert({1,14,25});
        mesh->insert({10,22,23});
        mesh->insert({15,22,23});
        mesh->insert({5,15,23});
        mesh->insert({5,16,21});
        mesh->insert({16,19,21});
        mesh->insert({4,19,21});
        mesh->insert({4,19,20});
        mesh->insert({17,19,20});
        mesh->insert({3,17,20});
        mesh->insert({3,17,18});
        mesh->insert({12,17,18});
        mesh->insert({2,12,18});
        mesh->insert({5,15,16});
        mesh->insert({13,15,16});
        mesh->insert({0,13,16});
        mesh->insert({2,12,14});
        mesh->insert({12,13,14});
        mesh->insert({1,13,14});
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

} // end namespace gamer