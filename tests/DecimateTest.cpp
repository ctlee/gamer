//
// Created by ishan on 12/3/18.
//



#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <array>
#include <memory>
#include "TetMesh.h"
#include "Vertex.h"
#include "gtest/gtest.h"

template <typename Complex>
struct Callback
{
    using SimplexSet = typename casc::SimplexSet<Complex>;
    using KeyType = typename Complex::KeyType;

    template <std::size_t k>
    int operator()(Complex& F,
                   const std::array<KeyType, k>& new_name,
                   const SimplexSet& merged){
            // std::cout << merged << " -> " << new_name << std::endl;
            return 0;
    }
};

class DecimateTest : public testing::Test {
protected:
    DecimateTest() {}
    ~DecimateTest() {}
    virtual void SetUp() {
        mesh = new SurfaceMesh();
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
        vectors.push_back(Vertex(-1.788854, 0.000000, -0.894427));
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

        std::vector<SurfaceMesh*> vec;
        vec.push_back(mesh);
        tetmesh = makeTetMesh(vec, "test");
    }
    virtual void TearDown() {}

    SurfaceMesh* mesh;
    std::unique_ptr<TetMesh> tetmesh;
};

TEST_F(DecimateTest, EdgeCollapseOp){
    auto s = (*mesh).get_simplex_up({0,1});
    //edgeCollapse<TetMesh, Callback>(tetmesh, s, 0, Callback<TetMesh>());
    EXPECT_EQ(0, 0);
}
