/**
 * This is a scratch space for testing new functions and such on the fly.
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "stringutil.h"
#include "SurfaceMeshOld.h"
#include "SurfaceMesh.h"
#include "TetMesh.h"
#include "Vertex.h"
#include "tensor.h"

#include <libraries/casc/include/CASCFunctions.h>
#include <libraries/casc/include/SimplexSet.h>
#include <libraries/casc/include/SimplexMap.h>
#include <libraries/casc/include/decimate.h>
#include <libraries/casc/include/typetraits.h>

int  main(int argc, char *argv[])
{
    // if(argc != 2)
    // {
    //     std::cerr << "Wrong arguments passed" << std::endl;
    //     return -1;
    // }
    // auto mesh = readOFF(argv[1]);
    // if(mesh == nullptr){
    //     std::cout << "Something bad happened..." << std::endl;
    //     exit(1);
    // }
    // std::cout << "Done reading" << std::endl;

    auto mesh = cube(2);

    compute_orientation(*mesh);
    if(getVolume(*mesh) < 0){
        for(auto &face : mesh->get_level<3>())
            face.orientation *= -1;
    }

    for(auto &face : mesh->get_level<3>())
        face.marker = 23;

    // coarseIT(*mesh, 0.016, 1, 0);
    // coarseIT(*mesh, 2.5, 0, 10);

    // smoothMesh(*mesh, 15, 165, 5, true);
    // writeOFF("tmp.off", *mesh);

    auto &global = *mesh->get_simplex_up();
    global.closed = true;
    global.ishole = true;
    global.marker = -1;

    std::cout << "General protein statistics: " << std::endl;
    std::cout << "# Vertices: " << mesh->size<1>() << std::endl;
    std::cout << "# Edges: " << mesh->size<2>() << std::endl;
    std::cout << "# Faces: " << mesh->size<3>() << std::endl;

    auto box = cube(2);    

    for(auto &face : box->get_level<3>())
        face.marker = 50;

    auto &global2 = *box->get_simplex_up();
    global2.closed = true;
    global2.ishole = false;
    global2.marker = 5;

    Vector center;
    double radius;
    std::tie(center, radius) = getCenterRadius(*mesh);
    scale(*box, radius*2);

    std::cout << "General box statistics: " << std::endl;
    std::cout << "# Vertices: " << box->size<1>() << std::endl;
    std::cout << "# Edges: " << box->size<2>() << std::endl;
    std::cout << "# Faces: " << box->size<3>() << std::endl;

    std::vector<std::unique_ptr<SurfaceMesh>> meshes;

    std::cout << "Volume of mesh: " << getVolume(*mesh) << std::endl;
    std::cout << "Volume of box: " << getVolume(*box) << std::endl;

    writeOFF("boxa.off",*mesh);
    writeOFF("boxb.off",*box);

    meshes.push_back(std::move(mesh));
    meshes.push_back(std::move(box));

    std::cout << "\n\n\nBEGINNING TETRAHEDRALIZATION..." << std::endl;
    auto tetmesh = makeTetMesh(meshes, "pq1.4zYAQ");

    int cnt = 0;
    int cnt2 = 0;
    for (auto face : tetmesh->get_level<3>()){
        if (face.marker == 23)
            cnt++;
        else if (face.marker == 50)
            cnt2++;
    }
    std::cout << "There are " << cnt << " faces with marker 23." << std::endl;
    std::cout << "There are " << cnt2 << " faces with marker 50." << std::endl;


    SurfaceMesh protein, boxE, jnk;

    for(auto face : tetmesh->get_level_id<3>()){
        auto fdata = *face;
        if (fdata.marker == 23){
            auto name = tetmesh->get_name(face);
            protein.insert(name);
            auto vertices = tetmesh->down(std::move(tetmesh->down(face)));

            for (auto vertex : vertices){
                auto &data  = *(protein.get_simplex_up(tetmesh->get_name(vertex)));
                data = *vertex;              
            }

        } else if (fdata.marker == 50){
            auto name = tetmesh->get_name(face);
            boxE.insert(name);
            auto vertices = tetmesh->down(std::move(tetmesh->down(face)));

            for (auto vertex : vertices){
                auto &data = *(boxE.get_simplex_up(tetmesh->get_name(vertex)));
                data = *vertex;              
            }
        } else{
            // std::cout << fdata.marker << std::endl;
            // auto tets = tetmesh->up(face);
            // for (auto tet : tets){
            //     std::cout << (*tet).marker << std::endl;
            // }
        }
    }

    compute_orientation(protein);
    compute_orientation(boxE);

    // writeOFF("proteinExtracted.off", protein);
    // writeOFF("boxExtracted.off", boxE);
    // writeOFF("jnk.off", jnk);

    // writeOFF("test.off", *mesh);
    // auto oldmesh = SurfaceMeshOld::readOFF(argv[1]);
    // //oldmesh->coarse(0.016,1,0,-1);
    // oldmesh->coarse(2.5,0,10,-1);
    // 
    std::cout << "EOF" << std::endl;
}
