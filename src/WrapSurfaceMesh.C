#include <iostream>
#include <memory>
#include "WrapSurfaceMesh.h"
#include "SurfaceMesh.h"
#include "Vertex.h"
#include "util.h"


using WSM = WrappedSurfaceMesh;

WSM::WrappedSurfaceMesh(){
    SurfaceMesh* mesh =  new SurfaceMesh();
    this->ptr_mesh = (void*) mesh;
}

WSM::~WrappedSurfaceMesh(){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    delete mesh;
}

void WSM::insertVertex(int n, double x, double y, double z){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    Vertex v = Vertex(x,y,z);
    mesh->insert({n},v);
}

void WSM::insertFace(int a, int b, int c){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    mesh->insert({a,b,c});
}

void WSM::removeVertex(int n){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    mesh->remove({n});
}

void WSM::removeFace(int a, int b, int c){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    mesh->remove({a,b,c});
}

void WSM::writeOFF(const std::string filename){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    ::writeOFF(filename, *mesh);
}

std::string WSM::as_string() const{
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh;
    std::string out = "";
    for(auto& x : mesh->get_level<1>()){
        out += x.to_string() + ", ";
    }
    return out;
}

WrappedSurfaceMesh& loadOFF(const std::string& filename){
    // std::pair<SurfaceMesh*, bool> result = readOFF(filename);
    //  if(result.second == false){
    //     std::cout << "Something bad happened...";
    //     exit(1);
    // }

    auto mesh = readOFF(filename);
    if(mesh == nullptr){
        std::cout << "Something bad happened...";
        exit(1);
    }

    WrappedSurfaceMesh* ptr  = new WrappedSurfaceMesh(mesh.release());
    return *ptr;
}

