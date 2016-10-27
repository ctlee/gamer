#include <iostream>
#include "WrapSurfaceMesh.h"
#include "SurfaceMesh.h"
#include "Vertex.h"
#include "Util.h"


using WSM = WrappedSurfaceMesh;

WSM::WrappedSurfaceMesh(){
    SurfaceMesh* mesh =  new SurfaceMesh();
    this->ptr_mesh = (void*) mesh;
}

WSM::~WrappedSurfaceMesh(){
    SurfaceMesh* mesh = (SurfaceMesh*) this->ptr_mesh; 
    delete mesh;
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
    std::pair<SurfaceMesh*, bool> result = readOFF(filename);
     if(result.second == false){
        std::cout << "Something bad happened...";
        exit(1);
    }
    WrappedSurfaceMesh* ptr  = new WrappedSurfaceMesh(result.first);
    return *ptr;
}

