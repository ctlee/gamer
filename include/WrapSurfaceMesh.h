#pragma once
#include <string>

class  WrappedSurfaceMesh{
public:
    WrappedSurfaceMesh();
    WrappedSurfaceMesh(void* ptr) : ptr_mesh(ptr) {};
    ~WrappedSurfaceMesh();
    void writeOFF(const std::string filename);
    std::string as_string() const;
private:
    void* ptr_mesh;
};

WrappedSurfaceMesh& loadOFF(const std::string& filename);
