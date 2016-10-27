#pragma once
#include <string>

class  WrappedSurfaceMesh{
public:
    WrappedSurfaceMesh();
    WrappedSurfaceMesh(void* ptr) : ptr_mesh(ptr) {};
    ~WrappedSurfaceMesh();
    void insertVertex(int n, double x, double y, double z);
    void insertFace(int a, int b, int c);
    void removeVertex(int n);
    void removeFace(int a, int b, int c);
    void writeOFF(const std::string filename);
    std::string as_string() const;
private:
    void* ptr_mesh;
};

WrappedSurfaceMesh& loadOFF(const std::string& filename);
