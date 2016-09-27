
#pragma once

class  WrappedSimplicialComplex{
public:
    WrappedSimplicialComplex();
    ~WrappedSimplicialComplex();
    void insertVertex(size_t n, double x, double y, double z);
    void insertFace(size_t a, size_t b, size_t c);
    int numVertices();
    int numFaces();
private:
    void* ptr_asc;
};
