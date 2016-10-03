#pragma once

class  WrappedSimplicialComplex{
public:
    WrappedSimplicialComplex();
    ~WrappedSimplicialComplex();
    void insertVertex(size_t n, double x, double y, double z);
    void insertFace(size_t a, size_t b, size_t c);
    int numVertices() const;
    int numFaces() const;
    void removeVertex(size_t n);
    void removeFace(size_t a, size_t b, size_t c);
    std::string as_string() const;
    void print();
private:
    void* ptr_asc;
};
