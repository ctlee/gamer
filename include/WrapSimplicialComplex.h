
#pragma once

class  WrappedSimplicialComplex{
public:
    WrappedSimplicialComplex();
    ~WrappedSimplicialComplex();
    template<typename KeyType> 
    void insertVertex(KeyType n, double x, double y, double z);
    template<typename KeyType>
    void insertFace(KeyType a, KeyType b, KeyType c);
    int numVertices();
    int numFaces();
private:
    void* ptr_asc;
};
