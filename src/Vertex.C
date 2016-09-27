#include "vertex.h"

Vertex operator+(const Vertex& A, const Vertex& B){
    Vertex rval(A);
    rval += B;
    return rval;
}

Vertex operator-(const Vertex& A, const Vertex& B){
    Vertex rval(A);
    rval -= B;
    return rval; 
}

Vertex operator*(double x, const Vertex& A){
    Vertex rval(A);
    rval *= x;
    return rval;
}

Vertex operator/(const Vertex& A, double x){
    Vertex rval(A);
    rval /= x;
    return rval;
}
