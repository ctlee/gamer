#define _USE_MATH_DEFINES
#include <cmath>
#include "vertex.h"

Vertex operator+(const Vertex& A, const Vector& B){
    Vertex rval(A);
    rval += B;
    return rval;
}

Vertex operator-(const Vertex& A, const Vector& B){
    Vertex rval(A);
    rval -= B;
    return rval; 
}

Vector operator-(const Vertex& A, const Vertex& B){
    return A.position - B.position;
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

double magnitude(const Vector& A){
    return  std::sqrt(A|A);
}

double distance(const Vertex& A, const Vertex& B){
    return magnitude(A - B);
}

double angle(const Vertex& A, const Vertex& B, const Vertex& C){
    Vector AB = A-B;
    Vector AC = A-C;
    Vector BC = B-C;

    double AB2 = AB|AB;
    double AC2 = AC|AC;
    double BC2 = BC|BC;

    if( AB2 == 0 || BC2 == 0){
        std::cerr << "Some length == 0, can't compute angle. Exiting..." << std::endl;
        exit(1);
    }       
    
    return std::acos(0.5*(AB2+BC2-AC2)/std::sqrt(AB2*BC2)) * 180 / M_PI;
}