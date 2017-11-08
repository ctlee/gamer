#define _USE_MATH_DEFINES
#include <cmath>
#include "Vertex.h"

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

/**
 * @brief     Compute the angle betweeen three vertices. 
 *
 * @param[in]  A     Vertex A
 * @param[in]  B     Vertex B is in the middle
 * @param[in]  C     Vertex C
 *
 * @return     The angle in degrees
 */
double angle(const Vertex& A, const Vertex& B, const Vertex& C){
    Vector AB = A-B;
    Vector CB = C-B;
    return angle(AB,CB);
}

double angle(const Vector& AB, const Vector& CB){
    auto ab = AB;
    auto cb = CB;
    double lenAB = magnitude(ab); 
    double lenCB = magnitude(cb);
    if (lenAB == 0 || lenCB == 0){
        std::cerr << "Some length == 0, can't compute angle." << std::endl;
        return -1;
    } 
    ab /= lenAB;
    cb /= lenCB;
    return std::acos(ab|cb)*180/M_PI;
}