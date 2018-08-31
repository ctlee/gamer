/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2018
 * by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
 *    and Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * ***************************************************************************
 */


#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>

#include "Vertex.h"

Vertex operator+(const Vertex& A, const Vector& B){
    Vertex rval(A);
    rval += B;
    return rval;
}

Vector operator+(const Vertex& A, const Vertex& B){
    return A.position + B.position;
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

