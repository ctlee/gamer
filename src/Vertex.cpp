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

#include "gamer/Vertex.h"

/// Namespace for all things gamer
namespace gamer
{

Vertex operator+(const Vertex &A, const Vector &B)
{
    Vertex rval(A);
    rval += B;
    return rval;
}

Vertex operator-(const Vertex &A, const Vector &B)
{
    Vertex rval(A);
    rval -= B;
    return rval;
}

Vector operator-(const Vertex &A, const Vertex &B)
{
    return A.position - B.position;
}

Vertex operator*(REAL x, const Vertex &A)
{
    Vertex rval(A);
    rval *= x;
    return rval;
}

Vertex operator*(const Vertex &A, REAL x)
{
    return x*A;
}

Vertex operator/(const Vertex &A, REAL x)
{
    Vertex rval(A);
    rval /= x;
    return rval;
}

REAL distance(const Vertex &A, const Vertex &B)
{
    return length(A - B);
}

REAL angle(const Vertex &A, const Vertex &B, const Vertex &C)
{
    Vector AB(A-B);
    Vector CB(C-B);
    return angle(AB, CB);
}

REAL angle(const Vector &AB, const Vector &CB)
{

    return angleRad(AB, CB)*180/M_PI;
}

REAL angleRad(const Vertex &A, const Vertex &B, const Vertex &C)
{
    Vector AB(A-B);
    Vector CB(C-B);
    return angleRad(AB, CB);
}

REAL angleRad(const Vector &AB, const Vector &CB)
{
    Vector ab(AB);
    Vector cb(CB);
    REAL   lenAB = length(ab);
    REAL   lenCB = length(cb);
    if (lenAB == 0 || lenCB == 0)
    {
        throw std::runtime_error("Some length == 0, cannot compute angle.");
    }
    ab /= lenAB;
    cb /= lenCB;
    return std::acos(ab|cb);
}

REAL angleRadTan(const Vertex &A, const Vertex &B, const Vertex &C){
    Vector AB(A-B);
    Vector CB(C-B);
    return angleRadTan(AB, CB);
}

REAL angleRadTan(const Vector &v1, const Vector &v2){
    return std::atan2(length(cross(v1,v2)), dot(v1,v2));
}

REAL signed_angle(const Vector& v1, const Vector& v2, const Vector& reference)
{
    return std::atan2(length(dot(cross(v1, v2), reference)), dot(v1, v2));
}

} // end namespace gamer
