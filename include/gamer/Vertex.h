/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2019
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

/**
 * @file  Vertex.h
 * @brief Vertex class definitions
 */

#pragma once

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gamer/gamer.h"
#include "gamer/tensor.h"

/// Namespace for all things gamer
namespace gamer
{
/**
 * @brief      Vertex struct represents a general vertex
 */
struct Vertex
{
    Vector position;            /**< @brief a 3 tensor for x, y, z */
    int    marker   = 0;        /**< @brief Boundary marking ID */
    bool   selected = false;    /**< @brief Selection flag */

    /**
     * @brief      Default constructor with x,y,z = 0
     */
    Vertex() : Vertex(0, 0, 0) {}

    /**
     * @brief      Constructor with initialized position
     *
     * @param[in]  x     x-position of the vertex
     * @param[in]  y     y-position of the vertex
     * @param[in]  z     z-position of the vertex
     */
    Vertex(REAL x, REAL y, REAL z) : Vertex(x, y, z, -1, false){}

    /**
     * @brief      Constructor with initialized position, marker, and selection
     *
     * @param[in]  x     x-position of the vertex
     * @param[in]  y     y-position of the vertex
     * @param[in]  z     z-position of the vertex
     * @param[in]  m     marker ID
     * @param[in]  sel   selection flag
     */
    Vertex(REAL x, REAL y, REAL z, int m, bool sel)
    {
        position[0] = x;
        position[1] = y;
        position[2] = z;
        marker   = m;
        selected = sel;
    }

    /**
     * @brief      Constructor seeded from a Vector
     *
     * @param[in]  v     Vector position
     */
    Vertex(Vector &v) : position(v), marker(-1), selected(false){}

    /**
     * @brief      Move construct a vertex from vector
     *
     * @param      v     Vector position
     */
    Vertex(Vector &&v) : position(std::move(v)), marker(-1), selected(false) {}

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  x     Vertex to copy
     */
    Vertex(const Vertex &x) : position(x.position), marker(x.marker), selected(x.selected){}

    /**
     * @brief      Move Constructor
     *
     * @param[in]  x     Vertex to move
     */
    Vertex(const Vertex &&x) : position(std::move(x.position)), marker(std::move(x.marker)), selected(std::move(x.selected)) {}

    /**
     * @brief      Implicit cast of Vertex to Vector type
     */
    operator Vector() const {
        return position;
    }

    /**
     * @brief      Print operator overload
     *
     * @param      output  stream to print to
     * @param[in]  v       Vertex to print
     *
     * @return     the stream
     */
    friend std::ostream &operator<<(std::ostream &output, const Vertex &v)
    {
        output << "Vertex(x:" << v[0]
               << ",y:" << v[1]
               << ",z:" << v[2]
               << ";m:" << v.marker
               << ";sel:" << v.selected
               << ")";
        return output;
    }

    /**
     * @brief      Returns a string representation of the object.
     *
     * @return     String representation of the object.
     */
    std::string to_string() const
    {
        std::ostringstream output;
        output << *this;
        return output.str();
    }

    /**
     * @brief      Const operator[] overload allows easy access to x, y, z using
     *             intuitive syntax
     *
     * @param[in]  index  Index to access
     *
     * @return     Reference to the value at the index
     */
    const REAL &operator[](std::size_t index) const
    {
        return position[index];
    }

    /**
     * @brief      Operator[] overload allows easy access to x, y, z using
     *intuitive syntax
     *
     * @param[in]  index  Index to access
     *
     * @return     Reference to the value at the index
     */
    REAL &operator[](std::size_t index)
    {
        return position[index];
    }

    /**
     * @brief      Assignment operator overload
     *
     * @param[in]  v     vertex to assign
     */
    void operator=(const Vertex &v)
    {
        position = v.position;
        marker   = v.marker;
        selected = v.selected;
    }

    /**
     * @brief      Equivalence operator
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     True if all values are equal
     */
    bool operator==(const Vertex &rhs) const
    {
        Vertex temp(rhs);
        if (position != temp.position) return false;
        if (marker != temp.marker) return false;
        if (selected != temp.selected) return false;
        return true;
    }

    /**
     * @brief      Inequivalence operator
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     True if not equal
     */
    bool operator!=(const Vertex &rhs) const
    {
        return !(*this == rhs);
    }

    /**
     * @brief      Add a vector to the vertex
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     Vertex with sum of positions
     */
    Vertex &operator+=(const Vector &rhs)
    {
        // retains the marker of the lhs
        position += rhs;
        return *this;
    }

    /**
     * @brief      Subtracts a vector from a vertex
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     Vertex with difference of positions
     */
    Vertex &operator-=(const Vector &rhs)
    {
        // retains the marker of the lhs
        position -= rhs;
        return *this;
    }

    /**
     * @brief      Multiply the Vertex by a scalar
     *
     * @param[in]  x     The scalar to multiply by
     *
     * @return     Post multiplied vertex
     */
    Vertex &operator*=(const REAL x)
    {
        position *= x;
        return *this;
    }

    /**
     * @brief     Divide the Vertex by a scalar
     *
     * @param[in]  x     Scalar to divide by
     *
     * @return     Post divided vertex
     */
    Vertex &operator/=(const REAL x)
    {
        position /= x;
        return *this;
    }
};

/**
 * @brief      Translate a Vertex by vector
 *
 * @param[in]  A     Vertex base
 * @param[in]  B     Vector translation
 *
 * @return     Resulting Vertex with new position
 */
Vertex operator+(const Vertex &A, const Vector &B);

/**
 * @brief      Translate a Vertex by difference
 *
 * @param[in]  A     Vertex base
 * @param[in]  B     Vector translation
 *
 * @return     Resulting Vertex with new position
 */
Vertex operator-(const Vertex &A, const Vector &B);

/**
 * @brief      Compute vector between two Vertices
 *
 * @param[in]  A     First Vertex
 * @param[in]  B     Second Vertex
 *
 * @return     Vector from Vertex A to Vertex B
 */
Vector operator-(const Vertex &A, const Vertex &B);

/**
 * @brief      Scale Vertex position
 *
 * @param[in]  x     Scalar multiplicand
 * @param[in]  A     Vertex to scale
 *
 * @return     Vertex with scaled position
 */
Vertex operator*(REAL x, const Vertex &A);

/**
 * @brief      Scale Vertex position
 *
 * @param[in]  A     Vertex to scale
 * @param[in]  x     Scalar multiplicand
 *
 * @return     Vertex with scaled position
 */
Vertex operator*(const Vertex &A, REAL x);

/**
 * @brief      Scale Vertex position by scalar division
 *
 * @param[in]  A     Vertex to scale
 * @param[in]  x     Scalar denominator
 *
 * @return     Vertex with scaled position
 */
Vertex operator/(const Vertex &A, REAL x);

/**
 * @brief      Get the distance between two Vertices
 *
 * @param[in]  A     The first Vertex
 * @param[in]  B     The other Vertex
 *
 * @return     Scalar distance between vertices
 */
REAL distance(const Vertex &A, const Vertex &B);

/**
 * @brief     Compute the angle between three vertices.
 *
 * @param[in]  A     Vertex A
 * @param[in]  B     Vertex B is in the middle
 * @param[in]  C     Vertex C
 *
 * @return     The angle in radians
 */
REAL angle(const Vertex &A, const Vertex &B, const Vertex &C);

/**
 * @brief      Compute angle between two vectors
 *
 * @param[in]  AB    First Vector
 * @param[in]  CB    Second Vector
 *
 * @return     The angle in radians
 */
REAL angle(const Vector &AB, const Vector &CB);

// REAL signed_angle(const Vertex &A, const Vertex &B, const Vertex &C);
REAL signed_angle(const Vector &v1, const Vector &v2, const Vector &reference);

/**
 * @brief     Compute the angle between three vertices.
 *
 * @param[in]  A     Vertex A
 * @param[in]  B     Vertex B is in the middle
 * @param[in]  C     Vertex C
 *
 * @return     The angle in degrees
 */
REAL angleDeg(const Vertex &A, const Vertex &B, const Vertex &C);

/**
 * @brief      Compute angle between two vectors
 *
 * @param[in]  AB    First Vector
 * @param[in]  CB    Second Vector
 *
 * @return     The angle in degrees
 */
REAL angleDeg(const Vector &AB, const Vector &CB);


/**
 * @brief      Get the length of a vector
 *
 * @param[in]  A     Vector of interest
 *
 * @return     length of the vector
 */
inline REAL length(const Vector &A)
{
    return std::sqrt(A|A);
}

/**
 * @brief      Normalize a vector
 *
 * @param      A     Vector of interest
 */
inline void normalize(Vector &A)
{
    REAL mag = length(A);
    if (mag == 0)
        throw std::runtime_error("Cannot normalize a vector with length of 0.");
    A /= mag;
}

} //end namespace gamer
