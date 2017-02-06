#pragma once

#include <iostream>
#include <sstream>
#include "tensor.h"


using Vector = tensor<double,3,1>;

/**
 * @brief      Vertex class
 */
struct Vertex
{
    Vector position;    /**< @brief a 3 tensor for x, y, z */
    int marker = 0;                 /**< @brief Boundary marking ID */
    bool selected = false;          /**< @brief Selection flag */
    
    /**
     * @brief      Default constructor with x,y,z = 0
     */
    Vertex(): Vertex(0,0,0) {}

    /**
     * @brief      Constructor
     *
     * @param[in]  x     x-position of the vertex
     * @param[in]  y     y-position of the vertex
     * @param[in]  z     z-position of the vertex
     */
    Vertex(double x, double y, double z): Vertex(x,y,z, 0, false){}

    /**
     * @brief      Constructor
     *
     * @param[in]  x     x-position of the vertex
     * @param[in]  y     y-position of the vertex
     * @param[in]  z     z-position of the vertex
     * @param[in]  m     marker ID
     * @param[in]  sel   selection flag
     */
    Vertex(double x, double y, double z, size_t m, bool sel){
        position[0] = x;
        position[1] = y;
        position[2] = z;
        marker = m;
        selected = sel;
    }
    /**
     * @brief      Copy Constructor
     *
     * @param[in]  x     Vertex to copy
     */
    Vertex(const Vertex& x) : position(x.position), marker(x.marker), selected(x.selected){}

    /**
     * @brief      Move Constructor
     *
     * @param[in]  x     Vertex to move
     */
    Vertex(const Vertex&& x) : position(std::move(x.position)), marker(std::move(x.marker)), selected(std::move(x.selected)) {}

    operator Vector() const {
        return position;
    }

    /**
     * @brief      Operator<< overload
     *
     * @param      output  stream to print to
     * @param[in]  v       Vertex to print
     *
     * @return     the stream
     */
    friend std::ostream& operator<<(std::ostream& output, const Vertex& v){
        output  << "Vertex(x:" << v[0]
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
    std::string to_string() const{
        std::ostringstream output;
        output  << "Vertex(x:" << position[0]
                << ",y:" << position[1] 
                << ",z:" << position[2]
                << ";m:" << marker
                << ";sel:" << selected
                << ")";
        return output.str();
    }

    /**
     * @brief      Const operator[] overload allows easy access to x, y, z using
     *             intuitive syntax
     *
     * @param[in]  index  Index to access
     *
     * @return     return the value of the tensor at index
     */
    const double& operator[](std::size_t index) const{
        return position[index];
    }
    
    /**
     * @brief      Operator[] overload allows easy access to x, y, z using intuitive syntax
     *
     * @param[in]  index  Index to access
     *
     * @return     the value of the tensor at index
     */
    double& operator[](std::size_t index){
        return position[index];
    }
    
    /**
     * @brief      Assignment operator overload
     *
     * @param[in]  v     vertex to assign
     */
    void operator=(const Vertex& v){
        position = v.position;
        marker = v.marker;
        selected = v.selected;
    }
   
    /**
     * @brief      Equivalence operator
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     True if all values are equal
     */
    bool operator==(const Vertex& rhs) const{
        Vertex temp(rhs); 
        if(position != temp.position) return false;
        if(marker != temp.marker) return false;
        if(selected != temp.selected) return false;
        return true;
    }

    /**
     * @brief      Inequivalence operator
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     True if not equal
     */
    bool operator!=(const Vertex& rhs) const{
        return !(*this == rhs);
    }

    /**
     * @brief      Add a vector to the vertex
     *
     * @param[in]  rhs   The right hand side
     *
     * @return     Vertex with sum of positions
     */
    Vertex& operator+=(const Vector& rhs){
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
    Vertex& operator-=(const Vector& rhs){
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
    Vertex& operator*=(double x){
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
    Vertex& operator/=(double x){
        position /=x;
        return *this;
    }
};

Vertex operator+(const Vertex& A, const Vector& B);
Vertex operator-(const Vertex& A, const Vector& B);
Vector operator-(const Vertex& A, const Vertex& B);
Vertex operator*(double x, const Vertex& A);
Vertex operator/(const Vertex& A, double x);
double distance(const Vertex& A, const Vertex& B);
double angle(const Vertex& A, const Vertex& B, const Vertex& C);
double angle(const Vector& AB, const Vector& CB);
double magnitude(const Vector& A);
