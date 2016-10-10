#pragma once

#include "util.h"
#include "SimplicialComplex.h"
#include "Orientable.h"
#include "Vertex.h"
#include <string>
#include <iostream>

struct FaceProperties
{
    int  m;   /**< @brief Marker */
    bool sel; /**< @brief selection flag */
};

struct Face : Orientable, FaceProperties
{
    Face() {}
    Face(Orientable orient, FaceProperties prop)
        : Orientable(orient), FaceProperties(prop)
    {}
};

struct Global
{
    bool  closed;                /**< @brief is the surface mesh closed or not */
    int   _marker;               /**< @brief doman marker, to be used when tetrahedralizing */
    float volume_constraint;     /**< @brief volume constraint of the tetrahedralized domain */
    bool  use_volume_constraint; /**< @brief flag that determines if the volume constraint is used */

    float min[3];                /**< @brief minimal coordinate of nodes */
    float max[3];                /**< @brief maximal coordinate of nodes */

    float avglen;                /**< @brief average edge length */

    bool hole;                   /**< @brief flag that determines if the mesh is a hole or not */
};

struct complex_traits
{
    using KeyType = int;
    using NodeTypes = util::type_holder<Global,Vertex,void,Face>;
    using EdgeTypes = util::type_holder<Orientable,Orientable,Orientable>;
};

using SurfaceMesh_ASC = simplicial_complex<complex_traits>;
using ASC = simplicial_complex<complex_traits>;

class SurfaceMesh : public ASC {
public:
    SurfaceMesh() : ASC(){};
    ~SurfaceMesh(){};
    void print() const{
        for(auto x: this->get_level<1>()) {
            std::cout << x << ", ";
        }
        std::cout << std::endl;
    } 
};

