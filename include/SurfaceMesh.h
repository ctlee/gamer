#pragma once

#include <iostream>
#include <string>
#include <utility>
#include "util.h"
#include "SimplicialComplex.h"
#include "Orientable.h"
#include "Vertex.h"

/**
 * @brief      Properties that Faces should have
 */
struct FaceProperties
{
    int  marker;   /**< @brief Marker */
    bool selected; /**< @brief selection flag */
};

/**
 * @brief      Face object
 */
struct Face : Orientable, FaceProperties
{
    Face() {}
    Face(Orientable orient, FaceProperties prop)
        : Orientable(orient), FaceProperties(prop)
    {}
};

/**
 * @brief      Type for containing root metadata
 */
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

/**
 * @brief      A helper struct containing the traits/types in the simplicial
 *             complex
 */
struct complex_traits
{
    using KeyType = int;                                                    /**< @brief the index type */
    using NodeTypes = util::type_holder<Global,Vertex,void,Face>;           /**< @brief the types of each Node */
    using EdgeTypes = util::type_holder<Orientable,Orientable,Orientable>;  /**< @brief the types of each Edge */
};

// This alias is for legacy purposes...
using SurfaceMesh_ASC = simplicial_complex<complex_traits>;
using ASC = simplicial_complex<complex_traits>; // Alias for the lazy
using SurfaceMesh = simplicial_complex<complex_traits>;

/**
 * @brief      Reads off.
 *
 * @param[in]  filename  The filename
 *
 * @return     { description_of_the_return_value }
 */
std::pair<SurfaceMesh*, bool> readOFF(const std::string& filename);

/**
 * @brief      Writes off.
 *
 * @param[in]  filename  The filename
 * @param[in]  mesh      The mesh
 */
void writeOFF(const std::string& filename, const SurfaceMesh& mesh);

/**
 * @brief      Convenience function to print the vertices
 *
 * @param[in]  mesh  The mesh
 */
void print_vertices(const SurfaceMesh& mesh);


/**
 * @brief      Convenience function to print the faces
 *
 * @param[in]  mesh  The mesh
 */
void print_faces(const SurfaceMesh& mesh);
void translate(SurfaceMesh& mesh, Vector v);
void translate(SurfaceMesh& mesh, double dx, double dy, double dz);
void scale(SurfaceMesh& mesh, Vector v);
void scale(SurfaceMesh& mesh, double sx, double sy, double sz);
void scale(SurfaceMesh& mesh, double s);
void generateHistogram(const SurfaceMesh& mesh);
bool smoothMesh(const SurfaceMesh& mesh, std::size_t minAngle, std::size_t maxAngle, std::size_t maxIter, bool preserveRidges);
void edgeFlip(SurfaceMesh& mesh, SurfaceMesh::NodeID<2> edgeID, bool preserveRidges);
int getValence(SurfaceMesh& mesh, SurfaceMesh::NodeID<1> nodeID);
