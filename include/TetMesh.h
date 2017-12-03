/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2017
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

#pragma once

#include <iostream>
#include <libraries/casc/include/SimplicialComplex.h>
#include <libraries/casc/include/CASCTraversals.h>
#include <libraries/casc/include/util.h>
#include <libraries/casc/include/Orientable.h>
#include <libraries/casc/include/stringutil.h>

#define TETLIBRARY

#include <libraries/tetgen/tetgen.h>

#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include "SurfaceMesh.h"
#include "Vertex.h"

namespace tetmesh_detail{
/**
 * @brief      Properties that Faces should have
 */
struct FaceProperties
{
    int  marker;   /**< @brief Marker */
    bool selected; /**< @brief Selection flag */
};

/**
 * @brief      Face object
 */
struct Face : FaceProperties
{
    /// Default constructor
    Face() {}
    /**
     * @brief      Constructor
     *
     * @param[in]  prop    Properties of a face
     */
    Face(FaceProperties prop)
        : FaceProperties(prop)
    {}
};

struct CellProperties
{
    int id;
    int grpID;
    int material;
};

struct Cell : casc::Orientable, CellProperties
{
    Cell() {}
    Cell(Orientable orient, CellProperties prop)
        : Orientable(orient), CellProperties(prop)
    {}
};


/**
 * @brief      Type for containing root metadata
 */
struct Global
{
    bool higher_order;
};

/**
 * @brief      A helper struct containing the traits/types in the simplicial
 *             complex
 */
struct complex_traits
{
    /// The index type
    using KeyType = int;
    /// The types of each node
    using NodeTypes = util::type_holder<Global, Vertex, void, Face, Cell>;
    /// The types of each edge
    using EdgeTypes = util::type_holder<casc::Orientable, casc::Orientable, casc::Orientable, casc::Orientable>;
};
}

using TetMesh = casc::simplicial_complex<tetmesh_detail::complex_traits>;


/// Forward class declaration
class tetgenio;


// TODO: Change this to use a more flexible vector of surface meshes...
template <size_t n>
std::unique_ptr<TetMesh> makeTetMesh(
        const std::array<std::unique_ptr<SurfaceMesh>, n> &surfmeshes, 
        char *tetgen_params){

    // Create new tetmesh object
    std::unique_ptr<TetMesh> tetmesh(new TetMesh);

    size_t nVertices = 0, nFaces = 0, nRegions = 0, nHoles = 0;

    for (int i = 0; i < n; ++i){
        auto surfmesh = surfmeshes[i];
        
        size_t nverts = surfmesh->template size<1>();
        size_t nfaces = surfmesh->template size<3>();

        if(nverts == 0 || nfaces == 0){
            std::cerr << "The " << i << "th surface mesh is empty."
                      << " Abandoning tetrahedralization." << std::endl;
            tetmesh.reset();
            return tetmesh;
        }

        auto metadata = *surfmesh.get_simplex_up();

        if (!metadata.closed){
            std::cerr << "The " << i << "th surface mesh is not closed."
                    << " Abandoning tetrahedralization." << std::endl;
            tetmesh.reset();
            return tetmesh;
        }
        metadata.ishole ? ++nRegions : ++nHoles;

        nVertices += nverts;
        nFaces += nfaces;
    }

    if (nRegions < 1){
        std::cerr << "Error(makeTetMesh): Expected at least one non-hole SurfaceMesh." << std::endl;
        tetmesh.reset();
        return tetmesh;
    }

    tetgenio in, out;

    tetgenio::facet *f;
    tetgenio::polygon *p;

    // Allocate memory for tetgenio object arrays
    in.numberofpoints   = nVertices;
    in.pointlist        = new REAL[nVertices * 3];

    in.numberoffacets   = nFaces;
    in.facetlist        = new tetgenio::facet[nFaces];
    in.facetmarkerlist  = new int[nFaces];

    in.numberofregions  = nRegions;
    in.regionlist       = new REAL[nRegions * 5];

    if (nHoles > 0){
        in.numberofholes    = nHoles;
        in.holelist         = new REAL[nHoles * 3];
    }

    // Reset counters
    nVertices = nFaces = nRegions = nHoles = 0;
    typename TetMesh::KeyType cnt = 0;

    for (auto surfmesh : surfmeshes){
        std::map<typename TetMesh::KeyType, typename TetMesh::KeyType> sigma;

        // Assign vertex information
        for(const auto vertexID : surfmesh->template get_level_id<1>())
        {
            sigma[surfmesh->get_name(vertexID)[0]] = cnt;

            auto vertex = *vertexID;
            auto idx = cnt*3;
            in.pointlist[idx] = vertex[0];
            in.pointlist[idx+1] = vertex[1];
            in.pointlist[idx+2] = vertex[2];

            ++cnt;
        }

        // Assign face information
        for (const auto faceID : surfmesh.template get_level_id<3>()){
            auto w = surfmesh.get_name(faceID);

            f                       = &in.facetlist[nFaces];
            f->holelist             = (REAL *)NULL;
            f->numberofholes        = 0;
            f->numberofpolygons     = 1;
            f->polygonlist          = new tetgenio::polygon[f->numberofpolygons];
            p                       = &f->polygonlist[0];
            p->numberofvertices     = 3;
            p->vertexlist           = new int[p->numberofvertices];
            p->vertexlist[0]        = sigma[w[0]];
            p->vertexlist[1]        = sigma[w[1]];
            p->vertexlist[2]        = sigma[w[2]];

            in.facetmarkerlist[nFaces] = *faceID.marker == -1? 0 : *faceID.marker;
            ++nFaces;
        }

        auto metadata = *surfmesh.get_simplex_up();

        // Pick a point inside the region
        auto faceID = *surfmesh.template get_level_id<3>().begin();
        Vector normal = getNormal(*surfmesh, faceID);
        normal /= std::sqrt(normal|normal);

        auto fname = surfmesh->get_name(faceID);
        Vertex a = *surfmesh->get_simplex_up({fname[0]});
        Vertex b = *surfmesh->get_simplex_up({fname[1]});
        Vertex c = *surfmesh->get_simplex_up({fname[2]});

        Vector d = a-b;
        double weight = std::sqrt(d|d); 
        normal *= weight;

        Vector midpoint = (a+b+c)/3 * weight;

        if(!metadata.ishole){
            auto idx                = nRegions*5;
            in.regionlist[idx]      = midpoint[0];
            in.regionlist[idx+1]    = midpoint[1];
            in.regionlist[idx+2]    = midpoint[2];

            in.regionlist[idx+3]    = metadata.marker;

            if(metadata.useVolumeConstraint){
                in.regionlist[idx+4]    = metadata.volumeConstraint; 
            }
            else{
                in.regionlist[idx+4]    = -1; 
            }
        }
        else{
            auto idx            = nHoles*3;
            in.holelist[idx]    = midpoint[0];
            in.holelist[idx+1]  = midpoint[1];
            in.holelist[idx+2]  = midpoint[2];
            ++nHoles;
        }
        ++nRegions;
    } // endif for surfmesh :surfmeshes

    // Add oundary marker on each node
    // TODO: Why? aren't the markers set on the generated mesh? (from old notes)
    in.pointmarkerlist = new int[in.numberofpoints];
    for (int i = 0; i < in.numberofpoints; ++i){
        in.pointmakerlist[j] = 1;
    }
   
    // Casting away const is an evil thing to do, however, tetgen has not yet
    // conformed... 
    auto plc = const_cast<char*>("plc");
    in.save_nodes(plc);
    in.save_poly(plc);

    // Call TetGen
    tetrahedralize(tetgen_params, &in, &out, NULL);

    auto result = const_cast<char*>("result");
    out.save_nodes(result);
    out.save_elements(result);
    out.save_faces(result);

    tetmesh.swap(tetgenToTetMesh(out));

    return tetmesh;
}


std::unique_ptr<TetMesh> tetgenToTetMesh(tetgenio &tetio);



void writeOFF(const std::string& filename, const TetMesh &mesh);
void writeMCSF(const std::string &filename, const TetMesh &mesh);
void writeDolfin(const std::string &filename, const TetMesh &mesh);

//void writeDiffPack
//void writeCARP



