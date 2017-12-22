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

#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>

#include <libraries/casc/casc>
//#include <libraries/casc/include/typetraits.h>

#include "TetMesh.h"

std::unique_ptr<TetMesh> makeTetMesh(
        const std::vector<std::unique_ptr<SurfaceMesh>> &surfmeshes, 
        std::string tetgen_params){

    // Create new tetmesh object
    std::unique_ptr<TetMesh> tetmesh(new TetMesh);

    size_t nVertices = 0, nFaces = 0, nRegions = 0, nHoles = 0;
    int i = 0;

    for (auto &surfmesh : surfmeshes){
        size_t nverts = surfmesh->template size<1>();
        size_t nfaces = surfmesh->template size<3>();

        if(nverts == 0 || nfaces == 0){
            std::cerr << "The " << i << "th surface mesh is empty."
                      << " Abandoning tetrahedralization." << std::endl;
            tetmesh.reset();
            return tetmesh;
        }

        auto metadata = *surfmesh->get_simplex_up();

        if (!metadata.closed){
            std::cerr << "The " << i << "th surface mesh is not closed."
                    << " Abandoning tetrahedralization." << std::endl;
            tetmesh.reset();
            return tetmesh;
        }
        metadata.ishole ? ++nHoles : ++nRegions;

        nVertices += nverts;
        nFaces += nfaces;
        ++i;
    }

    if (nRegions < 1){
        std::cerr << "Error(makeTetMesh): Expected at least one non-hole SurfaceMesh." << std::endl;
        tetmesh.reset();
        return tetmesh;
    }

    std::cout << "Number of vertices: " << nVertices << std::endl;
    std::cout << "Number of Faces: " << nFaces << std::endl;
    std::cout << "Number of Regions: " << nRegions << std::endl;
    std::cout << "Number of Holes: " << nHoles << std::endl;

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
    nFaces = nRegions = nHoles = 0;
    typename TetMesh::KeyType cnt = 0;

    for (auto &surfmesh : surfmeshes){
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
        for (const auto faceID : surfmesh->template get_level_id<3>()){
            auto w = surfmesh->get_name(faceID);

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

            in.facetmarkerlist[nFaces] = (*faceID).marker == -1? 0 : (*faceID).marker;
            ++nFaces;
        }

        auto metadata = *surfmesh->get_simplex_up();

        // Pick a point inside the region
        auto faceID = *surfmesh->template get_level_id<3>().begin();
        Vector normal = getNormal(*surfmesh, faceID);
        normal /= std::sqrt(normal|normal);

        auto fname = surfmesh->get_name(faceID);
        Vertex a = *surfmesh->get_simplex_up({fname[0]});
        Vertex b = *surfmesh->get_simplex_up({fname[1]});
        Vertex c = *surfmesh->get_simplex_up({fname[2]});

        Vector d = a-b;
        double weight = std::sqrt(d|d);

        // flip normal and scale by weight
        normal *= weight;
        Vector regionPoint = (a+b+c)/3 - normal;

        if(metadata.ishole){
            auto idx            = nHoles*3;
            in.holelist[idx]    = regionPoint[0];
            in.holelist[idx+1]  = regionPoint[1];
            in.holelist[idx+2]  = regionPoint[2];
            ++nHoles;
        }
        else
        {
            auto idx                = nRegions*5;
            in.regionlist[idx]      = regionPoint[0];
            in.regionlist[idx+1]    = regionPoint[1];
            in.regionlist[idx+2]    = regionPoint[2];
            in.regionlist[idx+3]    = metadata.marker;
            // std::cout << "Region marker: " << metadata.marker << std::endl;

            if(metadata.useVolumeConstraint){
                in.regionlist[idx+4]    = metadata.volumeConstraint; 
            }
            else{
                in.regionlist[idx+4]    = -1; 
            }
            ++nRegions; 
        }
    } // endif for surfmesh :surfmeshes

    // Add oundary marker on each node
    // TODO: Why? aren't the markers set on the generated mesh? (from old notes)
    in.pointmarkerlist = new int[in.numberofpoints];
    for (int i = 0; i < in.numberofpoints; ++i){
        in.pointmarkerlist[i] = 1;
    }
   
    // Casting away const is an evil thing to do, however, tetgen has not yet
    // conformed... 
    auto plc = const_cast<char*>("plc");
    in.save_nodes(plc);
    in.save_poly(plc);

    // Call TetGen
    tetrahedralize(tetgen_params.c_str(), &in, &out, NULL);

    auto result = const_cast<char*>("result");
    out.save_nodes(result);
    out.save_elements(result);
    out.save_faces(result);

    return tetgenToTetMesh(out);
}

std::unique_ptr<TetMesh> tetgenToTetMesh(tetgenio &tetio){
    std::unique_ptr<TetMesh> mesh(new TetMesh);

    if (tetio.mesh_dim == 2){
        std::cerr << "ERROR(tetgenToTetMesh): Expected a tetrahedral mesh from tetgen." << std::endl;
        mesh.reset();
        return mesh;
    }

    auto &metadata = *mesh->get_simplex_up();

    // Check for higher order cells
    const bool higher_order = tetio.numberofcorners == 10;
   
    metadata.higher_order = higher_order;

    std::set<int> vertices;

    // std::cout << "Number of tetrahedron attributes: " << tetio.numberoftetrahedronattributes
    //     << std::endl;

    // Copy over tetrahedron data
    // std::cout << "Copying over tetrahedron data..." << std::endl;
    for (int i = 0; i < tetio.numberoftetrahedra; ++i){
        // Set marker
        int marker = 0;

        if(tetio.numberoftetrahedronattributes > 0) 
            marker = (int) tetio.tetrahedronattributelist[i * tetio.numberoftetrahedronattributes];

        // Get vertex id's
        int *ptr = &tetio.tetrahedronlist[i*tetio.numberofcorners];

        vertices.insert({ptr[0], ptr[1], ptr[2], ptr[3]});

        // TODO: (0) Do we need to set the orientation?
        // std::cout << casc::to_string(std::array<int,4>({ptr[0], ptr[1], ptr[2], ptr[3]})) << std::endl;
        mesh->insert<4>({ptr[0], ptr[1], ptr[2], ptr[3]}, 
                tetmesh::Cell(casc::Orientable{0}, tetmesh::CellProperties{marker, 0, 0}));
    }

    // std::cout << "Number of vertices: " << mesh->size<1>() << std::endl;
    // std::cout << "Number of edges: " << mesh->size<2>() << std::endl;
    // std::cout << "Number of faces: " << mesh->size<3>() << std::endl;
    // std::cout << "Number of tetrahedra: " << mesh->size<4>() << std::endl;

    // Copy over vertex data
    // std::cout << "Copying over vertex data..." << std::endl;
    for (auto i : vertices){
        double *ptr = &tetio.pointlist[i*3];
        auto vertex = mesh->get_simplex_up({i});
        if (vertex != nullptr){
            auto &vdata = *vertex;
            // std::cout << casc::to_string(std::array<double,3>({ptr[0],ptr[1],ptr[2]})) << std::endl;
            vdata = Vertex(ptr[0], ptr[1], ptr[2], tetio.pointmarkerlist[i], false);
        }
    }

    // Go over faces and copy over marker information
    // std::cout << "Copying over face data..." << std::endl;
    for (int i = 0; i < tetio.numberoftrifaces; ++i){
    	int *ptr = &tetio.trifacelist[i*3];
        auto face = mesh->get_simplex_up({ptr[0], ptr[1], ptr[2]});
        if (face != nullptr){
            auto &fdata = *face;
            fdata.marker = tetio.trifacemarkerlist[i];
        }
    }

    // Copy over edge markers
    // std::cout << "Copying over edge data..."  << std::endl;
    for (int i = 0; i < tetio.numberofedges; ++i){
        int *ptr = &tetio.edgelist[i*2];

        auto edgeID = mesh->get_simplex_up({ptr[0], ptr[1]});
        if(edgeID != nullptr){
            auto &edata = *edgeID;
            edata.marker = tetio.edgemarkerlist[i];

            if (higher_order){
                double *pos = &tetio.pointlist[tetio.o2edgelist[i]*3];
                edata.position = Vector({pos[0], pos[1], pos[2]});                
            }
        }
    }
    compute_orientation(*mesh);
    return mesh;
}

void writeVTK(const std::string& filename, const TetMesh &mesh){
    std::ofstream fout(filename);
    if(!fout.is_open())
    {
        std::cerr   << "File '" << filename 
                    << "' could not be writen to." << std::endl;
        exit(1); 
    }

    fout << "# vtk DataFile Version 2.0\n"
         << "Unstructured Grid\n"
         << "ASCII\n"  // BINARY
         << "DATASET UNSTRUCTURED_GRID\n";

    std::map<typename TetMesh::KeyType,typename TetMesh::KeyType> sigma;
    typename TetMesh::KeyType cnt = 0;

    // Output vertices
    fout << "POINTS " << mesh.size<1>() << " double" << std::endl;
    for (auto vertexID : mesh.get_level_id<1>()){
        sigma[mesh.get_name(vertexID)[0]] = cnt++;
        auto vertex = *vertexID;
        fout    << std::setprecision(17) << vertex[0] << " "
                << vertex[1] << " " 
                << vertex[2] << "\n";
    }
    fout << "\n";

    bool orientationError = false;

    fout << "CELLS " << mesh.size<4>() << " " << mesh.size<4>()*(4+1) << "\n";
    for (auto cellID : mesh.get_level_id<4>()){
        auto w = mesh.get_name(cellID);
        auto orientation = (*cellID).orientation;

        if (orientation == 1){
            fout << "4 " << std::setw(4) << sigma[w[3]] << " " 
                         << std::setw(4) << sigma[w[1]] << " " 
                         << std::setw(4) << sigma[w[2]] << " " 
                         << std::setw(4) << sigma[w[0]] << "\n";
        }
        else if(orientation == -1){
            fout << "4 " << std::setw(4) << sigma[w[0]] << " " 
                         << std::setw(4) << sigma[w[1]] << " " 
                         << std::setw(4) << sigma[w[2]] << " " 
                         << std::setw(4) << sigma[w[3]] << "\n";
        }
        else{
            orientationError = true;
            fout << "4 " << std::setw(4) << sigma[w[0]] << " " 
                         << std::setw(4) << sigma[w[1]] << " " 
                         << std::setw(4) << sigma[w[2]] << " " 
                         << std::setw(4) << sigma[w[3]] << "\n";
        }
    }
    fout << "\n";

    fout << "CELL_TYPES " << mesh.size<4>() << "\n";
    for (int i = 0; i < mesh.size<4>(); ++i){
        fout << "10\n";
    }
    fout << "\n";

    fout << "CELL_DATA " << mesh.size<4>() << "\n";
    fout << "SCALARS cell_scalars int 1\n";
    fout << "LOOKUP_TABLE default\n";
    // This should output in the same order...
    for (auto cell : mesh.get_level<4>()){
        fout << cell.marker << "\n";
    }
    fout << "\n";

    if(orientationError){
        std::cerr << "WARNING(writeVTK): The orientation of one or more faces "
                  << "is not defined. Did you run compute_orientation()?" 
                  << std::endl;
    }
    fout.close();
}

void writeOFF(const std::string& filename, const TetMesh &mesh){
    std::ofstream fout(filename);
    if(!fout.is_open())
    {
        std::cerr   << "File '" << filename 
                    << "' could not be writen to." << std::endl;
        exit(1); 
    }

    fout << "OFF\n";
    fout << mesh.size<1>() << " " 
         << mesh.size<4>() << " " 
         << mesh.size<2>() << "\n";

    std::map<typename TetMesh::KeyType,typename TetMesh::KeyType> sigma;
    typename TetMesh::KeyType cnt = 0;

    fout.precision(10); 
    for(const auto vertexID : mesh.get_level_id<1>()){
        sigma[mesh.get_name(vertexID)[0]] = cnt++;
        auto vertex = *vertexID;

        fout    << vertex[0] << " "
                << vertex[1] << " " 
                << vertex[2] << " " 
                << "\n";
    }

    bool orientationError = false;

    for (auto cellID : mesh.get_level_id<4>()){
        auto w = mesh.get_name(cellID);
        auto orientation = (*cellID).orientation;

        if (orientation == 1){
            fout << "4 " << std::setw(4) << sigma[w[3]] << " " 
                         << std::setw(4) << sigma[w[1]] << " " 
                         << std::setw(4) << sigma[w[2]] << " " 
                         << std::setw(4) << sigma[w[0]] << "\n";
        }
        else if(orientation == -1){
            fout << "4 " << std::setw(4) << sigma[w[0]] << " " 
                         << std::setw(4) << sigma[w[1]] << " " 
                         << std::setw(4) << sigma[w[2]] << " " 
                         << std::setw(4) << sigma[w[3]] << "\n";
        }
        else{
            orientationError = true;
            fout << "4 " << std::setw(4) << sigma[w[0]] << " " 
                         << std::setw(4) << sigma[w[1]] << " " 
                         << std::setw(4) << sigma[w[2]] << " " 
                         << std::setw(4) << sigma[w[3]] << "\n";
        }
    }

    if(orientationError){
        std::cerr << "WARNING(writeOFF): The orientation of one or more cells "
                  << "is not defined. Did you run compute_orientation()?" 
                  << std::endl;
    }
    fout.close();
}

