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



#include <cmath>
#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <iostream>
#include <vector>

#include "gamer/stringutil.h"
#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{
//https://en.wikipedia.org/wiki/Wavefront_.obj_file
//http://paulbourke.net/dataformats/obj/
std::unique_ptr<SurfaceMesh> readOBJ(const std::string &filename)
{
    // Instantiate mesh!
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::ifstream                fin(filename);
    if (!fin.is_open())
    {
        std::cerr << "Read Error: File '" << filename << "' could not be read." << std::endl;
        return mesh;
    }

    int                      i = 0; // index of vertices
    std::string              line;
    std::vector<std::string> arr;

    // while the file isn't empty
    while (fin.peek() != -1)
    {
        getline(fin, line);
        stringutil::trim(line);

        if (line.empty())
            continue;       // skip empty lines
        if (line[0] == '#')
            continue;       // skip comments

        if (line[0] == 'v') // Vertex
        {
            // if not space see what kind of data it is
            if (!std::isspace(line[1]))
            {
                // ignore normals "vn" and textures "vt"
                // also ignore potentially diabolical filenames
                continue;
            }
            else
            {
                // List of geometric vertices, with (x,y,z[,w]) coordinates, w
                // is optional and defaults to 1.0.
                arr = stringutil::split(line, {' '});
                SMVertex v = SMVertex(
                    std::stod(arr[1]),
                    std::stod(arr[2]),
                    std::stod(arr[3])
                    );
                // ignore possible w for now...
                mesh->insert<1>({++i}, v);
            }
        }

        // Faces can be a pain also arbitrary dimension
        // f v1 v2 v3 ....
        // f v1/vt1 v2/vt2 v3/vt3 ...
        // f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...
        // f v1//vn1 v2//vn2 v3//vn3 ...
        if (line[0] == 'f')
        {
            arr = stringutil::split(line, {' '});
            if (arr.size() != 4)
            {
                std::cerr << "Unsupported: Found face that is not a triangle!" << std::endl;
                mesh.reset();
                return mesh;
            }
            // again we're going to ignore textures and normals
            for (auto it = arr.begin(); it != arr.end(); ++it)
            {
                *it = stringutil::split(*it, {'/'})[0];
            }
            mesh->insert<3>({
                    std::stoi(arr[1]),
                    std::stoi(arr[2]),
                    std::stoi(arr[3])});
        }

        // everything else is ignored for now also
    }
    return mesh;
}

void writeOBJ(const std::string &filename, const SurfaceMesh &mesh)
{
    std::ofstream fout(filename);
    if (!fout.is_open())
    {
        std::cerr << "File '" << filename
                  << "' could not be writen to." << std::endl;
        exit(1);
    }

    std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType> sigma;
    typename SurfaceMesh::KeyType cnt = 1;
    for (const auto &x : mesh.get_level_id<1>())
    {
        sigma[mesh.get_name(x)[0]] = cnt++;
    }

    fout.precision(10);
    // Get the vertex data directly
    for (const auto &vertex : mesh.get_level<1>())
    {
        fout << "v "
             << vertex[0] << " "
             << vertex[1] << " "
             << vertex[2] << " "
             << "\n";
    }

    // Get the face nodes
    for (auto faceNodeID : mesh.get_level_id<3>())
    {
        auto w = mesh.get_name(faceNodeID);

        auto orientation = (*faceNodeID).orientation;
        if (orientation == 1)
        {
            fout << "f " << sigma[w[2]] << " " << sigma[w[1]] << " " << sigma[w[0]] << "\n";
        }
        else if (orientation == -1)
        {
            fout << "f " << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";

        }
        else
        {
            std::cerr << "Warning: Orientation undefined..." << std::endl;
            fout << "f " << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";
        }
    }
    fout.close();

}
} // end namespace gamer
