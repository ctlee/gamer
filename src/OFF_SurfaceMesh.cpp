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


#include "gamer/SurfaceMesh.h"
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <ostream>
#include <regex>
#include <cmath>
#include <vector>

/// Namespace for all things gamer
namespace gamer
{
/**
 * @brief      Converts the color from float(0-1) to a marker value
 *
 * @param[in]  r     value of red
 * @param[in]  g     value of green
 * @param[in]  b     value of blue
 *
 * @return     The value of the marker.
 */
int get_marker(float r, float g, float b)
{
    if (r < 0 || r > 1 || g < 0 || g > 1 || b < 0 || b > 1)
    {
        std::cerr << "Expected individual RGB value to be betwen 0 and 1." << std::endl;
        exit(1);
    }
    return static_cast<int>(round(r * 10) * 121 + round(g * 10) * 11 + round(b * 10));
}

//http://www.geomview.org/docs/html/OFF.html
std::unique_ptr<SurfaceMesh> readOFF(const std::string &filename)
{
    // Instantiate mesh!
    //auto mesh = new SurfaceMesh();
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::ifstream                fin(filename);
    if (!fin.is_open())
    {
        std::cerr << "Read Error: File '" << filename << "' could not be read." << std::endl;
        mesh.reset();
        return mesh;
    }

    std::string line;
    // Parse the first line:
    // [ST][C][N][4][n]OFF  # Header keyword
    // we only read OFF's in 3 space... simplicial_complex only does triangles
    getline(fin, line);
    // Assume the first line must end with 'OFF\n'
    if (!(line.find("OFF", line.size()-4) == std::string::npos))
    {
        std::cerr << "File Format Error: File '" << filename << "' does not look like a valid OFF file." << std::endl;
        std::cerr << "Expected 'OFF' at end of line, found: '" << line << "'." << std::endl;
        mesh.reset();
        return mesh;
    }

    // Have the support for reading in various things. Currently we are ignoring
    // them though...
    bool textureCoords = false;
    if (!(line.find("ST") == std::string::npos))
    {
        std::cout << "Found vertex texture coordinates flag." << std::endl;
        textureCoords = true;
    }
    bool vertexColors = false;
    if (!(line.find("C") == std::string::npos))
    {
        std::cout << "Found vertex colors flag." << std::endl;
        vertexColors = true;
    }
    bool vertexNormals = false;
    if (!(line.find("N") == std::string::npos))
    {
        std::cout << "Found vertex normals flag." << std::endl;
        vertexNormals = true;
    }
    int dimension = 3;
    if (!(line.find("4") == std::string::npos))
    {
        std::cout << "Found dimension flag." << std::endl;
        dimension = 4;
    }

    if (!(line.find("n") == std::string::npos))
    {
        getline(fin, line);  // Assume no comments here yet...
        int tempDim = std::stoi(line);
        if (dimension == 4)
            dimension = tempDim + 1;
        else
            dimension = tempDim;
    }

    // Lambda function to split the string
    auto split = [](const std::string &cstr, std::vector<char> delim = {' ', '\t'}) -> std::vector<std::string> {
            std::string str = cstr;
            for (auto i = 1; i < delim.size(); i++)
            {
                std::replace(str.begin(), str.end(), delim[i], delim[0]);
            }
            std::vector<std::string> result;
            auto                     begin = str.begin();
            do
            {
                auto end = begin;
                while (*end != delim[0] && end != str.end())
                    end++;
                if (end != begin)
                    result.push_back(std::string(begin, end));
                begin = end;
            }
            while (begin++ != str.end());
            return result;
        };

    std::vector<std::string> arr;
    // Allow some comments denoted by #...
    while (getline(fin, line))
    {
        if (line.find("#") == 0)
        {
            continue;
        }
        arr = split(line, {'#'});
        if (arr[0].length() != 0) // Assume that comments are over
            break;
    }

    // Parse the second line:
    // NVertices  NFaces  NEdges
    arr = split(arr[0]);
    int numVertices = std::stoi(arr[0]);
    //std::cout << "#Vertices: " << numVertices << std::endl;

    int numFaces = std::stoi(arr[1]);
    //std::cout << "#Faces: " << numFaces << std::endl;

    // numEdges is ignored
    //int numEdges = std::stoi(arr[2]);
    //std::cout << "#Edges: " << numEdges << std::endl;

    // NOTE:: Assume there are no more comments...

    // Parse the vertices
    /*
       x[0]  y[0]  z[0]
     # Vertices, possibly with normals,
     # colors, and/or texture coordinates, in that order,
     # if the prefixes N, C, ST
     # are present.
     # If 4OFF, each vertex has 4 components,
     # including a final homogeneous component.
     # If nOFF, each vertex has Ndim components.
     # If 4nOFF, each vertex has Ndim+1 components.
     */
    for (int i = 0; i < numVertices; i++)
    {
        getline(fin, line);
        //std::cout << line << std::endl;
        arr = split(line);
        if (arr.size() < dimension)
        {
            std::cerr << "Parse Error: Vertex line has fewer dimensions than expected (" << dimension << ")" << std::endl;
            mesh.reset();
            return mesh;
        }
        double x = std::stod(arr[0]);
        double y = std::stod(arr[1]);
        double z = std::stod(arr[2]);
        mesh->insert<1>({i}, SMVertex(x, y, z));
    }

    // Parse Faces
    /*
     # Faces
     # Nv = # vertices on this face
     # v[0] ... v[Nv-1]: vertex indices
     #       in range 0..NVertices-1

        TODO: This should parse the number of vertices in each face (10)
        TODO: (10) Get orientation from the OFF file
     */
    for (int i = 0; i < numFaces; i++)
    {
        getline(fin, line);
        arr = split(line);
        if (std::stoi(arr[0]) != 3 && arr.size() < dimension+1)
        {
            std::cerr << "Unsupported: Found face that is not a triangle!" << std::endl;
            mesh.reset();
            return mesh;
        }
        else if (arr.size() == dimension+1)
        {
            auto v0 = std::stoi(arr[1]);
            auto v1 = std::stoi(arr[2]);
            auto v2 = std::stoi(arr[3]);
            mesh->insert<3>({v0, v1, v2});
        }
        else if (arr.size() == dimension+5)
        {
            auto v0 = std::stoi(arr[1]);
            auto v1 = std::stoi(arr[2]);
            auto v2 = std::stoi(arr[3]);
            // parse for marker/color data
            auto r = std::stod(arr[4]);
            auto g = std::stod(arr[5]);
            auto b = std::stod(arr[6]);
            //auto k = std::stod(arr[7]);
            mesh->insert<3>({v0, v1, v2}, SMFace(get_marker(r, g, b), false));
        }
        else
        {
            std::cerr << "Parse Error: Couldn't interpret face: '" << line << "'." << std::endl;
            mesh.reset();
            return mesh;
        }
    }
    fin.close();
    compute_orientation(*mesh);
    return mesh;
}

void writeOFF(const std::string &filename, const SurfaceMesh &mesh)
{
    std::ofstream fout(filename);
    if (!fout.is_open())
    {
        std::cerr << "File '" << filename
                  << "' could not be writen to." << std::endl;
        exit(1);
    }

    fout << "OFF" << std::endl;

    std::size_t numVertices = mesh.size<1>();
    std::size_t numFaces = mesh.size<3>();
    std::size_t numEdges = mesh.size<2>();
    fout << numVertices << " "
         << numFaces << " "
         << numEdges << "\n";

    std::map<typename SurfaceMesh::KeyType, typename SurfaceMesh::KeyType> sigma;
    typename SurfaceMesh::KeyType cnt = 0;

    fout.precision(10);
    // Get the vertex data directly
    for (const auto vertexID : mesh.get_level_id<1>())
    {
        sigma[mesh.get_name(vertexID)[0]] = cnt++;
        auto vertex = *vertexID;

        fout << vertex[0] << " "
             << vertex[1] << " "
             << vertex[2] << " "
             << "\n";
    }

    bool orientationError = false;

    // Get the face nodes
    for (auto faceNodeID : mesh.get_level_id<3>())
    {
        auto w = mesh.get_name(faceNodeID);

        auto orientation = (*faceNodeID).orientation;
        if (orientation == 1)
        {
            fout << "3 " << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";
        }
        else if (orientation == -1)
        {
            fout << "3 " << sigma[w[2]] << " " << sigma[w[1]] << " " << sigma[w[0]] << "\n";

        }
        else
        {
            orientationError = true;
            fout << "3 " << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";
        }
    }
    if (orientationError)
    {
        std::cerr << "WARNING(writeOFF): The orientation of one or more faces "
                  << "is not defined. Did you run compute_orientation()?"
                  << std::endl;
    }
    fout.close();
}
} // end namespace gamer
