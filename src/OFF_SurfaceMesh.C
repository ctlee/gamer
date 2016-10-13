#include "SurfaceMesh.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <regex>
#include <cmath>

int get_marker(float r, float g, float b)
{
    if (r < 0 || r > 1 || g < 0 || g > 1 || b < 0 || b > 1)
    {
        std::cerr << "Expected individual RGB value to be betwen 0 and 1." << std::endl;
        exit(1);
    }
    return round(r * 10) * 121 + round(g * 10) * 11 + round(b * 10);
}


SurfaceMesh* readOFF(const std::string filename)
{
	auto F = new SurfaceMesh();
    int n, m;
    int a, b, c;
    float x, y, z, color_r, color_g, color_b, color_a;
    unsigned int          v_n, t_n, e_n;
    char character;
    std::string  line;
    bool         has_marker = false;
    const char* input_name = filename.c_str();


    std::fstream fin(filename);
    if (!fin.is_open())
    {
        std::cerr << "Read error. File \'" << input_name << "\' could not be read." << std::endl;
        return nullptr;
    }

    while (getline(fin, line))
    {
        if ((line[0] == 'O') && (line[1] == 'F') && (line[2] == 'F'))
        {
            break;
        }
    }

    if (!getline(fin, line) || scanf(line.c_str(), "%d %d %d\n", &v_n, &t_n, &e_n) != 3)
    {

        std::cerr << "Read error. Expected 3 integer number for the number of vertices and simplices." << std::endl;
        return NULL;
    }

    std::cout << "   vertices: " << v_n << " --- simplices: " << t_n << std::endl;

    // Read vertex coordinates
    for (n = 0; n < v_n; n++) // surfmesh->num_vertices; n++) {
    {
        if (!getline(fin, line) || scanf(line.c_str(), "%f %f %f\n", &x, &y, &z) != 3)
        {
            std::cerr << "Read error. Expected 3 floats for the coordinates of vertex:" << n << std::endl;
            return NULL;
        }

        //        surfmesh->vertex.push_back(FLTVECT(x,y,z));

        Vertex v{x,y,z};  
        F->insert<1>({n}, v);
    }

    // Read the number of vertices per simplex from the first line of simplices
    getline(fin, line);
    n = scanf(line.c_str(), "%d", &m);

    std::cout << m << std::endl;

    // Check format of the read simplices
    if ((m != 3) && (m != 4))
    {
        std::cerr << "Read error. Expected a 3 or 4 for the first value in the first simplex line." << std::endl;
        return NULL;
    }

    // Input is surface mesh
    else if (m == 3)
    {
        std::cout << "   Input is surface mesh." << std::endl;

        // Read the rest of the first line
        if (!getline(fin, line) || scanf(line.c_str(), "%d %d %d", &a, &b, &c) != 3)
        {
            std::cerr << "Read error. Expected 3 integers for the first simplex." << std::endl;
            return NULL;
        }

        // Read any blanks
        do {
            fin.get(character);
        } while(character == ' ');

        // If we do not have en end of line character we have marker
        if (character != '\n')
        {
            // Set marker flag
            has_marker = true;

            // Read first rgb values
            if (!getline(fin, line) || scanf(line.c_str(), "%f %f %f %f", &color_r, &color_g, &color_b, &color_a) != 4)
            {
                std::cerr << "Read error. Expected 4 floats for the RGBA values of the first face." << std::endl;
                return NULL;
            }

            while (fin.get(character) && character != '\n')
            {}

        	F->insert<3>({a,b,c}, Face(Orientable{0}, FaceProperties{get_marker(color_r, color_g, color_b),0}));
        }
        else
        {
	        F->insert<3>({a,b,c});
        }

        // Read the rest of the simplices
        for (n = 1; n < t_n; n++)
        {
            if (!getline(fin, line) || scanf(line.c_str(), "%d %d %d %d", &m, &a, &b, &c) != 4)
            {
                std::cerr << "Read error. Expected 4 integers for simplex " << n << std::endl;
                return nullptr;
            }

            // If we have face markers
            if (has_marker)
            {
                // Read first rgb values
                if (!getline(fin, line) || scanf(line.c_str(), "%f %f %f %f", &color_r, &color_g, &color_b,
                           &color_a) != 4)
                {
                    std::cerr << "Read error. Expected 4 floats for the RGBA values of face:" << n << std::endl;
                    return NULL;
                }

                // Get the other face markers
                //surfmesh->face[n].m = get_marker(color_r, color_g, color_b);
	        	F->insert<3>({a,b,c}, Face(Orientable{0}, FaceProperties{get_marker(color_r, color_g, color_b),0}));
            }
            else
            {
		        F->insert<3>({a,b,c});
            }

            // Skip any additional character on this line
            while (fin.get(character) && character != '\n')
            {}
        }

    }

    fin.close();

    return F;
}

void writeOFF(const std::string& filename, const SurfaceMesh& mesh){

    std::ofstream fout(filename);
    if(!fout.is_open())
    {
        std::cerr   << "File '" << filename 
                    << "' could not be writen to." << std::endl;
        exit(1); 
    }

    fout << "OFF" << std::endl;

    int numVertices = mesh.size<1>();
    int numFaces = mesh.size<3>();
    fout << numVertices << numFaces << numVertices + numFaces - 2 << "\n"; 

    fout.precision(10); 
    // Get the vertex data directly 
    // TODO: this actually has to print out the vertices in order of the index...
    for(auto& vertex : mesh.get_level<1>()){
        fout    << vertex[0] << " "
                << vertex[1] << " " 
                << vertex[2] << " " 
                << "\n";
    }

    // Get the face nodes
    for(auto faceNodeID : mesh.get_level_id<3>()){
        auto w = mesh.get_name(faceNodeID);
        fout << "3 " << w[0] << " " << w[1] << " " << w[2] << "\n";
    }
    fout.close();
}

