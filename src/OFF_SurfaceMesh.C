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

//http://www.geomview.org/docs/html/OFF.html
std::pair<SurfaceMesh*, bool> readOFF(const std::string& filename)
{
    std::ifstream fin(filename);
    if(!fin.is_open())
    {
        std::cerr << "Read Error: File '" << filename << "' could not be read." << std::endl;
        return std::make_pair(nullptr, false);
    }
  
    std::string line;
    // Parse the first line:
    // [ST][C][N][4][n]OFF  # Header keyword
    // we only read OFF's in 3 space... simplicial_complex only does triangles
    getline(fin, line);  
    std::cout << line << std::endl;
    if (!line.compare("OFF\n")){
        std::cerr << "File Format Error: File '" << filename << "' does not look like a valid OFF file" << std::endl;
        return std::make_pair(nullptr, false);
    }

    int dimension = 3;

    // Lambda function to split the string
    auto split = [](const std::string& str, char delim = ' ') -> std::vector<std::string>{
        std::vector<std::string> result;
        auto begin = str.begin();
        do{
            auto end = begin;
            while(*end != delim && end != str.end())
                end++;
            if(end != begin) 
                result.push_back(std::string(begin,end));
            begin = end;
        } while (begin++ != str.end());  
        return result;
    };

    std::vector<std::string> arr; 
    // Allow some comments denoted by #...
    while(getline(fin,line)){
        if (line.find("#") == 0) {
            continue;
        }
        arr = split(line, '#');
        if(arr[0].length() != 0) // Assume that comments are over
            break;
    }
    
    // Parse the second line:
    // NVertices  NFaces  NEdges
    arr = split(arr[0]);
    int numVertices = std::stoi(arr[0]);
    std::cout << "#Vertices: " << numVertices << std::endl;

    int numFaces = std::stoi(arr[1]);
    std::cout << "#Faces: " << numFaces << std::endl;

    int numEdges = std::stoi(arr[2]);
    std::cout << "#Edges: " << numEdges << std::endl;

    // NOTE:: Assume there are no more comments...

    // Instantiate mesh!
    auto mesh = new SurfaceMesh();
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
    for(int i=0; i < numVertices; i++){
        getline(fin, line);
        std::cout << line << std::endl;
        arr = split(line, ' ');
        if(arr.size() < dimension){
            std::cerr << "Parse Error: Vertex line has fewer dimensions than expected (" << dimension << ")" << std::endl;
            delete mesh;
            return std::make_pair(nullptr, false);
        }
        double x = std::stod(arr[0]);
        double y = std::stod(arr[1]);
        double z = std::stod(arr[2]);
        mesh->insert<1>({i},Vertex(x,y,z));
    }

    // Parse Faces
    /*
        # Faces
        # Nv = # vertices on this face
        # v[0] ... v[Nv-1]: vertex indices
        #       in range 0..NVertices-1
    */
    for(int i=0; i < numFaces; i++){
        getline(fin, line);
        std::cout << line << std::endl;
        arr = split(line, ' ');
        if(std::stoi(arr[0]) != 3 && arr.size() < dimension+1){
            std::cerr << "Unsupported: Found face that is not a triangle!" << std::endl;
            delete mesh;
            return std::make_pair(nullptr, false);
        }
        auto v0 = std::stoi(arr[1]);
        auto v1 = std::stoi(arr[2]);
        auto v2 = std::stoi(arr[3]);
        mesh->insert<3>({v0,v1,v2});
    }
    fin.close();
    return std::make_pair(mesh, true);

/*    
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
*/
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

