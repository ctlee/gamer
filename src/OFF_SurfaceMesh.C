#include "SurfaceMesh.h"
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <ostream>
#include <regex>
#include <cmath>
#include <vector>

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
    return round(r * 10) * 121 + round(g * 10) * 11 + round(b * 10);
}

//http://www.geomview.org/docs/html/OFF.html
std::unique_ptr<SurfaceMesh> readOFF(const std::string& filename)
{
    // Instantiate mesh!
    //auto mesh = new SurfaceMesh();
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    std::ifstream fin(filename);
    if(!fin.is_open())
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
    if(!(line.find("OFF", line.size()-4) == std::string::npos)){
        std::cerr << "File Format Error: File '" << filename << "' does not look like a valid OFF file." << std::endl;
        std::cerr << "Expected 'OFF' at end of line, found: '" << line << "'." << std::endl;
        mesh.reset();
        return mesh;
    }

    // Have the support for reading in various things. Currently we are ignoring them though...
    bool textureCoords = false;
    if(!(line.find("ST") == std::string::npos)){
        std::cout << "Found vertex texture coordinates flag." << std::endl;
        textureCoords = true;
    }
    bool vertexColors = false;
    if(!(line.find("C") == std::string::npos)){
        std::cout << "Found vertex colors flag." << std::endl;
        vertexColors = true;
    }
    bool vertexNormals = false;
    if(!(line.find("N") == std::string::npos)){
        std::cout << "Found vertex normals flag." << std::endl;
        vertexNormals = true;
    }
    int dimension = 3;
    if(!(line.find("4") == std::string::npos)){
        std::cout << "Found dimension flag." << std::endl;
        dimension = 4;
    }

    if(!(line.find("n") == std::string::npos)){
        getline(fin,line);  // Assume no comments here yet...
        int tempDim = std::stoi(line);
        if (dimension == 4)
            dimension = tempDim + 1;
        else
            dimension = tempDim;
    }
  
    // Lambda function to split the string
    auto split = [](const std::string& cstr, std::vector<char> delim = {' ','\t'}) -> std::vector<std::string>{
        std::string str = cstr;
        for (auto i=1; i< delim.size(); i++){
            std::replace(str.begin(), str.end(),delim[i], delim[0]);
        }
        std::vector<std::string> result;
        auto begin = str.begin();
        do{
            auto end = begin;
            while(*end != delim[0] && end != str.end())
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
        arr = split(line, {'#'});
        if(arr[0].length() != 0) // Assume that comments are over
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
    for(int i=0; i < numVertices; i++){
        getline(fin, line);
        //std::cout << line << std::endl;
        arr = split(line);
        if(arr.size() < dimension){
            std::cerr << "Parse Error: Vertex line has fewer dimensions than expected (" << dimension << ")" << std::endl;
            mesh.reset();
            return mesh;
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
        arr = split(line);
        if(std::stoi(arr[0]) != 3 && arr.size() < dimension+1){
            std::cerr << "Unsupported: Found face that is not a triangle!" << std::endl;
            mesh.reset();
            return mesh;
        }
        else if(arr.size() == dimension+1){
            auto v0 = std::stoi(arr[1]);
            auto v1 = std::stoi(arr[2]);
            auto v2 = std::stoi(arr[3]);
            mesh->insert<3>({v0,v1,v2});
        }
        else if(arr.size() == dimension+5){
            auto v0 = std::stoi(arr[1]);
            auto v1 = std::stoi(arr[2]);
            auto v2 = std::stoi(arr[3]);
            // parse for marker/color data
            auto r = std::stod(arr[4]);
            auto g = std::stod(arr[5]);
            auto b = std::stod(arr[6]);
            //auto k = std::stod(arr[7]);
            mesh->insert<3>({v0,v1,v2},Face(Orientable{0}, FaceProperties{get_marker(r,g,b),0}));
        }
        else {
            std::cerr << "Parse Error: Couldn't interpret face: '" << line << "'." << std::endl;
            mesh.reset();
            return mesh;
        }
    }
    fin.close();
    return mesh;
}

void writeOFF(const std::string& filename, const SurfaceMesh& mesh){
    std::ofstream fout(filename);
    if(!fout.is_open())
    {
        std::cerr   << "File '" << filename 
                    << "' could not be writen to." << std::endl;
        exit(1); 
    }

    std::map<typename SurfaceMesh::KeyType,typename SurfaceMesh::KeyType> sigma;
    typename SurfaceMesh::KeyType cnt = 0;
    for(const auto& x : mesh.get_level_id<1>())
    {
        sigma[mesh.get_name(x)[0]] = cnt++;
    }

    fout << "OFF" << std::endl;

    int numVertices = mesh.size<1>();
    int numFaces = mesh.size<3>();
    int numEdges = mesh.size<2>();
    fout    << numVertices << " " 
            << numFaces << " "
            << numEdges << "\n"; 

    fout.precision(10); 
    // Get the vertex data directly 
    // TODO: this actually has to print out the vertices in order of the index...
    for(const auto& vertex : mesh.get_level<1>()){
        fout    << vertex[0] << " "
                << vertex[1] << " " 
                << vertex[2] << " " 
                << "\n";
    }

    // Get the face nodes
    for(auto faceNodeID : mesh.get_level_id<3>()){
        auto w = mesh.get_name(faceNodeID);

        auto orientation = (*faceNodeID).orientation;
        if (orientation == 1){
            fout << "3 " << sigma[w[2]] << " " << sigma[w[1]] << " " << sigma[w[0]] << "\n";
        }
        else if(orientation == -1){
            fout << "3 " << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";

        }
        else{
            std::cerr << "~~Orientation undefined..." << std::endl;
            fout << "3 " << sigma[w[0]] << " " << sigma[w[1]] << " " << sigma[w[2]] << "\n";
        }
    }
    fout.close();
}

