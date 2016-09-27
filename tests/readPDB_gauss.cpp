
#include "gamer.h"
#include <iostream>
#include <regex>

int main(int argc, char *argv[])
{
    std::string filename(argv[1]);

    float  max_density;
    float *dataset;
    int    xdim, ydim, zdim;
    std::vector<ATOM>  atom_list;
    float  min[3], max[3];

    float max_density_correct = 11.02886;
    float min_correct[3] = {-40.1917, 5.502299, -7.506701};
    float max_correct[3] = {21.3867, 68.6437, 74.137703};
    float tol = 1e-4;

    max_density = PDB2Volume(filename, &dataset, &xdim, &ydim, &zdim, min, max,
                             atom_list);

    if(fabs(max_density - max_density_correct) > tol)
    {
        return -1;
    }
    if(fabs(min[0] - min_correct[0]) > tol)
    {
        return -1;
    }
    if(fabs(min[1] - min_correct[1]) > tol)
    {
        return -1;
    }
    if(fabs(min[2] - min_correct[2]) > tol)
    {
        return -1;
    }
    if(fabs(max[0] - max_correct[0]) > tol)
    {
        return -1;
    }
    if(fabs(max[1] - max_correct[1]) > tol)
    {
        return -1;
    }
    if(fabs(max[2] - max_correct[2]) > tol)
    {
        return -1;
    }

    std::cout << atom_list.size() << std::endl;



    SurfaceMesh_ASC* mesh = std::get<1>(SurfaceMesh::readPDB_gauss(filename.c_str(), 0.0, 2.5));
/*
    mesh->writeOFF(const_cast<char*>("1CID_01.off"));

    if(mesh->numVertices() != 45532)
    {
        return -1;
    }
    if(mesh->numFaces() != 91068)
    {
        return -1;
    }
*/
    double hash = 0.0;
    double hash_correct = 1.15136e8;
    for(auto face : mesh->get_level_id<3>())
    {
        double x = 1;
        double y = 1;
        double z = 1;
        for(auto a : mesh->get_name(face))
        {
            auto vertex = mesh->get<1>({a});
            x *= vertex.position.get(0,0);
            y *= vertex.position.get(1,0);
            z *= vertex.position.get(2,0);
        }
        hash += sqrt(sqrt(fabs(x))*sqrt(fabs(y))*sqrt(fabs(z)));
    }

    std::cout << fabs(hash - hash_correct) / (hash_correct + 1e-16) << std::endl;
    if(fabs(hash - hash_correct) / (hash_correct + 1e-16) > tol)
    {
        std::cout << "Hash incorrect" << std::endl;
//        return -1;
    }


    SurfaceMesh* smesh = std::get<0>(SurfaceMesh::readPDB_gauss(filename.c_str(), 0.0, 2.5));
    smesh->writeOFF(const_cast<char*>("1CID_01.off"));
/*
    for(int i = 0; i < atom_num; ++i)
    {
        std::cout << atom_list[i].radius << std::endl;
    }

    for(int i = 0; i < xdim; ++i)
    {
        for(int j = 0; j < ydim; ++j)
        {
            for(int k = 0; k < zdim; ++k)
            {
                std::cout << dataset[i + xdim*(j + ydim*k)] << std::endl;
            }
        }
    }

    std::cout << xdim << " " << ydim << " " << zdim << std::endl;

    std::cout << max_density << std::endl;

//    std::string atom_str("ATOM(?:\\s*([^\\s]*))*");
    const char* line = "ATOM 0 1 2 3 4 5 6 7 ";

    //std::string atom_str("ATOM(\\s+([^\\s]+))*.*");
    //std::string atom_str("ATOM(?:\\s+([^\\s]+))(?:\\s+([^\\s]+))");
    std::regex atom("ATOM\\s+(.*)", std::regex::icase);
    std::cmatch match;

    if(std::regex_match(line, match, atom))
    {
        std::regex strseparator("\\s+", std::regex::icase);
        std::string data = match.str(1);

        std::sregex_token_iterator next(data.begin(), data.end(), strseparator, -1);
        std::sregex_token_iterator end;

        std::vector<std::string> tokens;

        std::copy(next, end, std::back_inserter(tokens));

        for(std::string s : tokens)
        {
            std::cout << s << std::endl;
        }
    }

    std::string test1(line+5, line+8);
    std::cout << "> " << test1 << std::endl;    

    while (next != end)
    {
        std::cout << *next << std::endl;
        std::smatch match = *next;
        std::cout << line << std::endl;
        std::cout << ">" << match.str(1) << "|" << match.str(2) << "<" << std::endl;
        std::cout << match.size() << std::endl;
        std::cout << match.str(6) << " " << match.str(7) << " " << match.str(8) << std::endl;
        ++next;
    }

//    SurfaceMesh* mesh = SurfaceMesh::readPDB_gauss(filename, -0.2, 2.5);

//    std::cerr << mesh->numVertices() << std::endl;
*/
    return 0;
}
