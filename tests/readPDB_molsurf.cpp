
#include "gamer.h"
#include <iostream>
#include <regex>

int main(int argc, char *argv[])
{
    std::string filename(argv[1]);

    SurfaceMesh* mesh = std::get<0>(SurfaceMesh::readPDB_molsurf(filename.c_str()));

    mesh->writeOFF(const_cast<char*>("1CID_01.off"));

    std::cout << mesh->numVertices() << std::endl;
    std::cout << mesh->numFaces() << std::endl;

    if(mesh->numVertices() != 48734)
    {
        return -1;
    }
    if(mesh->numFaces() != 97464)
    {
        return -1;
    }

    SurfaceMesh_ASC* amesh = std::get<1>(SurfaceMesh::readPDB_molsurf(filename.c_str()));

    std::cout << amesh->size<1>() << std::endl;
    std::cout << amesh->size<3>() << std::endl;

    if(amesh->size<1>() != 48734)
    {
        return -1;
    }
    if(amesh->size<3>() != 97464)
    {
        return -1;
    }

/*
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
            x *= vertex.x;
            y *= vertex.y;
            z *= vertex.z;
        }
        hash += sqrt(sqrt(fabs(x))*sqrt(fabs(y))*sqrt(fabs(z)));
    }

    std::cout << fabs(hash - hash_correct) / (hash_correct + 1e-16) << std::endl;
    if(fabs(hash - hash_correct) / (hash_correct + 1e-16) > tol)
    {
        std::cout << "Hash incorrect" << std::endl;
        return -1;
    }
    */
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
