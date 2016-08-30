
#include "gamer.h"
#include <iostream>
#include <regex>

int main(int argc, char *argv[])
{
    char* filename = argv[1];

    float  max_density;
    float *dataset;
    int    xdim, ydim, zdim;
    int    atom_num  = 0;
    ATOM  *atom_list = NULL;
    float  min[3], max[3];
    char   IsXYZR = 0;

    max_density = PDB2Volume(filename, &dataset, &xdim, &ydim, &zdim, min, max,
                             &atom_list, &atom_num, IsXYZR);


    std::cout << atom_num << std::endl;


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
*/
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
        /*
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
        */
//    SurfaceMesh* mesh = SurfaceMesh::readPDB_gauss(filename, -0.2, 2.5);

//    std::cerr << mesh->numVertices() << std::endl;

    return -1;
}
