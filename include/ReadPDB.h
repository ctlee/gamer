#include <algorithm>
#include <regex>
#include <iostream>
#include <string>
#include <fstream>
#include <array>


namespace detail
{
    const float EPSILON    = 1e-3f;
//    const float BLOBBYNESS = 1e-3f;
    static const std::regex PDB(".*.pdb", std::regex::icase | std::regex::optimize);
    static const std::regex PQR(".*.pqr", std::regex::icase | std::regex::optimize);
    static const std::regex XYZR(".*.xyzr", std::regex::icase | std::regex::optimize);
    static const std::regex atom("ATOM.*\n*", std::regex::optimize);

    struct PDBelementInformation 
    {
      const char*           atomName;
      const char*           residueName;
      float         radius;
      float         red;
      float         green;
      float         blue;
      int           hydrophobicity;
      unsigned char         residueIndex;
    };

    static const std::size_t MAX_BIOCHEM_ELEMENTS = 167;
    static PDBelementInformation PDBelementTable[MAX_BIOCHEM_ELEMENTS] =
    {
         {" N  ", "GLY", 1.625f, 0.0f, 0.0f, 1.0f,  1, 10 },
         {" CA ", "GLY", 1.750f, 0.3f, 0.3f, 0.3f, -1, 10 },
         {" C  ", "GLY", 1.875f, 0.3f, 0.3f, 0.3f,  1, 10 },
         {" O  ", "GLY", 1.480f, 1.0f, 0.0f, 0.0f,  1, 10 },
         {" N  ", "ALA", 1.625f, 0.0f, 0.0f, 1.0f,  1,  1 },
         {" CA ", "ALA", 1.750f, 0.3f, 0.3f, 0.3f, -1,  1 },
         {" C  ", "ALA", 1.875f, 0.3f, 0.3f, 0.3f,  1,  1 },
         {" O  ", "ALA", 1.480f, 1.0f, 0.0f, 0.0f,  1,  1 },
         {" CB ", "ALA", 1.750f, 0.3f, 0.3f, 0.3f, -1,  1 },
         {" N  ", "VAL", 1.625f, 0.0f, 0.0f, 1.0f,  1, 23 },
         {" CA ", "VAL", 1.750f, 0.3f, 0.3f, 0.3f, -1, 23 },
         {" C  ", "VAL", 1.875f, 0.3f, 0.3f, 0.3f,  1, 23 },
         {" O  ", "VAL", 1.480f, 1.0f, 0.0f, 0.0f,  1, 23 },
         {" CB ", "VAL", 1.750f, 0.3f, 0.3f, 0.3f, -1, 23 },
         {" CG1", "VAL", 1.750f, 0.3f, 0.3f, 0.3f, -1, 23 },
         {" CG2", "VAL", 1.750f, 0.3f, 0.3f, 0.3f, -1, 23 },
         {" N  ", "LEU", 1.625f, 0.0f, 0.0f, 1.0f,  1, 13 },
         {" CA ", "LEU", 1.750f, 0.3f, 0.3f, 0.3f, -1, 13 },
         {" C  ", "LEU", 1.875f, 0.3f, 0.3f, 0.3f,  1, 13 },
         {" O  ", "LEU", 1.480f, 1.0f, 0.0f, 0.0f,  1, 13 },
         {" CB ", "LEU", 1.750f, 0.3f, 0.3f, 0.3f, -1, 13 },
         {" CG ", "LEU", 1.750f, 0.3f, 0.3f, 0.3f, -1, 13 },
         {" CD1", "LEU", 1.750f, 0.3f, 0.3f, 0.3f, -1, 13 },
         {" CD2", "LEU", 1.750f, 0.3f, 0.3f, 0.3f, -1, 13 },
         {" N  ", "ILE", 1.625f, 0.0f, 0.0f, 1.0f,  1, 12 },
         {" CA ", "ILE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 12 },
         {" C  ", "ILE", 1.875f, 0.3f, 0.3f, 0.3f,  1, 12 },
         {" O  ", "ILE", 1.480f, 1.0f, 0.0f, 0.0f,  1, 12 },
         {" CB ", "ILE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 12 },
         {" CG1", "ILE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 12 },
         {" CG2", "ILE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 12 },
         {" CD1", "ILE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 12 },
         {" N  ", "MET", 1.625f, 0.0f, 0.0f, 1.0f,  1, 15 },
         {" CA ", "MET", 1.750f, 0.3f, 0.3f, 0.3f, -1, 15 },
         {" C  ", "MET", 1.875f, 0.3f, 0.3f, 0.3f,  1, 15 },
         {" O  ", "MET", 1.480f, 1.0f, 0.0f, 0.0f,  1, 15 },
         {" CB ", "MET", 1.750f, 0.3f, 0.3f, 0.3f, -1, 15 },
         {" CG ", "MET", 1.750f, 0.3f, 0.3f, 0.3f, -1, 15 },
         {" SD ", "MET", 1.775f, 1.0f, 1.0f, 0.0f,  1, 15 },
         {" CE ", "MET", 1.750f, 0.3f, 0.3f, 0.3f, -1, 15 },
         {" N  ", "PRO", 1.625f, 0.0f, 0.0f, 1.0f, -1, 17 },
         {" CA ", "PRO", 1.750f, 0.3f, 0.3f, 0.3f, -1, 17 },
         {" C  ", "PRO", 1.875f, 0.3f, 0.3f, 0.3f,  1, 17 },
         {" O  ", "PRO", 1.480f, 1.0f, 0.0f, 0.0f,  1, 17 },
         {" CB ", "PRO", 1.750f, 0.3f, 0.3f, 0.3f, -1, 17 },
         {" CG ", "PRO", 1.750f, 0.3f, 0.3f, 0.3f, -1, 17 },
         {" CD ", "PRO", 1.750f, 0.3f, 0.3f, 0.3f, -1, 17 },
         {" N  ", "PHE", 1.625f, 0.0f, 0.0f, 1.0f,  1, 16 },
         {" CA ", "PHE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" C  ", "PHE", 1.875f, 0.3f, 0.3f, 0.3f,  1, 16 },
         {" O  ", "PHE", 1.480f, 1.0f, 0.0f, 0.0f,  1, 16 },
         {" CB ", "PHE", 1.750f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" CG ", "PHE", 1.775f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" CD1", "PHE", 1.775f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" CD2", "PHE", 1.775f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" CE1", "PHE", 1.775f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" CE2", "PHE", 1.775f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" CZ ", "PHE", 1.775f, 0.3f, 0.3f, 0.3f, -1, 16 },
         {" N  ", "TRP", 1.625f, 0.0f, 0.0f, 1.0f,  1, 20 },
         {" CA ", "TRP", 1.750f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" C  ", "TRP", 1.875f, 0.3f, 0.3f, 0.3f,  1, 20 },
         {" O  ", "TRP", 1.480f, 1.0f, 0.0f, 0.0f,  1, 20 },
         {" CB ", "TRP", 1.750f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CG ", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CD1", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CD2", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" NE1", "TRP", 1.625f, 0.2f, 0.2f, 1.0f,  1, 20 },
         {" CE2", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CE3", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CZ2", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CZ3", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" CH2", "TRP", 1.775f, 0.3f, 0.3f, 0.3f, -1, 20 },
         {" N  ", "SER", 1.625f, 0.0f, 0.0f, 1.0f,  1, 18 },
         {" CA ", "SER", 1.750f, 0.3f, 0.3f, 0.3f, -1, 18 },
         {" C  ", "SER", 1.875f, 0.3f, 0.3f, 0.3f,  1, 18 },
         {" O  ", "SER", 1.480f, 1.0f, 0.0f, 0.0f,  1, 18 },
         {" CB ", "SER", 1.750f, 0.3f, 0.3f, 0.3f, -1, 18 },
         {" OG ", "SER", 1.560f, 1.0f, 0.0f, 0.0f,  1, 18 },
         {" N  ", "THR", 1.625f, 0.0f, 0.0f, 1.0f,  1, 19 },
         {" CA ", "THR", 1.750f, 0.3f, 0.3f, 0.3f, -1, 19 },
         {" C  ", "THR", 1.875f, 0.3f, 0.3f, 0.3f,  1, 19 },
         {" O  ", "THR", 1.480f, 1.0f, 0.0f, 0.0f,  1, 19 },
         {" CB ", "THR", 1.750f, 0.3f, 0.3f, 0.3f, -1, 19 },
         {" OG1", "THR", 1.560f, 1.0f, 0.0f, 0.0f,  1, 19 },
         {" CG2", "THR", 1.750f, 0.3f, 0.3f, 0.3f, -1, 19 },
         {" N  ", "ASN", 1.625f, 0.0f, 0.0f, 1.0f,  1,  3 },
         {" CA ", "ASN", 1.750f, 0.3f, 0.3f, 0.3f, -1,  3 },
         {" C  ", "ASN", 1.875f, 0.3f, 0.3f, 0.3f,  1,  3 },
         {" O  ", "ASN", 1.480f, 1.0f, 0.0f, 0.0f,  1,  3 },
         {" CB ", "ASN", 1.750f, 0.3f, 0.3f, 0.3f, -1,  3 },
         {" CG ", "ASN", 1.875f, 0.3f, 0.3f, 0.3f,  1,  3 },
         {" OD1", "ASN", 1.480f, 1.0f, 0.0f, 0.0f,  1,  3 },
         {" ND2", "ASN", 1.625f, 0.2f, 0.2f, 1.0f,  1,  3 },
         {" N  ", "GLN", 1.625f, 0.0f, 0.0f, 1.0f,  1,  7 },
         {" CA ", "GLN", 1.750f, 0.3f, 0.3f, 0.3f, -1,  7 },
         {" C  ", "GLN", 1.875f, 0.3f, 0.3f, 0.3f,  1,  7 },
         {" O  ", "GLN", 1.480f, 1.0f, 0.0f, 0.0f,  1,  7 },
         {" CB ", "GLN", 1.750f, 0.3f, 0.3f, 0.3f, -1,  7 },
         {" CG ", "GLN", 1.750f, 0.3f, 0.3f, 0.3f, -1,  7 },
         {" CD ", "GLN", 1.875f, 0.3f, 0.3f, 0.3f,  1,  7 },
         {" OE1", "GLN", 1.480f, 1.0f, 0.0f, 0.0f,  1,  7 },
         {" NE2", "GLN", 1.625f, 0.2f, 0.2f, 1.0f,  1,  7 },
         {" N  ", "TYR", 1.625f, 0.0f, 0.0f, 1.0f,  1, 21 },
         {" CA ", "TYR", 1.750f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" C  ", "TYR", 1.875f, 0.3f, 0.3f, 0.3f,  1, 21 },
         {" O  ", "TYR", 1.480f, 1.0f, 0.0f, 0.0f,  1, 21 },
         {" CB ", "TYR", 1.750f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" CG ", "TYR", 1.775f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" CD1", "TYR", 1.775f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" CD2", "TYR", 1.775f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" CE1", "TYR", 1.775f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" CE2", "TYR", 1.775f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" CZ ", "TYR", 1.775f, 0.3f, 0.3f, 0.3f, -1, 21 },
         {" OH ", "TYR", 1.535f, 1.0f, 0.0f, 0.0f,  1, 21 },
         {" N  ", "CYS", 1.625f, 0.0f, 0.0f, 1.0f,  1,  6 },
         {" CA ", "CYS", 1.750f, 0.3f, 0.3f, 0.3f, -1,  6 },
         {" C  ", "CYS", 1.875f, 0.3f, 0.3f, 0.3f,  1,  6 },
         {" O  ", "CYS", 1.480f, 1.0f, 0.0f, 0.0f,  1,  6 },
         {" CB ", "CYS", 1.750f, 0.3f, 0.3f, 0.3f, -1,  6 },
         {" SG ", "CYS", 1.775f, 1.0f, 1.0f, 0.0f,  1,  6 },
         {" N  ", "LYS", 1.625f, 0.0f, 0.0f, 1.0f,  1, 14 },
         {" CA ", "LYS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 14 },
         {" C  ", "LYS", 1.875f, 0.3f, 0.3f, 0.3f,  1, 14 },
         {" O  ", "LYS", 1.480f, 1.0f, 0.0f, 0.0f,  1, 14 },
         {" CB ", "LYS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 14 },
         {" CG ", "LYS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 14 },
         {" CD ", "LYS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 14 },
         {" CE ", "LYS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 14 },
         {" NZ ", "LYS", 1.625f, 0.2f, 0.2f, 1.0f,  1, 14 },
         {" N  ", "ARG", 1.625f, 0.0f, 0.0f, 1.0f,  1,  2 },
         {" CA ", "ARG", 1.750f, 0.3f, 0.3f, 0.3f, -1,  2 },
         {" C  ", "ARG", 1.875f, 0.3f, 0.3f, 0.3f,  1,  2 },
         {" O  ", "ARG", 1.480f, 1.0f, 0.0f, 0.0f,  1,  2 },
         {" CB ", "ARG", 1.750f, 0.3f, 0.3f, 0.3f, -1,  2 },
         {" CG ", "ARG", 1.750f, 0.3f, 0.3f, 0.3f, -1,  2 },
         {" CD ", "ARG", 1.750f, 0.3f, 0.3f, 0.3f, -1,  2 },
         {" NE ", "ARG", 1.625f, 0.2f, 0.2f, 1.0f,  1,  2 },
         {" CZ ", "ARG", 1.125f, 0.3f, 0.3f, 0.3f,  1,  2 },
         {" NH1", "ARG", 1.625f, 0.2f, 0.2f, 1.0f,  1,  2 },
         {" NH2", "ARG", 1.625f, 0.2f, 0.2f, 1.0f,  1,  2 },
         {" N  ", "HIS", 1.625f, 0.0f, 0.0f, 1.0f,  1, 11 },
         {" CA ", "HIS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 11 },
         {" C  ", "HIS", 1.875f, 0.3f, 0.3f, 0.3f,  1, 11 },
         {" O  ", "HIS", 1.480f, 1.0f, 0.0f, 0.0f,  1, 11 },
         {" CB ", "HIS", 1.750f, 0.3f, 0.3f, 0.3f, -1, 11 },
         {" CG ", "HIS", 1.775f, 0.3f, 0.3f, 0.3f, -1, 11 },
         {" ND1", "HIS", 1.625f, 0.2f, 0.2f, 1.0f,  1, 11 },
         {" CD2", "HIS", 1.775f, 0.3f, 0.3f, 0.3f, -1, 11 },
         {" CE1", "HIS", 1.775f, 0.3f, 0.3f, 0.3f,  1, 11 },
         {" NE2", "HIS", 1.625f, 0.2f, 0.2f, 1.0f,  1, 11 },
         {" N  ", "ASP", 1.625f, 0.0f, 0.0f, 1.0f,  1,  4 },
         {" CA ", "ASP", 1.750f, 0.3f, 0.3f, 0.3f, -1,  4 },
         {" C  ", "ASP", 1.875f, 0.3f, 0.3f, 0.3f,  1,  4 },
         {" O  ", "ASP", 1.480f, 1.0f, 0.0f, 0.0f,  1,  4 },
         {" CB ", "ASP", 1.750f, 0.3f, 0.3f, 0.3f, -1,  4 },
         {" CG ", "ASP", 1.875f, 0.3f, 0.3f, 0.3f,  1,  4 },
         {" OD1", "ASP", 1.480f, 1.0f, 1.0f, 1.0f,  1,  4 },
         {" OD2", "ASP", 1.480f, 1.0f, 0.0f, 0.0f,  1,  4 },
         {" N  ", "GLU", 1.625f, 0.0f, 0.0f, 1.0f,  1,  8 },
         {" CA ", "GLU", 1.750f, 0.3f, 0.3f, 0.3f, -1,  8 },
         {" C  ", "GLU", 1.875f, 0.3f, 0.3f, 0.3f,  1,  8 },
         {" O  ", "GLU", 1.480f, 1.0f, 0.0f, 0.0f,  1,  8 },
         {" CB ", "GLU", 1.750f, 0.3f, 0.3f, 0.3f, -1,  8 },
         {" CG ", "GLU", 1.750f, 0.3f, 0.3f, 0.3f, -1,  8 },
         {" CD ", "GLU", 1.875f, 0.3f, 0.3f, 0.3f,  1,  8 },
         {" OE1", "GLU", 1.480f, 1.0f, 0.0f, 0.0f,  1,  8 },
         {" OE2", "GLU", 1.480f, 1.0f, 0.0f, 0.0f,  1,  8 }
    };

    static std::map<unsigned long, std::map<unsigned long, PDBelementInformation>> PDBelementMap;
}

static unsigned long str2long(std::string s)
{
    unsigned long rval = 0;

    for(char c : s)
    {
        rval |= c;
        rval = rval << 8;
    }

    return rval;
}

/*
 * ***************************************************************************
 * Routine:  getMinMax    < ... >
 *
 * Authors:   Zeyun Yu (zeyun.yu@gmail.com)
 *            John Moody (brogan@gmail.com)
 *
 * Purpose:  Calculate the minimum and maximum range of the blurred atoms
 * ***************************************************************************
 */
template <typename Iterator>
void getMinMax(Iterator begin,
               Iterator end,
               float min[3],
               float max[3])
{
    float maxRad = 0.0;
    float tempRad;

    min[0] = min[1] = min[2] = std::numeric_limits<float>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<float>::infinity();
    maxRad = 0;

    for (auto curr = begin; curr != end; ++curr)
    {
        float x = curr->x;
        float y = curr->y;
        float z = curr->z;

        min[0] = std::min(min[0], x);
        max[0] = std::max(max[0], x);

        min[1] = std::min(min[1], y);
        max[1] = std::max(max[1], y);

        min[2] = std::min(min[2], z);
        max[2] = std::max(max[2], z);

        tempRad = curr->radius * sqrt(1.0 + log(detail::EPSILON) / BLOBBYNESS);

        maxRad = std::max(maxRad, tempRad);
    }

    for (int i = 0; i < 3; i++)
    {
        min[i] -= maxRad;
        max[i] += maxRad;
    }
}


/*
 * ***************************************************************************
 * Routine:  readPDB    < ... >
 *
 * Author:   John Moody (brogan@gmail.com)
 *
 * Purpose:  Read a PDB file and extract the Atom records.
 *
 * Notes:    Requires strict adherence to the PDB file format.
 * ***************************************************************************
 */
template <typename Inserter>
void readPDB(std::string filename, Inserter inserter)
{
    std::ifstream infile(filename);
    std::string line;

    if(detail::PDBelementMap.size() == 0)
    {
        for(int i = 0; i < detail::MAX_BIOCHEM_ELEMENTS; ++i)
        {
            const char* atomName    = detail::PDBelementTable[i].atomName;
            const char* residueName = detail::PDBelementTable[i].residueName;

            detail::PDBelementMap[str2long(atomName)][str2long(residueName)] = detail::PDBelementTable[i];
        }
    }

    if(infile.is_open())
    {
        while (std::getline(infile, line))
        {
            std::smatch match;
            if (std::regex_match(line, match, detail::atom))
            {
                ATOM new_atom;
                // see Format_v33_Letter.pdf in gamer/doc
                new_atom.x = std::atof(line.substr(30,8).c_str());
                new_atom.y = std::atof(line.substr(38,8).c_str());
                new_atom.z = std::atof(line.substr(46,8).c_str());

                new_atom.radius = 1.0f; // default radius
                std::string atomName = line.substr(12,4);
                std::string residueName = line.substr(17,3);

                auto eInfo = detail::PDBelementMap[str2long(line.substr(12,4))][str2long(line.substr(17,3))];
                new_atom.radius = eInfo.radius;

                *inserter++ = new_atom;
            }
        }
    }
    else
    {
        std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
    }
}


/*
 * ***************************************************************************
 * Routine:  readPQR    < ... >
 *
 * Author:   John Moody (brogan@gmail.com)
 *
 * Purpose:  Read a PDB file and extract the Atom records.
 *
 * Notes:    Requires strict adherence to the PDB file format.
 * ***************************************************************************
 */
template <typename Inserter>
void readPQR(std::string filename, Inserter inserter)
{
    char   line[256];
    char   string[8];
    int   k, m, n;
    FILE  *fp, *fout;
    char  file_name[256];
    sprintf(file_name, "%s.xyzr", filename.c_str());

    if ((fout = fopen(file_name, "wb")) == NULL)
    {
        printf("write error...\n");
        exit(0);
    }

    if ((fp = fopen(filename.c_str(), "r")) == NULL)
    {
        printf("read error...\n");
        exit(0);
    }

    m = 0;
    while (fgets(line, 256, fp) != NULL)
    {
        if ((line[0] == 'A') && (line[1] == 'T') && (line[2] == 'O') && (line[3] == 'M'))
        {
            /* more general format, could be used for pqr format */
            k = 30;

            ATOM new_atom;

            while (line[k] == ' ')
            {
                k++;
            }
            n = 0;

            while (line[k] != ' ')
            {
                string[n] = line[k];
                n++;
                k++;

                if (line[k] == '-')
                {
                    break;
                }
            }
            string[n]      = '\0';
            new_atom.x = atof(string);

            while (line[k] == ' ')
            {
                k++;
            }
            n = 0;

            while (line[k] != ' ')
            {
                string[n] = line[k];
                n++;
                k++;

                if (line[k] == '-')
                {
                    break;
                }
            }
            string[n]      = '\0';
            new_atom.y = atof(string);

            while (line[k] == ' ')
            {
                k++;
            }
            n = 0;

            while (line[k] != ' ')
            {
                string[n] = line[k];
                n++;
                k++;

                if (line[k] == '-')
                {
                    break;
                }
            }
            string[n]      = '\0';
            new_atom.z = atof(string);

            // skip whitespace
            while (line[k] == ' ')
            {
                k++;
            }
            n = 0;

            // read field
            while (line[k] != ' ')
            {
                string[n] = line[k];
                n++;
                k++;

                // This makes no sense
                if (line[k] == '-')
                {
                    break;
                }
            }

            // skip whitespace
            while (line[k] == ' ')
            {
                k++;
            }
            n = 0;

            while (line[k] != ' ' && line[k] != '\n' && line[k] != '\0')
            {
                string[n] = line[k];
                n++;
                k++;
            }
            string[n]           = '\0';
            new_atom.radius = atof(string);

            // new_atom.radius += 1.0;

            if (new_atom.radius < 1.0)
            {
                new_atom.radius = 1.0;
            }

            *inserter++ = new_atom;
            fprintf(fout, "%f %f %f %f\n", new_atom.x, new_atom.y, new_atom.z, new_atom.radius);
            m++;
        }
    }
    fclose(fp);
    fclose(fout);
}


/*
 * ***************************************************************************
 * Routine:  readXYZR    < ... >
 *
 * Author:   John Moody (brogan@gmail.com)
 *
 * Purpose:  Read a PDB file and extract the Atom records.
 *
 * Notes:    Requires strict adherence to the PDB file format.
 * ***************************************************************************
 */
template <typename Inserter>
void readXYZR(std::string filename, Inserter inserter)
{
    float x, y, z, radius;
    ATOM new_atom;
    FILE *fp;
    int   m;

    if ((fp = fopen(filename.c_str(), "r")) == NULL)
    {
        printf("read error...\n");
        exit(0);
    }

    fscanf(fp, "%d\n", &m);

    for (int n = 0; n < m; n++)
    {
        fscanf(fp, "%f %f %f %f\n", &x, &y, &z, &radius);

        new_atom.x      = x;
        new_atom.y      = y;
        new_atom.z      = z;
        new_atom.radius = radius;

        *inserter++ = new_atom;
    }

    printf("number of atoms: %d\n", m);
}


/*
 * ***************************************************************************
 * Routine:  readAtomFile    < ... >
 *
 * Author:   John Moody (brogan@gmail.com)
 *
 * Purpose:  Read a PDB file and extract the Atom records.
 *
 * Notes:    Requires strict adherence to the PDB file format.
 * ***************************************************************************
 */
template <typename Inserter>
void readAtomFile(std::string filename, Inserter inserter)
{
    // XYZR file
    if (std::regex_match(filename, detail::XYZR))
    {
        readXYZR(filename, std::move(inserter));
    }
    // PDB file
    if (std::regex_match(filename, detail::PDB))
    {
        readPDB(filename, std::move(inserter));
    }
    // PQR file
    else if (std::regex_match(filename, detail::PQR))
    {
        readPQR(filename, std::move(inserter));
    }
    else
    {
        printf("Input file name end with PDB/PQR/XYZR.\n");
        exit(0);
    }
}
