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
 * Lesser General Public License for more pdbreader_details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * ***************************************************************************
 */

/*
 * Parts of this are adapted from the PDBParser developed in Chandrajit
 * Bajaj's group at The University of Texas at Austin for GAMer by Zeyun Yu.
 */

#pragma once

#include <algorithm>
#include <regex>
#include <iostream>
#include <string>
#include <fstream>
#include <array>

#include "gamer/gamer.h"
#include "gamer/Vertex.h"


/// Namespace for all things gamer
namespace gamer
{

/// @cond detail
namespace pdbreader_detail
{
const double            EPSILON = 1e-3;
static const std::regex PDB(".*.pdb", std::regex::icase | std::regex::optimize);
static const std::regex PQR(".*.pqr", std::regex::icase | std::regex::optimize);
static const std::regex XYZR(".*.xyzr", std::regex::icase | std::regex::optimize);
static const std::regex atom("ATOM.*\n*", std::regex::optimize);

/**
 * @brief      PDB element information
 */
struct PDBelementInformation
{
    const char  * atomName;
    const char  * residueName;
    float         radius;
    float         red;
    float         green;
    float         blue;
    int           hydrophobicity;
    unsigned char residueIndex;
};

/// Total number of elements in the table
static const std::size_t     MAX_BIOCHEM_ELEMENTS = 167;

/// Lookup table
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
    // {"SI  ", "UNL", 1.875f, 1.0f, 1.0f, 1.0f,  1, 27 },
    // {" O  ", "UNL", 1.480f, 1.0f, 0.0f, 0.0f,  1, 27 }
};

/// Map of residueName, atomName to PDBelementInformation
static std::map<std::string, std::map<std::string, PDBelementInformation> > PDBelementMap;
} // End namespace pdbreader_detail
/// @endcond


/**
 * @brief      Basic atom containing position and radius
 */
struct Atom {
    f3Vector pos;    /**< @brief position */
    double   radius; /**< @brief radius */
};

/**
 * @brief      Initialize the element map
 */
static void initElementMap()
{
    // If the element map is empty build it...
    if (pdbreader_detail::PDBelementMap.size() == 0)
    {
        for (int i = 0; i < pdbreader_detail::MAX_BIOCHEM_ELEMENTS; ++i)
        {
            const std::string atomName = pdbreader_detail::PDBelementTable[i].atomName;
            const std::string residueName = pdbreader_detail::PDBelementTable[i].residueName;
            pdbreader_detail::PDBelementMap[residueName][atomName] = pdbreader_detail::PDBelementTable[i];
        }
    }
}

/**
 * @brief      Read a PDB file
 *
 * @param[in]  filename  The filename
 * @param[in]  inserter  The inserter
 *
 * @tparam     Inserter  Inserter
 *
 * @return     True if successfully read
 */
template <typename Inserter>
bool readPDB(const std::string &filename, Inserter inserter)
{
    std::ifstream infile(filename);
    std::string   line;

    initElementMap(); // Init the element map if it

    if (infile.is_open())
    {
        while (std::getline(infile, line))
        {
            std::smatch match;
            if (std::regex_match(line, match, pdbreader_detail::atom))
            {
                Atom  atom;
                // See PDB file formatting guidelines
                float x = std::atof(line.substr(30, 8).c_str());
                float y = std::atof(line.substr(38, 8).c_str());
                float z = std::atof(line.substr(46, 8).c_str());
                atom.pos = Vector({x, y, z});

                atom.radius = 1.0f; // default radius
                std::string atomName = line.substr(12, 4);
                std::string residueName = line.substr(17, 3);

                auto        innerMapIT = pdbreader_detail::PDBelementMap.find(residueName);

                if (innerMapIT != pdbreader_detail::PDBelementMap.end())
                {
                    auto innerMap = innerMapIT->second;
                    auto typeIT   = innerMap.find(atomName);
                    if (typeIT != innerMap.end())
                    {
                        atom.radius = typeIT->second.radius;
                    }
                    else
                    {
                        std::cout << "Could not find atomtype of '"
                                  << atomName << "' in residue '"
                                  << residueName << "'. "
                                  << "Using default radius." << std::endl;
                    }
                }
                else
                {
                    std::cout << "Could not find ResidueName '"
                              << residueName << "' in table. "
                              << "Using default radius." << std::endl;
                }
                *inserter++ = atom;
            }
        }
        return true;
    }
    else
    {
        std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
        return false;
    }
}

/**
 * @brief      Reads a pqr.
 *
 * @param[in]  filename  The filename
 * @param[in]  inserter  The inserter
 *
 * @tparam     Inserter  { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename Inserter>
bool readPQR(const std::string &filename, Inserter inserter)
{
    std::ifstream infile(filename);
    std::string   line;

    initElementMap(); // Init the element map if it

    if (infile.is_open())
    {
        while (std::getline(infile, line))
        {
            std::smatch match;
            if (std::regex_match(line, match, pdbreader_detail::atom))
            {
                Atom  atom;
                // See PDB file formatting guidelines
                float x = std::atof(line.substr(30, 8).c_str());
                float y = std::atof(line.substr(38, 8).c_str());
                float z = std::atof(line.substr(46, 8).c_str());
                atom.pos = Vector({x, y, z});
                atom.radius = std::atof(line.substr(62, 7).c_str());
                *inserter++ = atom;
            }
        }
        return true;
    }
    else
    {
        std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
        return false;
    }
}

template <typename Iterator, typename BlurFunc>
void getMinMax(Iterator begin, Iterator end, f3Vector &min, f3Vector &max, BlurFunc &&f)
{
    float maxRad = 0.0;
    float tmpRad;

    min[0] = min[1] = min[2] = std::numeric_limits<float>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<float>::infinity();

    for (auto curr = begin; curr != end; ++curr)
    {
        float x = curr->pos[0];
        float y = curr->pos[1];
        float z = curr->pos[2];

        if (min[0] > x)
            min[0] = x;
        if (max[0] < x)
            max[0] = x;

        if (min[1] > y)
            min[1] = y;
        if (max[1] < y)
            max[1] = y;

        if (min[2] > z)
            min[2] = z;
        if (max[2] < z)
            max[2] = z;

        tmpRad = f(curr->radius); // * sqrt(1.0 + log(pdbreader_detail::EPSILON)
                                  // / blobbyness);
        if (maxRad < tmpRad)
            maxRad = tmpRad;
    }

    min -= f3Vector({maxRad, maxRad, maxRad});
    max += f3Vector({maxRad, maxRad, maxRad});
}

/**
 * @brief      Apply a gaussian blur to a list of atoms
 *
 * @param[in]  begin       The begin
 * @param[in]  end         The end
 * @param      dataset     The dataset
 * @param[in]  min         The minimum
 * @param[in]  maxMin      The maximum minimum
 * @param[in]  dim         The dim
 * @param[in]  blobbyness  The blobbyness
 *
 * @tparam     Iterator    Typename of the iterator
 */
template <typename Iterator>
void blurAtoms(Iterator begin, Iterator end,
               float* dataset,
               const f3Vector &min,
               const f3Vector &maxMin,
               const i3Vector &dim,
               float blobbyness)
{

    // Functor to calculate gaussian blur
    auto evalDensity = [blobbyness](const Atom &atom, f3Vector &pnt, float maxRadius)
        -> float {
            double   expval;

            f3Vector tmp = atom.pos - pnt;
            double   r   = tmp|tmp;
            double   r0  = atom.radius*atom.radius;

            // expval = BLOBBYNESS*(r/r0 - 1.0);
            expval = blobbyness*(r-r0);

            // Truncate gaussian
            if (sqrt(r) > maxRadius)
            {
                return 0.0;
            }
            return (float) exp(expval);
        };

    f3Vector span;
    span = (maxMin).ElementwiseDivision(static_cast<f3Vector>((dim - i3Vector({1, 1, 1}))));

    float radFactor = sqrt(1.0 + log(pdbreader_detail::EPSILON)/(2.0 * blobbyness));

    for (auto curr = begin; curr != end; ++curr)
    {
        float    maxRad = curr->radius * radFactor;
        // compute the dataset coordinates of the atom's center
        f3Vector tmpVec = (curr->pos-min).ElementwiseDivision(span);
        i3Vector c;
        std::transform(tmpVec.begin(), tmpVec.end(), c.begin(), [](float v) -> int {
                return round(v);
            });

        // std::cout << "Max Radius: " << maxRad << std::endl;

        // compute the bounding box of the atom (maxRad^3)
        i3Vector amin;
        i3Vector amax;
        for (int j = 0; j < 3; ++j)
        {
            int   tmp;
            float tmpRad = maxRad/span[j];

            tmp = (int)(c[j] - tmpRad - 1);
            amin[j] = (tmp < 0) ? 0 : tmp; // check if tmp is < 0
            tmp = (int)(c[j] + tmpRad + 1);
            amax[j] = (tmp > (dim[j] - 1)) ? (dim[j] - 1) : tmp;
        }

        // std::cout << amin << " " << amax << std::endl;

        // Blur kernel in bounding box
        for (int k = amin[2]; k <= amax[2]; k++)
        {
            for (int j = amin[1]; j <= amax[1]; j++)
            {
                for (int i = amin[0]; i <= amax[0]; i++)
                {
                    f3Vector pnt = min + f3Vector({static_cast<float>(i),
                                                   static_cast<float>(j),
                                                   static_cast<float>(k)}).ElementwiseProduct(span);
                    dataset[Vect2Index(i, j, k, dim)] += evalDensity(*curr, pnt, maxRad);
                }
            }
        }
    }
}

/**
 * @brief      Compute the grid based Solvent Accessible Area.
 *
 * Assumes that the domain is in the {+,+,+} octant.
 *
 * @param[in]  begin       The begin
 * @param[in]  end         The end
 * @param[in]  dim         The dim
 * @param      dataset     The atom index
 *
 * @tparam     Iterator    { description }
 */
template <typename Iterator>
void gridSAS(const Iterator begin, const Iterator end, const i3Vector &dim, float* dataset)
{
    // For atom in atoms :
    for (auto curr = begin; curr != end; ++curr)
    {
        float    radius = curr->radius;
        f3Vector pos = curr->pos;
        // compute the dataset coordinates of the atom's center
        i3Vector c;
        std::transform(pos.begin(), pos.end(), c.begin(), [](float v) -> int {
                return round(v);
            });

        // compute bounding box for atom
        i3Vector amin;
        i3Vector amax;
        for (int j = 0; j < 3; ++j)
        {
            int tmp;
            tmp = (int)(c[j] - radius - 1);
            amin[j] = (tmp < 0) ? 0 : tmp; // check if tmp is < 0
            tmp = (int)(c[j] + radius + 1);
            amax[j] = (tmp > (dim[j] - 1)) ? (dim[j] - 1) : tmp;
        }

        // Blur kernel in bounding box
        for (int k = amin[2]; k <= amax[2]; k++)
        {
            for (int j = amin[1]; j <= amax[1]; j++)
            {
                for (int i = amin[0]; i <= amax[0]; i++)
                {
                    f3Vector coord = f3Vector({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
                    coord -= curr->pos;
                    float    dist = -(std::sqrt(coord|coord)-radius); // inside
                                                                      // is
                                                                      // positive
                    int      idx  = Vect2Index(i, j, k, dim);

                    if (dist > dataset[idx])
                    {
                        dataset[idx] = dist;
                    }
                }
            }
        }
    }
}

template <typename Iterator>
void gridSES(const Iterator begin, const Iterator end, const i3Vector &dim,
             float* dataset, const float radius)
{
    for (auto curr = begin; curr != end; ++curr)
    {
        f3Vector pos = (*curr).position;

        // compute bounding box for atom
        i3Vector amin;
        i3Vector amax;
        for (int i = 0; i < 3; ++i)
        {
            int tmp;
            tmp = (int)(pos[i] - radius - 1);
            amin[i] = (tmp < 0) ? 0 : tmp; // check if tmp is < 0
            tmp = (int)(pos[i] + radius + 1);
            amax[i] = (tmp > (dim[i] - 1)) ? (dim[i] - 1) : tmp;
        }

        // Blur kernel in bounding box
        for (int k = amin[2]; k <= amax[2]; k++)
        {
            for (int j = amin[1]; j <= amax[1]; j++)
            {
                for (int i = amin[0]; i <= amax[0]; i++)
                {
                    f3Vector coord = f3Vector({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
                    coord -= pos;
                    float    dist = -(std::sqrt(coord|coord)-radius);
                    int      idx  = Vect2Index(i, j, k, dim);
                    if (dist > dataset[idx])
                    {
                        dataset[idx] = dist;
                    }
                }
            }
        }
    }
}

std::unique_ptr<SurfaceMesh> readPDB_molsurf(const std::string &filename);
std::unique_ptr<SurfaceMesh> readPDB_gauss(const std::string &filename, float blobbyness, float isovalue);
std::unique_ptr<SurfaceMesh> readPDB_distgrid(const std::string &filename, const float radius);

std::unique_ptr<SurfaceMesh> readPQR_gauss(const std::string &filename, float blobbyness, float isovalue);

} // end namespace gamer
