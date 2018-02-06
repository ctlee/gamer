/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2017
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

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <ostream>
#include <regex>
#include <cmath>
#include <vector>

#include "SurfaceMesh.h"
#include "PDBReader.h"


std::unique_ptr<SurfaceMesh> readPDB_gauss(const std::string& filename, 
 			double blobbyness, 
			float iso_value){
	std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

	std::ifstream fin(filename);

	std::vector<AtomType> atomTypes;
    if(!readPDB(filename, std::back_inserter(atomTypes)))
    {
        return mesh;
    }

    std::cout << "Begin blurring coordinates..." << std::endl;

    std::cout << "Done blurring coords" << std::endl;
}