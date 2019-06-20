/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2019
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gamer/SurfaceMesh.h"

/// Namespace for all things gamer
namespace gamer
{

namespace py = pybind11;

void init_SMFunctions(py::module& mod){
    mod.def("cube", &cube,
        py::arg("order"),
        R"delim(
            Construct a triangulated cube mesh

            Args:
                order (int): Number of subdivisions

            Returns:
                :py:class:`SurfaceMesh`: Mesh cube
        )delim"
    );


    mod.def("sphere", &sphere,
        py::arg("order"),
        R"delim(
            Construct a triangulated icosahedral sphere mesh

            Args:
                order (int): Number of subdivisions

            Returns:
                :py:class:`SurfaceMesh`: Spherical mesh
        )delim"
    );
}

} // end namespace gamer