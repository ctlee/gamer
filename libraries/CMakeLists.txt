# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2018
# by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
#    and Michael Holst

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# ***************************************************************************

include(FetchContent)

FetchContent_Declare(
    casc
    GIT_REPOSITORY git@github.com:ctlee/casc.git
    GIT_TAG        origin/development
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/casc
)
FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/v2.2.4.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pybind11
)

FetchContent_GetProperties(casc)
if(NOT casc_POPULATED)
  FetchContent_Populate(casc)
endif()

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
  FetchContent_Populate(pybind11)
  add_subdirectory(${pybind11_SOURCE_DIR})
endif()

add_subdirectory(tetgen)
# add_subdirectory(triangle)