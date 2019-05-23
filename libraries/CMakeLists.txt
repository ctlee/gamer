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

###########################
# GET CASC
###########################
FetchContent_Declare(
    casclib
    GIT_REPOSITORY  https://github.com/ctlee/casc.git
    GIT_TAG         origin/development
    SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/casc-src"
    BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/casc-build"
)
FetchContent_GetProperties(casclib)
if(NOT casclib_POPULATED)
  FetchContent_Populate(casclib)
  add_subdirectory(${casclib_SOURCE_DIR} ${casclib_BINARY_DIR})
endif()

###########################
# GET PYBIND11 2.2.4
###########################
FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/v2.2.4.tar.gz
    SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/pybind11-src"
    BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/pybind11-build"
)
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
  FetchContent_Populate(pybind11)
  add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

if(NOT GETEIGEN)
    find_package (Eigen3 3.3 NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(FATAL_ERROR "Could not find required library Eigen. Please append -DGETEIGEN=true to your cmake call and I will download Eigen for you.")
    endif()
    add_library(eigen ALIAS Eigen3::Eigen)
endif()

if(GETEIGEN)
    ###########################
    # GET EIGEN 3.3.7
    ###########################
    FetchContent_Declare(
        eigen
        URL https://github.com/eigenteam/eigen-git-mirror/archive/3.3.7.tar.gz
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-build"
    )
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
        add_library (eigen INTERFACE)
        target_include_directories(eigen INTERFACE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)
    endif()
endif()

add_subdirectory(tetgen)
# add_subdirectory(triangle)
