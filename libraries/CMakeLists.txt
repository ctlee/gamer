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


find_package(Git)
if(NOT Git_FOUND)
  message(FATAL_ERROR "Git could not be found. The install process requires Git to download dependencies.")
endif()

include(FetchContent)

###########################
# GET CASC
###########################
FetchContent_Declare(
    casclib
    GIT_REPOSITORY  https://github.com/ctlee/casc.git
    GIT_TAG         v1.0.5
    GIT_SHALLOW     TRUE 
    SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/casc-src"
    BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/casc-build"
)
FetchContent_GetProperties(casclib)
if(NOT casclib_POPULATED)
  FetchContent_Populate(casclib)
  add_subdirectory(${casclib_SOURCE_DIR} ${casclib_BINARY_DIR})
endif()

###########################
# FIND/GET EIGEN
###########################
if(NOT GETEIGEN)
    find_package (Eigen3 3.3 NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(FATAL_ERROR "Could not find required library Eigen."
            "Please append -DGETEIGEN=on to your cmake call and I will download Eigen for you.")
    endif()
    target_compile_options(Eigen3::Eigen PRIVATE -w)
else()
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY  https://gitlab.com/libeigen/eigen.git
        GIT_TAG         3.3.9
        GIT_SHALLOW     TRUE 
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-build"
    )
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
        add_library (eigen INTERFACE)
        add_library (Eigen3::Eigen ALIAS eigen)
        target_include_directories(eigen INTERFACE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)
        target_compile_options(eigen INTERFACE -w)   # Suppress warnings
    endif()
endif()

###########################
# Tetgen
###########################
add_subdirectory(tetgen EXCLUDE_FROM_ALL)
# add_subdirectory(triangle EXCLUDE_FROM_ALL)
