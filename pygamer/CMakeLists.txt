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
# GET PYBIND11
###########################
if(NOT GETPYBIND11)
    find_package(pybind11 CONFIG REQUIRED)
    if(NOT pybind11_FOUND)
        message(FATAL_ERROR "Could not find required library pybind11."
            "Please append -DGETPYBIND11=ON to your cmake call and I will download pybind11 for you.")
    endif()
else()
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY  https://github.com/pybind/pybind11.git
        GIT_TAG         v2.7.1
        GIT_SHALLOW     TRUE 
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/pybind11-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/pybind11-build"
    )
    FetchContent_GetProperties(pybind11)
    if(NOT pybind11_POPULATED)
      FetchContent_Populate(pybind11)
      add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
endif()

list(APPEND PYGAMER_SOURCES
    "src/Vector.cpp"

    "src/SMGlobal.cpp"
    "src/SMVertex.cpp"
    "src/SMEdge.cpp"
    "src/SMFace.cpp"
    "src/SMSimplexID.cpp"
    "src/SMFunctions.cpp"
    "src/SurfaceMesh.cpp"

    "src/TMGlobal.cpp"
    "src/TMVertex.cpp"
    "src/TMEdge.cpp"
    "src/TMFace.cpp"
    "src/TMCell.cpp"
    "src/TMSimplexID.cpp"
    "src/TetMesh.cpp"

    "src/pygamer.cpp"
    )

# Create the target
pybind11_add_module(pygamer MODULE ${PYGAMER_SOURCES})
# Link pygamer to gamer static
target_link_libraries(pygamer PRIVATE gamerstatic)

file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
             INPUT ${CMAKE_SOURCE_DIR}/cmake/__init__.py.in)

if(SKBUILD)
    install(TARGETS pygamer LIBRARY DESTINATION pygamer)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION pygamer)
endif()
