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
# GET TUTORIALS
###########################
FetchContent_Declare(
    gamer_tutorials
    GIT_REPOSITORY  https://github.com/ctlee/gamer_tutorials.git
    GIT_TAG         origin/master
    SOURCE_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/src/tutorials"
)
FetchContent_GetProperties(gamer_tutorials)
if(NOT gamer_tutorials_POPULATED)
    FetchContent_Populate(gamer_tutorials)
endif()

#####################################################################
# Target to generate Doxygen documentation
#####################################################################
find_package(Doxygen)
if(DOXYGEN_FOUND)
    # message(STATUS "Doxygen Executable: ${DOXYGEN_EXECUTABLE}")
    message(STATUS "Build the Doxygen GAMer documentation using `make doxy_docs`")
    # Parse Doxyfile.in and replace CMake @macros@.
    configure_file(${CMAKE_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doxy_docs
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
else()
    message(STATUS "CMake could not find Doxygen. Please install Doxygen or help me find your Doxygen binary to build the documentation!")
endif(DOXYGEN_FOUND)

find_package(Sphinx)
find_package(Exhale)
find_package(Breathe)

SET(SPHINXOPTS CACHE STRING "Extra options to pass to Sphinx")

if(SPHINX_FOUND)
    set(DOXYGEN_TO_SPHINX False)
    if(DOXYGEN_FOUND AND EXHALE_FOUND AND BREATHE_FOUND)
        set(DOXYGEN_TO_SPHINX True)
    endif()

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
        @ONLY)

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/conf.py"
        @ONLY)

    add_custom_target(sphinx_docs
        ${SPHINX_EXECUTABLE}
            -M html
            "${CMAKE_CURRENT_SOURCE_DIR}/src"
            "${CMAKE_CURRENT_BINARY_DIR}"
            -c "${CMAKE_CURRENT_BINARY_DIR}"
            ${SPHINXOPTS}
        DEPENDS pygamer
        COMMENT "Building HTML documentation with Sphinx")
    message(STATUS "Build the Sphinx documentation using `make sphinx_docs`")
else()
    message(STATUS "CMake could not find Sphinx. Please install Sphinx or help me find your Sphinx binary to build the documentation!")
endif(SPHINX_FOUND)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/_cppapi"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/_pythonapi"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/_doxyoutput"
        "${CMAKE_CURRENT_BINARY_DIR}/doctrees"
        "${CMAKE_CURRENT_BINARY_DIR}/html")