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

cmake_minimum_required(VERSION 3.12)
# Disable in source builds
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

if(SKBUILD)
    message(STATUS "The project is being built using scikit-build")
endif()

# Look for version from GIT
include(GetGitRevisionDescription)
git_describe(VERSION --tags --dirty=-d)

#parse the version information into pieces.
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(-alpha|-beta|)-(.*)" "\\2" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

#####################################################################
# Project GAMer
#####################################################################
project(GAMer VERSION ${VERSION_SHORT})

message(STATUS "GAMer version: ${VERSION_SHORT}")

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules/version.cpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")
list(APPEND GAMER_SOURCES "${version_file}")

#####################################################################
# Options
#####################################################################
option(SINGLE "Use single precision floating point numbers?" OFF)

option(BUILD_PYGAMER "Build GAMer python extension?" OFF)
option(BUILD_BLENDER "Build the GAMer addon for Blender?" OFF)
option(BLENDER_PLUGIN_INSTALL "Have CMake install the Blender plugin?" OFF)
option(BLENDER_VERSION_STRICT "Have CMake verify compatibility of plugin with Blender?" OFF)

option(GAMER_TESTS "Build the GAMer tests?" OFF)
option(GETEIGEN "Download Eigen?" ON)

option(GAMER_CMAKE_VERBOSE "Print out information for debugging CMake configuration?")

#####################################################################
# Configuration
#####################################################################
# Require c++14 and standard libraries
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Define where to put the libraries and binaries
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Override rules for MSVC static compiler flags
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_LIST_DIR}/cmake-modules/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_LIST_DIR}/cmake-modules/cxx_flag_overrides.cmake)

# # Add -fPIC to all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set build type to Release if none was defined
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type was specified. Setting to 'Release'")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Required to tell MSVC to export symbols
if(MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Blender Addon requires Python extension
if(BUILD_BLENDER)
    set(BUILD_PYGAMER ON)
endif()

## Print compile flags
if(GAMER_CMAKE_VERBOSE)
    message(STATUS "CMAKE_C_FLAGS is: ${CMAKE_C_FLAGS}")
    message(STATUS "CMAKE_C_FLAGS_DEBUG is: ${CMAKE_C_FLAGS_DEBUG}")
    message(STATUS "CMAKE_C_FLAGS_RELEASE is: ${CMAKE_C_FLAGS_RELEASE}")
    message(STATUS "CMAKE_C_FLAGS_RELWITHDEBINFO is: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
    message(STATUS "CMAKE_C_FLAGS_MINSIZEREL is: ${CMAKE_C_FLAGS_MINSIZEREL}")
    message(STATUS "CMAKE_CXX_FLAGS is: ${CMAKE_CXX_FLAGS}")
    message(STATUS "CMAKE_CXX_FLAGS_DEBUG is: ${CMAKE_CXX_FLAGS_DEBUG}")
    message(STATUS "CMAKE_CXX_FLAGS_RELEASE is: ${CMAKE_CXX_FLAGS_RELEASE}")
    message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO is: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
    message(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL is: ${CMAKE_CXX_FLAGS_MINSIZEREL}")
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE})
endif(GAMER_CMAKE_VERBOSE)

# Add and configure library dependencies
add_subdirectory(libraries EXCLUDE_FROM_ALL)

list(APPEND GAMER_SOURCES
    "src/OFF_SurfaceMesh.cpp"
    "src/OBJ_SurfaceMesh.cpp"
    "src/SurfaceMesh.cpp"
    "src/SurfaceMeshDetail.cpp"
    "src/Vertex.cpp"
    "src/TetMesh.cpp"
    "src/PDBReader.cpp"
    "src/pdb2mesh.cpp"
)

#####################################################################
# LIBRARIES
#####################################################################
# OBJECT LIBRARY: compiles the sources only once
add_library(gamer_objlib OBJECT ${GAMER_SOURCES})
target_include_directories(gamer_objlib PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
target_link_libraries(gamer_objlib PUBLIC casc tetstatic eigen)
# target_compile_options(gamer_objlib PRIVATE -Werror -Wall -Weverything -Wextra -pedantic-errors -Wconversion -Wsign-conversion -ferror-limit=3 -Wno-c++17-extensions)
# -ferror-limit=3 -Wno-unused-variable -Wno-unused-parameter)

# SHARED LIBRARY
add_library(gamershared SHARED $<TARGET_OBJECTS:gamer_objlib>)
target_link_libraries(gamershared PUBLIC gamer_objlib)
set_target_properties(gamershared PROPERTIES OUTPUT_NAME gamer)

# STATIC LIBRARY
add_library(gamerstatic STATIC $<TARGET_OBJECTS:gamer_objlib>)
target_link_libraries(gamerstatic PUBLIC gamer_objlib)
if(NOT WIN32)
    # Shared and static libs will clobber each other on Windows
    set_target_properties(gamerstatic PROPERTIES OUTPUT_NAME gamer)
endif()

include(GNUInstallDirs)
# Install targets and headers
install(TARGETS gamershared gamerstatic
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/gamer DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")


#####################################################################
# PYTHON EXTENSION AND BLENDER ADDON
#####################################################################
if(BUILD_PYGAMER)
    add_subdirectory(pygamer)
endif()

if(BUILD_BLENDER)
    add_subdirectory(tools)
endif()

#####################################################################
# TESTING AND DOCUMENTATION
#####################################################################
if(GAMER_TEST)
    enable_testing()
    add_subdirectory(tests)
endif()

# Configure documentation builders
add_subdirectory(docs)
