# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2021
# by Christopher T. Lee and contributors

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see <http://www.gnu.org/licenses/>
# or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# ***************************************************************************

# 3.12 required for FetchContent and objectlib INTERFACE include inheritance
cmake_minimum_required(VERSION 3.12...3.18)

# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(SKBUILD)
    message(STATUS "The project is being built using scikit-build")
    # Currently, Scikit-build does not support FindPython, so we convert the
    # provided hints ourselves.
    set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
    set(Python_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
    set(Python_LIBRARY "${PYTHON_LIBRARY}")
    set(Python_VERSION "${PYTHON_VERSION_STRING}")
endif()

include(Macros)
set_default_build_type()
get_version_from_git()

#####################################################################
# Project GAMer
#####################################################################
project(GAMer VERSION ${VERSION_SHORT})
message(STATUS "GAMer version: ${VERSION}")

#####################################################################
# Options
#####################################################################
option(SINGLE "Use single precision floating point numbers?" OFF)

option(BUILD_PYGAMER "Build GAMer python extension?" OFF)
option(BUILD_BLENDGAMER "Build the GAMer addon for Blender?" OFF)
option(BLENDER_PLUGIN_INSTALL "Have CMake install the Blender plugin?" OFF)
option(BLENDER_VERSION_STRICT "Have CMake verify compatibility of plugin with Blender?" OFF)

option(GAMER_TESTS "Build the GAMer tests?" OFF)
option(GETEIGEN "Download Eigen?" ON)
option(GETPYBIND11 "Download pybind11?" ON)

option(GAMER_DOCS "Download and configure documentation?" OFF)

option(VECTORIZE "Enable vectorization?" OFF)

option(BLENDER_VERSION_OVERRIDE "Override the version number" "")
mark_as_advanced(BLENDER_VERSION_OVERRIDE)
option(GAMER_CMAKE_VERBOSE "Print out information for debugging CMake configuration?" OFF)
mark_as_advanced(GAMER_CMAKE_VERBOSE)

option(FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." OFF)
if (FORCE_COLORED_OUTPUT)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

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
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_LIST_DIR}/cmake/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_LIST_DIR}/cmake/cxx_flag_overrides.cmake)

if(VECTORIZE)
    message(WARNING "Vectorization in GAMer is experimental, use at your own risk. Enabling vectorization may make your libraries less portable across machines.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# Compile GAMer in memory where possible. Note that this can break compilation
# on remote VMs such at readthedocs.io
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)
cmake_policy(SET CMP0063 NEW)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# # Add -fPIC to all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Required to tell MSVC to export symbols
if(MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Blender Addon requires Python extension
if(BUILD_BLENDGAMER)
    set(BUILD_PYGAMER ON)
endif()

## Print compile flags
if(GAMER_CMAKE_VERBOSE)
    print_debug_messages()
endif(GAMER_CMAKE_VERBOSE)

if(BUILD_PYGAMER)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

# Add and configure library dependencies
add_subdirectory(libraries EXCLUDE_FROM_ALL)
add_subdirectory(include)
# Src must be traversed prior to appending to GAMER_SOURCES
add_subdirectory(src)
# Configure version.cpp to give access to version in code
configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
list(APPEND GAMER_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

#####################################################################
# LIBRARIES
#####################################################################
# OBJECT LIBRARY: compiles the sources only once
add_library(gamer_objlib OBJECT ${GAMER_INCLUDES} ${GAMER_SOURCES})
target_include_directories(gamer_objlib PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
target_link_libraries(gamer_objlib PUBLIC casc tetstatic Eigen3::Eigen)

# SHARED LIBRARY
add_library(gamershared SHARED $<TARGET_OBJECTS:gamer_objlib>)
target_link_libraries(gamershared PUBLIC gamer_objlib)
target_compile_definitions(gamershared PRIVATE GAMER_DLL)
# set_target_properties(gamershared PROPERTIES OUTPUT_NAME gamer)

# STATIC LIBRARY
add_library(gamerstatic STATIC $<TARGET_OBJECTS:gamer_objlib>)
target_link_libraries(gamerstatic PUBLIC gamer_objlib)
# if(NOT WIN32)
#     # Shared and static libs will clobber each other on Windows
#     set_target_properties(gamerstatic PROPERTIES OUTPUT_NAME gamer)
# endif()

# Alias library names
set_target_properties(gamershared PROPERTIES OUTPUT_NAME gamer)
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

if(BUILD_BLENDGAMER)
    add_subdirectory(tools)
endif()

#####################################################################
# TESTING AND DOCUMENTATION
#####################################################################
if(GAMER_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Configure documentation builders
if(GAMER_DOCS)
    add_subdirectory(docs)
endif()
