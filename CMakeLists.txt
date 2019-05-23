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

cmake_minimum_required(VERSION 3.11)
# Disable in source builds
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Add path to custom modules
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

# Override rules for MSVC static compiler flags
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_MODULE_PATH}/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_MODULE_PATH}/cxx_flag_overrides.cmake)

#####################################################################
# Project GAMer
#####################################################################
set(PRERELEASE ON)
set(PRE_VERSION beta)
project(GAMer VERSION 2.0.1)

# Require c++14 and standard libraries
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#####################################################################
# Options
#####################################################################
option(SINGLE "Use single precision floating point numbers?" OFF)
option(BUILD_STATIC_LIB "Build static libraries?" OFF)
option(BUILD_PYTHONEXT "Build GAMer python extension?" OFF)
option(BUILD_BLENDER "Build the GAMer addon for Blender?" OFF)
option(BLENDER_PLUGIN_INSTALL "Have CMake install the plugin?" OFF)
option(BUILD_TESTS "Build the GAMer tests?" OFF)
option(GETEIGEN "Download Eigen?" ON)

# Set conditional build
if(BUILD_BLENDER)
    set(BUILD_PYTHONEXT ON)
endif()

if(BUILD_PYTHONEXT)
    set(BUILD_STATIC_LIB ON)
endif()

# Add -fPIC to all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set build type to Release if none was defined
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

## Print debug flags
# message(STATUS "CMAKE_C_FLAGS is: ${CMAKE_C_FLAGS}")
# message(STATUS "CMAKE_C_FLAGS_DEBUG is: ${CMAKE_C_FLAGS_DEBUG}")
# message(STATUS "CMAKE_C_FLAGS_RELEASE is: ${CMAKE_C_FLAGS_RELEASE}")
# message(STATUS "CMAKE_C_FLAGS_RELWITHDEBINFO is: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
# message(STATUS "CMAKE_C_FLAGS_MINSIZEREL is: ${CMAKE_C_FLAGS_MINSIZEREL}")
# message(STATUS "CMAKE_CXX_FLAGS is: ${CMAKE_CXX_FLAGS}")
# message(STATUS "CMAKE_CXX_FLAGS_DEBUG is: ${CMAKE_CXX_FLAGS_DEBUG}")
# message(STATUS "CMAKE_CXX_FLAGS_RELEASE is: ${CMAKE_CXX_FLAGS_RELEASE}")
# message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO is: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
# message(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL is: ${CMAKE_CXX_FLAGS_MINSIZEREL}")
# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
# message(STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE})

# Define where to put the libraries and binaries
# set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# set(TEST_DATA_PATH ${PROJECT_SOURCE_DIR}/data)

add_subdirectory(libraries)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/libraries")

# List the include directories for debugging
# get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# foreach(dir ${dirs})
#     message(STATUS "GAMer include_dir='${dir}'")
# endforeach()

list(APPEND SOURCES
    "src/OFF_SurfaceMesh.cpp"
    "src/OBJ_SurfaceMesh.cpp"
    "src/SurfaceMesh.cpp"
    "src/SurfaceMeshDetail.cpp"
    "src/Vertex.cpp"
    "src/TetMesh.cpp"
    "src/PDBReader.cpp"
    "src/pdb2mesh.cpp"
)

add_library(gamershared SHARED ${SOURCES})
#message(STATUS "GAMer sources: ${SOURCES}")

target_link_libraries(gamershared tetshared casc eigen)
install(DIRECTORY include DESTINATION include FILES_MATCHING PATTERN "*.h")
install(TARGETS gamershared DESTINATION lib)

if(BUILD_STATIC_LIB)
    add_library(gamerstatic STATIC ${SOURCES})
    target_link_libraries(gamerstatic tetstatic casc eigen)
    install(TARGETS gamerstatic DESTINATION lib)
endif()

if(BUILD_PYTHONEXT)
    add_subdirectory(python)
    # add_subdirectory(swig)
endif()

if(BUILD_BLENDER)
    # Generate the name of the build
    set(ZIPNAME "gamer_addon-${GAMer_VERSION}-")
    if(PRERELEASE)
        string(CONCAT ZIPNAME ${ZIPNAME} "${PRE_VERSION}-")
    endif()
    string(CONCAT ZIPNAME ${ZIPNAME} "${CMAKE_HOST_SYSTEM_NAME}")

    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        if(WIN32)
            string(CONCAT ZIPNAME ${ZIPNAME} "64")
        elseif(CMAKE_SYSTEM_NAME STREQUAL Linux)
            string(CONCAT ZIPNAME ${ZIPNAME} "-x86_64")
        endif()
    else()
        if(WIN32)
            string(CONCAT ZIPNAME ${ZIPNAME} "32")
        elseif(CMAKE_SYSTEM_NAME STREQUAL Linux)
            string(CONCAT ZIPNAME ${ZIPNAME} "-i686")
        endif()
    endif()
    string(CONCAT ZIPNAME ${ZIPNAME} ".zip")

    configure_file(${CMAKE_SOURCE_DIR}/tools/blender_addon/__init__.in
                ${CMAKE_BINARY_DIR}/plugin/gamer_addon/__init__.py @ONLY)
    add_custom_target(zip_plugin
        ALL
        COMMAND ${CMAKE_COMMAND} -E tar "cfv"
            ../${ZIPNAME} --format=zip
            ./gamer_addon
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/plugin
        COMMENT "Zipping up addon components"
        VERBATIM
    )
    add_dependencies(zip_plugin stage_files)

    # Have CMake install Blender...
    if(BLENDER_PLUGIN_INSTALL)
        find_package(Blender REQUIRED)
        install(DIRECTORY ${CMAKE_BINARY_DIR}/plugin/ DESTINATION ${BLENDER_ADDON_PATH}/addons PATTERN "*")
    endif()
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Configure documentation building
add_subdirectory(docs)


