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

find_package(Threads REQUIRED)

include(ExternalProject)
ExternalProject_Add(
    gtest
    # URL https://github.com/google/googletest/archive/master.zip
    URL https://github.com/google/googletest/archive/release-1.8.0.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)
# message(STATUS "source_dir ${source_dir}")
# message(STATUS "binary_dir ${binary_dir}")

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(libgtest PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# Create a libgmock target to be used as a dependency by test programs
# add_library(libgmock IMPORTED STATIC GLOBAL)
# add_dependencies(libgmock gtest)

# Set libgmock properties
# set_target_properties(libgmock PROPERTIES
#     "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
#     "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
# )

# I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
include_directories("${source_dir}/googletest/include")
                    # "${source_dir}/googlemock/include")

add_executable(objecttests main.cpp VertexTest.cpp tensorTest.cpp SurfaceMeshTest.cpp)
target_link_libraries(objecttests gamershared libgtest)
add_test(objecttests ${EXECUTABLE_OUTPUT_PATH}/objecttests)


add_executable(scratchtest scratchtest.cpp)
target_link_libraries(scratchtest gamershared)
# Profiling with gperftools
# target_link_libraries(scratchtest gamershared profiler)

add_executable(masonscratch masonscratch.cpp)
set_target_properties(masonscratch PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(masonscratch gamershared)