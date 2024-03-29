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

if(${CMAKE_VERSION} VERSION_LESS 3.11)
    include(FetchContentLocal)
else()
    include(FetchContent)
endif()

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG origin/main
  GIT_SHALLOW TRUE
  SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/googletest-src" 
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
)

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(GAMER_TEST_SOURCES
    "main.cpp" 
    "VertexTest.cpp" 
    "tensorTest.cpp" 
    "SurfaceMeshTest.cpp"
    "tetrahedralizationTest.cpp"
)

add_executable(objecttests ${GAMER_TEST_SOURCES})
target_link_libraries(objecttests gamerstatic gtest_main)
# target_compile_options(objecttests PRIVATE -Werror -Wall -Weverything
#           -Wextra -pedantic-errors -Wconversion -Wsign-conversion
#           -ferror-limit=3 -Wno-c++17-extensions)
# -ferror-limit=3 -Wno-unused-variable -Wno-unused-parameter)


add_test(NAME ObjectTests COMMAND objecttests)

# Configure testing of Python module
find_package(pytest)
if(NOT PYTEST_FOUND AND BUILD_PYGAMER)
    message(WARNING "pytest could not be found! Please install pytest to run python tests.")
elseif(BUILD_PYGAMER)
    # Add a test for pytest
    add_test(NAME python-tests
      COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}
      WORKING_DIRECTORY $<TARGET_FILE_DIR:pygamer>
    )
endif()


add_executable(scratchtest scratchtest.cpp)
target_link_libraries(scratchtest gamerstatic)
# Profiling with gperftools
# target_link_libraries(scratchtest gamershared profiler)
