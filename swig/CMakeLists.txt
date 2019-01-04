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

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# First find the python interpreter
find_package(PythonInterp REQUIRED)
# message(STATUS "Found Python binary: ${PYTHON_EXECUTABLE}")
# message(STATUS "Python Version String: ${PYTHON_VERSION_STRING}")

# Then use it to inform pythonlibs
find_package(PythonLibs REQUIRED)
message(STATUS "PythonLibs Version String: ${PYTHONLIBS_VERSION_STRING}")
message(STATUS "Python Include Dir: ${PYTHON_INCLUDE_PATH}")
message(STATUS "Python Lib Dir: ${PYTHON_LIBRARIES}")

if(NOT ${PYTHON_VERSION_STRING} VERSION_EQUAL ${PYTHONLIBS_VERSION_STRING})
    message(FATAL_ERROR "Python interpreter and libs versions do not match")
endif()

execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.maxsize > 2**32)"
    OUTPUT_VARIABLE is_interp_x64
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT ${is_interp_x64})
    message(FATAL_ERROR "Python interpreter is not 64bits.")
endif()

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
# List the include directories...
# get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# foreach(dir ${dirs})
#     message(STATUS "SWIG include_dir='${dir}'")
# endforeach()

SET_SOURCE_FILES_PROPERTIES(pygamer.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(pygamer.i PROPERTIES SWIG_FLAGS "-O")

if(CMAKE_VERSION VERSION_LESS 3.8.0)
    SWIG_ADD_MODULE(pygamer python pygamer.i)
else()
    SWIG_ADD_LIBRARY(pygamer LANGUAGE python SOURCES pygamer.i)
endif()

# Handle some problems with Anaconda dylib linking in OSX
if(APPLE)
    set_target_properties(${SWIG_MODULE_pygamer_REAL_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
elseif(WIN32)
    swig_link_libraries(pygamer ${PYTHON_LIBRARIES})
else()
    # Nothing yet
endif()

SWIG_LINK_LIBRARIES(pygamer gamerstatic)

install(TARGETS _pygamer DESTINATION lib/site-packages)
install(FILES ${CMAKE_BINARY_DIR}/swig/pygamer.py DESTINATION lib/site-packages)


#####################################################################
# Have CMake install Blender...
#####################################################################
add_custom_target(stage_files
    COMMAND ${CMAKE_COMMAND} -E make_directory plugin/gamer
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/swig/__init__.py plugin/gamer/
    COMMAND ${CMAKE_COMMAND} -E copy swig/pygamer.py plugin/gamer/
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_pygamer> plugin/gamer/
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/tools/blender_addon/src ./plugin/gamer_addon
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/util.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/gameraddon_main.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/tetrahedralization.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/markers.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/util.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/versions.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/surfacemesh_ops.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/ui.py
    COMMENT "Staging the GAMer Blender addon files"
    VERBATIM
)
add_dependencies(stage_files _pygamer)