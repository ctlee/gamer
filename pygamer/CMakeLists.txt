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


execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "from distutils import sysconfig as s;import sys;import struct;
print(s.get_python_lib(plat_specific=True));
"
        RESULT_VARIABLE _PYTHON_SUCCESS
        OUTPUT_VARIABLE _PYTHON_VALUES
        ERROR_VARIABLE _PYTHON_ERROR_VALUE)

if(NOT _PYTHON_SUCCESS MATCHES 0)
    message(FATAL_ERROR
        "Python config failure:\n${_PYTHON_ERROR_VALUE}")
endif()


# Convert the process output into a list
string(REGEX REPLACE "\n" ";" _PYTHON_VALUES ${_PYTHON_VALUES})
list(GET _PYTHON_VALUES 0 PYTHON_SITE_PACKAGES)

file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
             INPUT ${CMAKE_SOURCE_DIR}/cmake-modules/__init__.py.in)

if(SKBUILD)
    install(TARGETS pygamer LIBRARY DESTINATION pygamer)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION pygamer)
else()
    message(STATUS "CMake install no longer installs PyGAMer to your Python site-packages. "
        "Please configure and install using 'python setup.py install' instead."
    )
    # Otherwise put in PREFIX/lib if conventionally installed
    # install(TARGETS pygamer LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pygamer")
    # install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION "${PYTHON_SITE_PACKAGES}/pygamer")
endif()