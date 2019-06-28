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


set(_BPY_STRICT FALSE)
if(BLENDER_PLUGIN_INSTALL OR BLENDER_VERSION_STRICT)
    set(_BPY_STRICT TRUE)

    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import sys;import struct;
print('.'.join(str(v) for v in sys.version_info[0:3]));
print(struct.calcsize('@P'));
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
    list(GET _PYTHON_VALUES 0 PYTHON_VERSION)
    list(GET _PYTHON_VALUES 1 PYTHON_SIZEOF_VOID_P)
endif()


# If strict then ensure system Python matches Blender's closely
if(_BPY_STRICT)
    find_package(Blender REQUIRED)
    set(_PYTHON_CMP_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
        )
    set(_BPY_CMP_VERSION "${BLENDER_PYTHON_VERSION_MAJOR}.${BLENDER_PYTHON_VERSION_MINOR}")
    if(NOT _PYTHON_CMP_VERSION VERSION_EQUAL _BPY_CMP_VERSION)
        message(FATAL_ERROR
            "Blender Addon Config Failure: "
            "Blender embedded Python version (${BLENDER_PYTHON_VERSION}) "
            "does not match the linked version of Python (${PYTHON_VERSION}).")
    elseif(NOT PYTHON_SIZEOF_VOID_P EQUAL BLENDER_PYTHON_SIZEOF_VOID_P)
        math(EXPR _PYTHON_BITS "${PYTHON_SIZEOF_VOID_P} * 8")
        math(EXPR _BLENDER_BITS "${BLENDER_PYTHON_SIZEOF_VOID_P} * 8")
        message(FATAL_ERROR
            "Blender Addon Config Failure: "
            "Linked Python is ${_PYTHON_BITS}-bit does not match "
            "${_BLENDER_BITS}-bit Blender).")
    endif()
else()
    find_package(Blender)
endif()

#####################################################################
# Stage files for Blender plugin packaging
#####################################################################
add_custom_target(stage_files
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pygamer> ${CMAKE_BINARY_DIR}/plugin/blendgamer/
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/tools/blender_addon/src ${CMAKE_BINARY_DIR}/plugin/blendgamer/
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/util.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/blendgamer.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/tetrahedralization.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/markers.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/util.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/versions.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/surfacemesh_ops.py
            ${CMAKE_SOURCE_DIR}/tools/blender_addon/src/ui.py
    COMMENT "Staging the GAMer Blender addon files"
    VERBATIM
)
add_dependencies(stage_files pygamer)

# Generate the name of the build
set(_ZIPNAME "blendgamer-${VERSION_SHORT}-")
if(VERSION_INFO)
    string(CONCAT _ZIPNAME ${_ZIPNAME} "${VERSION_INFO}-${VERSION_SHA1}-")
    if(VERSION_DIRTY)
        string(CONCAT _ZIPNAME ${_ZIPNAME} "${VERSION_DIRTY}-")
    endif()
endif()
string(CONCAT _ZIPNAME ${_ZIPNAME} "${CMAKE_HOST_SYSTEM_NAME}")

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    if(WIN32)
        string(CONCAT _ZIPNAME ${_ZIPNAME} "64")
    elseif(CMAKE_SYSTEM_NAME STREQUAL Linux)
        string(CONCAT _ZIPNAME ${_ZIPNAME} "-x86_64")
    endif()
else()
    if(WIN32)
        string(CONCAT _ZIPNAME ${_ZIPNAME} "32")
    elseif(CMAKE_SYSTEM_NAME STREQUAL Linux)
        string(CONCAT _ZIPNAME ${_ZIPNAME} "-i686")
    endif()
endif()
string(CONCAT _ZIPNAME ${_ZIPNAME} ".zip")

# Configure the Blender addon metadata
configure_file(${CMAKE_SOURCE_DIR}/tools/blender_addon/__init__.in
            ${CMAKE_BINARY_DIR}/plugin/blendgamer/__init__.py @ONLY)

# Target to zip up the plugin
add_custom_target(zip_plugin
    ALL
    COMMAND ${CMAKE_COMMAND} -E tar "cfv"
        ../${_ZIPNAME} --format=zip
        blendgamer
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/plugin
    COMMENT "Zipping up addon components"
    VERBATIM
)
add_dependencies(zip_plugin stage_files)

# Have CMake install Blender...
if(BLENDER_PLUGIN_INSTALL)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/plugin/ DESTINATION ${BLENDER_SCRIPT_PATH}/addons PATTERN "*")
endif()

### SAMPLE CODE TO VERIFY PLUGIN LOADS WITH BLENDER PYTHON
# execute_process(COMMAND "${BLENDER_PYTHON_EXECUTABLE}" "-c"
#     "import os; import sys; from importlib import util;
# sys.path.insert(0, os.path.abspath('${CMAKE_BINARY_DIR}/lib/'))
# import pygamer
# "
#         RESULT_VARIABLE _PYTHON_SUCCESS
#         OUTPUT_VARIABLE _PYTHON_VALUES
#         ERROR_VARIABLE _PYTHON_ERROR_VALUE)
