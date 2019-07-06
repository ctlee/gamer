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

# - Find blender
# This module attempts to find Blender software
# FindBlender provides the following variables:
#
# BLENDER_FOUND                 - has Blender been found
# BLENDER_EXECUTABLE            - Path to Blender executable
# BLENDER_VERSION               - Version of Blender
# BLENDER_SCRIPT_PATH           - Location of user script installation path
# BLENDER_PYTHON_VERSION        - Version of Blender embedded python
# BLENDER_PYTHON_EXECUTABLE     - Path to Blender embedded python
# BLENDER_PYTHON_SIZEOF_VOID_P  - Number of bits for a void_p

find_program(BLENDER_EXECUTABLE blender)

if(BLENDER_EXECUTABLE)
    # Get the python path and version from blender
    execute_process(COMMAND ${BLENDER_EXECUTABLE} -b -noaudio --factory-startup --python-expr
        "import sys;import struct;import bpy;
print(str(bpy.app.version[0]) + '.' + str(bpy.app.version[1]) + '.' + str(bpy.app.version[2]))
print(bpy.utils.script_path_user())
print('.'.join(str(v) for v in sys.version_info[0:3]));
print(bpy.app.binary_path_python)
print(struct.calcsize('@P'))
"
                    RESULT_VARIABLE _BLENDER_SUCCESS
                    OUTPUT_VARIABLE _BLENDER_VALUES
                    ERROR_VARIABLE  _BLENDER_ERROR_VALUE)

    if(_BLENDER_SUCCESS MATCHES 0)
        string(REGEX REPLACE "\n" ";" _BLENDER_VALUES ${_BLENDER_VALUES})
        list(GET _BLENDER_VALUES 0 BLENDER_VERSION)
        list(GET _BLENDER_VALUES 1 BLENDER_SCRIPT_PATH)
        list(GET _BLENDER_VALUES 2 BLENDER_PYTHON_VERSION)
        list(GET _BLENDER_VALUES 3 BLENDER_PYTHON_EXECUTABLE)
        list(GET _BLENDER_VALUES 4 BLENDER_PYTHON_SIZEOF_VOID_P)

        string(REPLACE "." ";" _BLENDER_VERSION "${BLENDER_VERSION}")
        list(GET _BLENDER_VERSION 0 BLENDER_VERSION_MAJOR)
        list(GET _BLENDER_VERSION 1 BLENDER_VERSION_MINOR)
        list(GET _BLENDER_VERSION 2 BLENDER_VERSION_PATCH)

        string(REPLACE "." ";" _BPY_VERSION "${BLENDER_PYTHON_VERSION}")
        list(GET _BPY_VERSION 0 BLENDER_PYTHON_VERSION_MAJOR)
        list(GET _BPY_VERSION 1 BLENDER_PYTHON_VERSION_MINOR)
        list(GET _BPY_VERSION 2 BLENDER_PYTHON_VERSION_PATCH)
    elseif(NOT blender_FIND_QUIETLY)
            message(WARNING "Blender config failure:\n${_BLENDER_ERROR_VALUE}")
    endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(BLENDER
        REQUIRED_VARS
            BLENDER_EXECUTABLE
            BLENDER_VERSION
            BLENDER_SCRIPT_PATH
            BLENDER_PYTHON_VERSION
            BLENDER_PYTHON_EXECUTABLE
            BLENDER_PYTHON_SIZEOF_VOID_P
       VERSION_VAR BLENDER_VERSION
       HANDLE_COMPONENTS
    )