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

project(Triangle VERSION 1.6)

# Set some default options
option(SINGLE "Use single precision" OFF)
option(BUILD_TRIANGLE_BIN "Build the triangle binaries" OFF)
option(BUILD_STATIC_LIB "Build the static library" OFF)
option(SELF_CHECK "Insert lots of self checks" OFF)

set(TARGET_BINS "")
set(TARGET_LIBS "")
set(LIB_COMPILE_DEFS "")
set(COMPILE_DEFS "")

#####################################################################
# Define library targets
#####################################################################
list(APPEND LIB_COMPILE_DEFS ANSI_DECLARATORS)  # ANSI Declarations
list(APPEND LIB_COMPILE_DEFS TRILIBRARY)        # Build library version
list(APPEND LIB_COMPILE_DEFS REDUCED)           # See triangle.c
list(APPEND LIB_COMPILE_DEFS CDT_ONLY)          # See triangle.c


add_library(trishared SHARED triangle.c)
list(APPEND ${TARGET_LIBS} trishared)
# install(FILES triangle.h DESTINATION include/triangle)

if(BUILD_BLENDER OR BUILD_STATIC_LIB)
    add_library(tristatic STATIC triangle.c)
    list(APPEND TARGET_LIBS tristatic)
endif(BUILD_BLENDER OR BUILD_STATIC_LIB)
install(TARGETS ${TARGET_LIBS} DESTINATION lib)

#####################################################################
# Setup binary targets
#####################################################################
if(BUILD_TRIANGLE_BIN)
    add_executable(triangle triangle.c)
    list(APPEND TARGET_BINS triangle)

    find_package(X11 QUIET)
    if(X11_FOUND)
        include_directories(${X11_INCLUDE_DIR})
        add_executable(showme showme.c)
        target_link_libraries(showme ${X11_LIBRARIES})
        set_target_properties(showme PROPERTIES COMPILE_FLAGS "-Wno-error")
        list(APPEND TARGET_BINS showme)
    else()
        message(WARNING "Could not find libX11 skipping building Triangle's \
                visualization package: showme")
    endif(X11_FOUND)
    install(TARGETS ${TARGET_BINS} DESTINATION bin)
endif(BUILD_TRIANGLE_BIN)


#####################################################################
# Define various compile defs and flags for proper function
#####################################################################
list(APPEND COMPILE_DEFS NO_TIMER) # Skip timing code
if(SINGLE)
    list(APPEND COMPILE_DEFS SINGLE)
endif(SINGLE)

if(SELF_CHECK)
    list(APPEND COMPILE_DEFS SELF_CHECK)
endif(SELF_CHECK)

# Set some definitions in accord with Shewchuk's recommendations
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    # Linux && gcc only... APPLE does not define fpu_control
    list(APPEND COMPILE_DEFS LINUX)
elseif(WIN32)
    list(APPEND COMPILE_DEFS CPU86)
endif()


foreach(LIB IN LISTS TARGET_LIBS)
    target_compile_definitions(${LIB} PRIVATE ${LIB_COMPILE_DEFS} ${COMPILE_DEFS})
endforeach(LIB IN LISTS TARGET_LIBS)
