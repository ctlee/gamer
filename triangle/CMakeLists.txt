# CMake List to build Triangle
#
#

cmake_minimum_required(VERSION 3.6)

add_executable(triangle triangle.c)
add_library(tri SHARED triangle.c)
add_library(tristatic STATIC triangle.c)

# Tell the libraries they're libraries
set_target_properties(tri tristatic PROPERTIES "COMPILE_DEFINITIONS" TRILIBRARY)

# Set some definitions in accord with Shewchuk's recommendations
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    set_target_properties(triangle tri tristatic PROPERTIES 
            "COMPILE_DEFINITIONS" LINUX)
elseif(WIN32)
    set_target_properties(triangle tri tristatic PROPERTIES
            "COMPILE_DEFINITIONS" CPU86)
endif(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)

if(NOT UNIX)
set_target_properties(triangle tri tristatic PROPERTIES 
        "COMPILE_DEFINITIONS" NO_TIMER)
endif(NOT UNIX)
