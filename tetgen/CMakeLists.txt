cmake_minimum_required(VERSION 3.6)

project(Tetgen VERSION 1.5.0)

# Set some default options
option(SINGLE "Use single precision" OFF)
option(BUILD_TETGEN_BIN "Build the tetgen binaries" OFF)
option(BUILD_STATIC_LIB "Build the static library" OFF)

set(TARGET_BINS "")

## Set the library targets
# The shared library is required for building standard GAMer.
add_library(tetshared SHARED tetgen.cxx predicates.cxx)
set(TARGET_LIBS tetshared)
install(FILES tetgen.h DESTINATION include/tetgen)

if(BUILD_BLENDER OR BUILD_STATIC_LIB)
    add_library(tetstatic STATIC tetgen.cxx predicates.cxx)
    list(APPEND TARGET_LIBS tetstatic)
endif(BUILD_BLENDER OR BUILD_STATIC_LIB)

install(TARGETS ${TARGET_LIBS} DESTINATION lib/tetgen)

## Set the binary targets
if(BUILD_TETGEN_BIN)
    add_executable(tetgen tetgen.cxx predicates.cxx)
    list(APPEND TARGET_BINS tetgen)
    install(TARGETS tetgen DESTINATION bin)
endif(BUILD_TETGEN_BIN)

list(APPEND ALL_TARGETS ${TARGET_LIBS})

# Tell the library targets they are libraries
set_target_properties(${TARGET_LIBS} PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)

if(SINGLE)
    set_target_properties(${TARGET_BINS} ${TARGET_LIBS} PROPERTIES 
            "COMPILE_DEFINITIONS" SINGLE)
endif(SINGLE)

# Set some definitions in accord with Shewchuk's recommendations
# TODO: Check if there are other cases... I suspect there are some.
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    set_target_properties(${TARGET_BINS} ${TARGET_LIBS} PROPERTIES 
            "COMPILE_DEFINITIONS" LINUX)
elseif(WIN32)
    set_target_properties(${TARGET_BINS} ${TARGET_LIBS} PROPERTIES
            "COMPILE_DEFINITIONS" CPU86)
    # TODO: define _MSC_VER  if microsoft visual C++ <2010 see tetgen.h
endif(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
