cmake_minimum_required(VERSION 3.0)

project(Tetgen VERSION 1.5.0)

set(TETGEN_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(TETGEN_MASTER_PROJECT ON)
endif()

# Set some default options
option(BUILD_TETGEN_BIN "Build the tetgen binary?" ${TETGEN_MASTER_PROJECT})
option(TET_STATIC "Build the tetgen static library?" OFF)
option(TET_SHARED "Build the tetgen shared libary?" OFF)
option(SINGLE "Use single precision floating point numbers?" OFF)

# If added by add_subdirectory(... EXCLUDE_FROM_ALL) automatically expose shared
# and static library for building as dependencies.
get_directory_property(EXCLUDE_FROM_ALL DIRECTORY . EXCLUDE_FROM_ALL)
if (NOT TETGEN_MASTER_PROJECT AND EXCLUDE_FROM_ALL)
    set(TET_STATIC on)
    set(TET_SHARED on)
endif()

####################################################################
# Define various compile defs and flags for proper function
#####################################################################
set(COMPILE_DEFS "")
if(SINGLE)
    list(APPEND COMPILE_DEFS SINGLE)
endif(SINGLE)

# Set some definitions in accord with Shewchuk's recommendations
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    # Linux & gcc only... APPLE does not define fpu_control
    list(APPEND COMPILE_DEFS LINUX)
    set_source_files_properties(predicates.cxx PROPERTIES COMPILE_FLAGS -O0)
elseif(WIN32)
    list(APPEND COMPILE_DEFS CPU86)
    list(APPEND COMPILE_DEFS _MSC_VER)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        list(APPEND COMPILE_DEFS _WIN64)
    endif()
    set_source_files_properties(predicates.cxx PROPERTIES COMPILE_FLAGS /Od)
endif()

#####################################################################
# Define library targets
#####################################################################
if(TET_SHARED)
    add_library(tetshared SHARED tetgen.cxx predicates.cxx)
    target_compile_definitions(tetshared PRIVATE TETLIBRARY ${COMPILE_DEFS})
    target_include_directories(tetshared INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
    install(TARGETS tetshared DESTINATION lib)
endif()

if(TET_STATIC)
    add_library(tetstatic STATIC tetgen.cxx predicates.cxx)
    target_compile_definitions(tetstatic PRIVATE TETLIBRARY ${COMPILE_DEFS})
    target_include_directories(tetstatic INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
    install(TARGETS tetstatic DESTINATION lib)
endif()

#####################################################################
# Setup binary targets
#####################################################################
if(BUILD_TETGEN_BIN)
    add_executable(tetgen tetgen.cxx predicates.cxx)
    target_compile_definitions(tetgen PRIVATE ${COMPILE_DEFS})
    install(TARGETS tetgen DESTINATION bin)
endif()