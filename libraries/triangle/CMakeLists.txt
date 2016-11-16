# CMake List to build trinagle
cmake_minimum_required(VERSION 3.6)

project(Triangle VERSION 1.6)

# Set some default options
option(SINGLE "Use single precision" OFF)
option(BUILD_TRIANGLE_BIN "Build the triangle binaries" OFF)
option(BUILD_STATIC_LIB "Build the static library" OFF)

set(TARGET_BINS "")
set(BIN_COMPILE_DEFS "")
set(LIB_COMPILE_DEFS "")
set(ALL_COMPILE_DEFS "")

## Set the library targets
# The shared library is required for building standard GAMer.
add_library(trishared SHARED triangle.c)
set(TARGET_LIBS trishared)
install(FILES triangle.h DESTINATION include/triangle)

if(BUILD_BLENDER OR BUILD_STATIC_LIB)
    add_library(tristatic STATIC triangle.c)
    list(APPEND TARGET_LIBS tristatic)
endif(BUILD_BLENDER OR BUILD_STATIC_LIB)

list(APPEND LIB_COMPILE_DEFS ANSI_DECLARATORS)
install(TARGETS ${TARGET_LIBS} DESTINATION lib/triangle)

## Set the binary targets
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


# Tell the libraries they're libraries
list(APPEND LIB_COMPILE_DEFS TRILIBRARY)

if(SINGLE)
    list(APPEND LIB_COMPILE_DEFS SINGLE)
endif(SINGLE)

# Set some definitions in accord with Shewchuk's recommendations
# TODO: Check if there are other cases... I suspect there are some.
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    list(APPEND LIB_COMPILE_DEFS LINUX)
elseif(WIN32)
    list(APPEND LIB_COMPILE_DEFS CPU86)
endif(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)

if(NOT UNIX)
    list(APPEND LIB_COMPILE_DEFS NO_TIMER)
endif(NOT UNIX)

foreach(LIB IN LISTS TARGET_LIBS)
    target_compile_definitions(${LIB} PRIVATE ${LIB_COMPILE_DEFS})
endforeach(LIB IN LISTS TARGET_LIBS)
