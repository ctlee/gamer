cmake_minimum_required(VERSION 3.6)

project(GAMer VERSION 1.0.0)

enable_testing()

option(BUILD_BLENDER "BUild the GAMer addon for Blender" OFF)
option(BUILD_STATIC_LIB "Build the static library" OFF)

# Set the flags for empty build and print all build flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Werror")
message(STATUS "CMAKE_C_FLAGS is: ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_DEBUG is: ${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE is: ${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CMAKE_C_FLAGS_RELWITHDEBINFO is: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_C_FLAGS_MINSIZEREL is: ${CMAKE_C_FLAGS_MINSIZEREL}")
message(STATUS "CMAKE_CXX_FLAGS is: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG is: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE is: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO is: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL is: ${CMAKE_CXX_FLAGS_MINSIZEREL}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Define where to put the libraries and binaries
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(tetgen)
add_subdirectory(triangle)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(include)
# Hack to include maloc for now... TODO: fix
include_directories(/Users/brogan/include)

# List the include directories...
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
    message(STATUS "dir='${dir}'")
endforeach()

# Change this to not use GLOB
file(GLOB SOURCES "src/*.C")

add_library(gamershared SHARED ${SOURCES})
set_target_properties(gamershared PROPERTIES COMPILE_FLAGS "-std=c++14 -Wno-error")
target_link_libraries(gamershared tetshared trishared)
# TODO: populate the list of headers to install...
#install(FILES gamer.h DESTINATION include/gamer)
install(TARGETS gamershared DESTINATION lib/gamer)

if(BUILD_BLENDER OR BUILD_STATIC_LIB)
    add_library(gamerstatic STATIC ${SOURCES})
    set_target_properties(gamerstatic PROPERTIES COMPILE_FLAGS "-std=c++14 -Wno-error")
    target_link_libraries(gamerstatic tetstatic tristatic)
    install(TARGETS gamerstatic DESTINATION lib/gamer)
endif(BUILD_BLENDER OR BUILD_STATIC_LIB)

if(BUILD_BLENDER)
    # TODO: actually build the swigs and blender stuff
endif(BUILD_BLENDER)

add_subdirectory(tools)

add_subdirectory(tests)
