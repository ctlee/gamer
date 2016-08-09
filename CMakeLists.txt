# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 3.6)

project(GAMer)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Werror")

message("CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")
message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
message("CMAKE_C_FLAGS_RELWITHDEBINFO is ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
message("CMAKE_C_FLAGS_MINSIZEREL is ${CMAKE_C_FLAGS_MINSIZEREL}")

message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message("CMAKE_CXX_FLAGS_MINSIZEREL is ${CMAKE_CXX_FLAGS_MINSIZEREL}")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(tetgen)
add_subdirectory(triangle)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(include)
# Hack to include maloc for now... TODO: fix
include_directories(/Users/ctlee/fetk/build/x86_64-apple-darwin15.6.0/include)

find_package(X11 QUIET)
if(X11_FOUND)
    include_directories(${X11_INCLUDE_DIR})
endif(X11_FOUND)

# List the include directories...
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
    message(STATUS "dir='${dir}'")
endforeach()

# Change this to not use GLOB
file(GLOB SOURCES "src/*.C")

add_library(gamer SHARED ${SOURCES})
set_target_properties(gamer PROPERTIES COMPILE_FLAGS "-std=c++14 -Wno-error")
target_link_libraries(gamer tet tri)

add_subdirectory(tools)
