# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 3.6)

# Include gamer headers
include_directories(${CMAKE_SOURCE_DIR}/include)
add_executable(generatemesh GenerateMesh.C)

set_target_properties(generatemesh PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(generatemesh gamershared)

