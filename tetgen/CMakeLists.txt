# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 3.6)

# Add an executable to the project using the specified source files.
add_executable(tetgen tetgen.cxx predicates.cxx)

#Add a library to the project using the specified source files. 
add_library(tet SHARED tetgen.cxx predicates.cxx)
add_library(tetstatic STATIC tetgen.cxx predicates.cxx)

# These targets are libraries
set_target_properties(tet tetstatic PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)
target_include_directories(tet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(FILES tetgen.h DESTINATION include/gamer)
install(TARGETS tet tetstatic DESTINATION lib/gamer)
install(TARGETS tetgen DESTINATION bin)
