# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 3.6)

# Add an executable to the project using the specified source files.
add_executable(tetgen tetgen.cxx predicates.cxx)

#Add a library to the project using the specified source files. 
add_library(tet SHARED tetgen.cxx predicates.cxx)
add_library(tetstatic STATIC tetgen.cxx predicates.cxx)

#Set properties on a target. 
#We use this here to set -DTETLIBRARY for when compiling the
#library
set_target_properties(tet PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)
set_target_properties(tetstatic PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)
