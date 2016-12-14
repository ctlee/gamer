cmake_minimum_required(VERSION 3.6)

project(GAMer VERSION 1.0.0)

option(BUILD_STATIC_LIB "Build the static libraries" OFF)
# BUILD_BLENDER will build static libraries as well
option(BUILD_BLENDER "Build the GAMer addon for Blender" OFF)

option(BUILD_TESTS "Build the test scripts" OFF)

if(BUILD_TESTS)
    enable_testing()
endif(BUILD_TESTS)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Print all build flags
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

#MESSAGE( STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )

# Define where to put the libraries and binaries
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(TEST_DATA_PATH ${PROJECT_SOURCE_DIR}/data)

add_subdirectory(libraries)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(include)

# List the include directories...
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
    message(STATUS "dir='${dir}'")
endforeach()

# Change this to not use GLOB
#file(GLOB SOURCES "src/*.C")
list(APPEND SOURCES 
    "src/CorrectNormals.C"
    "src/GemMesh.C"
    "src/GenerateSphereMesh.C"
    "src/MarchingCube.C"
    "src/OFF_SurfaceMesh.C"
    "src/OBJ_SurfaceMesh.C"
    "src/PDB2Mesh.C"
    "src/PDB2Volume.C"
    "src/ReadLattice.C"
    "src/ReadWriteMesh.C"
    "src/RegionGrowing.C"
    "src/SurfaceMesh.C"
    "src/SurfaceMeshOld.C"
    "src/SurfExtract.C"
    "src/SurfRefine.C"
    "src/SurfSmooth.C"
    "src/Triangulate.C"
    "src/Vertex.C"
    "src/WrapSimplicialComplex.C"
    "src/WrapSurfaceMesh.C"
)

add_library(gamershared SHARED ${SOURCES})

message(STATUS "SOURCES: ${SOURCES}")

set_target_properties(gamershared PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(gamershared tetshared trishared)
# TODO: populate the list of headers to install...
#install(FILES gamer.h DESTINATION include/gamer)
install(TARGETS gamershared DESTINATION lib/gamer)

if(BUILD_BLENDER OR BUILD_STATIC_LIB)
    add_library(gamerstatic STATIC ${SOURCES})
    
    set_target_properties(gamerstatic PROPERTIES COMPILE_FLAGS "-std=c++14")
    target_link_libraries(gamerstatic tetstatic tristatic)
    install(TARGETS gamerstatic DESTINATION lib/gamer)
endif(BUILD_BLENDER OR BUILD_STATIC_LIB)

if(BUILD_BLENDER)
    # build swig wrapper for python
    add_subdirectory(swig)
endif(BUILD_BLENDER)

add_subdirectory(tools)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)

message(STATUS "Doxygen ${DOXYGEN_EXECUTABLE}")
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
