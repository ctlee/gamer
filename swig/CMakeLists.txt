cmake_minimum_required(VERSION 3.6)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# First find the python interpreter
find_package(PythonInterp REQUIRED)
message(STATUS "Found Python binary: ${PYTHON_EXECUTABLE}")
message(STATUS "Python Version String: ${PYTHON_VERSION_STRING}")

# Then use it to inform pythonlibs
find_package(PythonLibs REQUIRED)
message(STATUS "PythonLibs Version String: ${PYTHONLIBS_VERSION_STRING}")
message(STATUS "Python Include Dir: ${PYTHON_INCLUDE_PATH}")
message(STATUS "Python Lib Dir: ${PYTHON_LIBRARIES}")

find_package(Numpy REQUIRED)

if(NOT ${PYTHON_VERSION_STRING} STREQUAL ${PYTHONLIBS_VERSION_STRING})
    message(FATAL_ERROR "Python interpreter and libs versions do not match")
endif()

find_package(Blender REQUIRED)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET(SWIG_MODULE_gamer_EXTRA_DEPS "../include/WrapSurfaceMesh.h" "../src/WrapSurfaceMesh.C")

SET_SOURCE_FILES_PROPERTIES(gamer.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(gamer.i PROPERTIES SWIG_FLAGS "-builtin")
SWIG_ADD_MODULE(gamer python gamer.i)
SWIG_LINK_LIBRARIES(gamer ${PYTHON_LIBRARIES} gamerstatic)

install(TARGETS _gamer DESTINATION lib/site-packages)
install(FILES ${CMAKE_BINARY_DIR}/swig/gamer.py DESTINATION lib/site-packages)
