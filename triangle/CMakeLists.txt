# CMake List to build trinagle
cmake_minimum_required(VERSION 3.6)

add_executable(triangle triangle.c)

find_package(X11 QUIET)
if(X11_FOUND)
    include_directories(${X11_INCLUDE_DIR})
    add_executable(showme showme.c)
    target_link_libraries(showme ${X11_LIBRARIES})
    set_target_properties(showme PROPERTIES COMPILE_FLAGS "-Wno-error")
    #MESSAGE ( STATUS " X11_INCLUDE_DIR: " ${X11_INCLUDE_DIR} )
    #MESSAGE ( STATUS " X11_LIBRARIES: " ${X11_LIBRARIES} )
endif(X11_FOUND)

add_library(tri SHARED triangle.c)
add_library(tristatic STATIC triangle.c)


set_target_properties(tri tristatic PROPERTIES 
        "COMPILE_DEFINITIONS" ANSI_DECLARATORS)

target_compile_definitions(tri PRIVATE REAL=double VOID=void)
target_compile_definitions(tristatic PRIVATE REAL=double VOID=void)
# Tell the libraries they're libraries
set_target_properties(tri tristatic PROPERTIES "COMPILE_DEFINITIONS" TRILIBRARY)

# Set some definitions in accord with Shewchuk's recommendations
# TODO: Check if there are other cases... I suspect there are some.
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    set_target_properties(triangle showme tri tristatic PROPERTIES 
            "COMPILE_DEFINITIONS" LINUX)
elseif(WIN32)
    set_target_properties(triangle showme tri tristatic PROPERTIES
            "COMPILE_DEFINITIONS" CPU86)
endif(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)

if(NOT UNIX)
set_target_properties(triangle showme tri tristatic PROPERTIES 
        "COMPILE_DEFINITIONS" NO_TIMER)
endif(NOT UNIX)

target_include_directories(tri PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(FILES triangle.h DESTINATION include/triangle)
install(TARGETS tri tristatic DESTINATION lib/triangle)
install(TARGETS triangle showme DESTINATION bin)
