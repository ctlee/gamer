cmake_minimum_required(VERSION 3.0)
project(Tetgen VERSION 1.5.0)

# Set some default options
option(SINGLE "Use single precision" OFF)
option(BUILD_TETGEN_BIN "Build the tetgen binaries" OFF)
option(BUILD_STATIC_LIB "Build the static library" OFF)

set(TARGET_BINS "")
set(TARGET_LIBS "")
set(COMPILE_DEFS "")

#####################################################################
# Define library targets
#####################################################################
add_library(tetshared SHARED tetgen.cxx predicates.cxx)
list(APPEND TARGET_LIBS tetshared)
# install(FILES tetgen.h DESTINATION include/tetgen)

if(BUILD_BLENDER OR BUILD_STATIC_LIB)
    add_library(tetstatic STATIC tetgen.cxx predicates.cxx)
    list(APPEND TARGET_LIBS tetstatic)
endif()
install(TARGETS ${TARGET_LIBS} DESTINATION lib)

#####################################################################
# Setup binary targets
#####################################################################
if(BUILD_TETGEN_BIN)
    add_executable(tetgen tetgen.cxx predicates.cxx)

    list(APPEND TARGET_BINS tetgen)
    install(TARGETS tetgen DESTINATION bin)
endif()

#####################################################################
# Define various compile defs and flags for proper function
#####################################################################
list(APPEND COMPILE_DEFS NDEBUG) # No internal debug for all targets
if(SINGLE)
    list(APPEND COMPILE_DEFS SINGLE)
endif(SINGLE)

# Set some definitions in accord with Shewchuk's recommendations
if(CMAKE_SYSTEM_NAME STREQUAL Linux AND CMAKE_COMPILER_IS_GNUCC)
    # Linux && gcc only... APPLE does not define fpu_control
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

# Set compile defs for binary targets
foreach(BIN IN LISTS TARGET_BINS)
    target_compile_definitions(${BIN} PRIVATE ${COMPILE_DEFS})
endforeach(BIN IN LISTS TARGET_BINS)

# Set compile defs for library targets
foreach(LIB IN LISTS TARGET_LIBS)
    target_compile_definitions(${LIB} PRIVATE TETLIBRARY ${COMPILE_DEFS})
endforeach(LIB IN LISTS TARGET_LIBS)