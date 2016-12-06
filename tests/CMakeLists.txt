find_package(Threads REQUIRED)
include(ExternalProject)

ExternalProject_Add(
    gtest
    # URL https://github.com/google/googletest/archive/master.zip
    URL https://github.com/google/googletest/archive/release-1.8.0.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

message(STATUS "source_dir ${source_dir}")
message(STATUS "binary_dir ${binary_dir}")
# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(libgtest PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# Create a libgmock target to be used as a dependency by test programs
add_library(libgmock IMPORTED STATIC GLOBAL)
add_dependencies(libgmock gtest)

# Set libgmock properties
set_target_properties(libgmock PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
include_directories("${source_dir}/googletest/include"
                    "${source_dir}/googlemock/include")




add_executable(objecttests main.cpp vertex_test.cpp tensor_test.cpp SurfaceMesh_test.cpp)
set_target_properties(objecttests PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(objecttests gamershared libgtest)

add_test(objecttests ${EXECUTABLE_OUTPUT_PATH}/objecttests)


add_executable(scratchtest scratchtest.cpp)
set_target_properties(scratchtest PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(scratchtest gamershared)

add_executable(scratchtest_jbm scratchtest_jbm.cpp)
set_target_properties(scratchtest_jbm PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(scratchtest_jbm gamershared)

add_executable(cmdtest cmdtest.cpp)
set_target_properties(cmdtest PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(cmdtest gamershared)

add_executable(cmdtest2 cmdtest2.cpp)
set_target_properties(cmdtest2 PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(cmdtest2 gamershared)
