add_executable(simplicial_complex_test simplicial_complex_test.cpp)
set_target_properties(simplicial_complex_test PROPERTIES COMPILE_FLAGS "-std=c++14")

add_test(simplicial_complex_test ${EXECUTABLE_OUTPUT_PATH}/simplicial_complex_test ${TEST_DATA_PATH}/klein.off)

target_link_libraries(simplicial_complex_test)