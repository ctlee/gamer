add_executable(simplicial_complex_test simplicial_complex_test.cpp)
set_target_properties(simplicial_complex_test PROPERTIES COMPILE_FLAGS "-std=c++14")

# add_test(simplicial_complex_test ${EXECUTABLE_OUTPUT_PATH}/simplicial_complex_test ${TEST_DATA_PATH}/klein.off)


add_executable(orientable_test orientable_test.cpp)
set_target_properties(orientable_test PROPERTIES COMPILE_FLAGS "-std=c++14")

add_test(orientable_test_klein   ${EXECUTABLE_OUTPUT_PATH}/orientable_test ${TEST_DATA_PATH}/klein.off   false)
add_test(orientable_test_icosa   ${EXECUTABLE_OUTPUT_PATH}/orientable_test ${TEST_DATA_PATH}/icosa.off   true)
add_test(orientable_test_tetra   ${EXECUTABLE_OUTPUT_PATH}/orientable_test ${TEST_DATA_PATH}/tetra.off   true)
add_test(orientable_test_protein ${EXECUTABLE_OUTPUT_PATH}/orientable_test ${TEST_DATA_PATH}/protein.off true)

add_executable(readPDB_gauss readPDB_gauss.cpp)
set_target_properties(readPDB_gauss PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(readPDB_gauss gamershared)

#add_test(readPDB_gauss ${EXECUTABLE_OUTPUT_PATH}/readPDB_gauss ${TEST_DATA_PATH}/1CID.pdb)


add_executable(readPDB_molsurf readPDB_molsurf.cpp)
set_target_properties(readPDB_molsurf PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(readPDB_molsurf gamershared)

add_test(readPDB_molsurf ${EXECUTABLE_OUTPUT_PATH}/readPDB_molsurf ${TEST_DATA_PATH}/1CID.pdb)


add_executable(cmdtest cmdtest.cpp)
set_target_properties(cmdtest PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(cmdtest gamershared)

add_test(cmdtest ${EXECUTABLE_OUTPUT_PATH}/cmdtest)


add_executable(multivector_test multivector_test.cpp)
set_target_properties(multivector_test PROPERTIES COMPILE_FLAGS "-std=c++14")
target_link_libraries(multivector_test gamershared)

add_test(multivector_test ${EXECUTABLE_OUTPUT_PATH}/multivector_test)
