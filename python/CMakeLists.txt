
PYBIND11_ADD_MODULE(pygamer MODULE
    src/SurfaceMesh.cpp)
target_link_libraries(pygamer PRIVATE gamerstatic)

add_subdirectory(docs)