
PYBIND11_ADD_MODULE(pygamer MODULE
    src/pygamer.cpp
    src/SurfaceMesh.cpp
    src/Vertex.cpp
    src/Face.cpp
    )
target_link_libraries(pygamer PRIVATE gamerstatic)