
list(APPEND PYBIND11_SOURCES
    "src/pygamer.cpp"
    "src/SurfaceMesh.cpp"
    "src/Vertex.cpp"
    "src/SMFace.cpp"
    "src/SMSimplexID.cpp"
    "src/SMFunctions.cpp"
    )

PYBIND11_ADD_MODULE(pygamer MODULE ${PYBIND11_SOURCES})
# Link pygamer to gamer static
target_link_libraries(pygamer PRIVATE gamerstatic)