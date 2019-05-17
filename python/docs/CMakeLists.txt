
add_custom_target(pygamer_docs
    COMMAND make html
    DEPENDS pygamer
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Sphinx"
    VERBATIM
    SOURCES ../src/SurfaceMesh.cpp
)