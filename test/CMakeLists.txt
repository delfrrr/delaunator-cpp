add_executable(delaunator-tests
    src/delaunator-test.cpp
    src/vertex_points.cpp
    src/vertex_points.h
    src/triangle_indices.cpp
    src/triangle_indices.h    
)

target_link_libraries(delaunator-tests
    PRIVATE
    delaunator
)

target_include_directories(delaunator-tests
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${CMAKE_SOURCE_DIR}/includes/catch/single_include/catch2/
)
