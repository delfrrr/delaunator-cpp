cmake_minimum_required(VERSION 3.0.0)
project(delaunator VERSION 0.1.0)
set (CMAKE_CXX_STANDARD 14)
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/includes")
    execute_process(COMMAND bash "-c" "(cd ${CMAKE_CURRENT_SOURCE_DIR} && ./fetch-includes.sh)")
endif()

#delaunator
add_library(delaunator src/delaunator.cpp)
target_include_directories (delaunator PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/includes/rapidjson/include")
target_include_directories (delaunator PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/includes/prettyprint")

#delaunator
add_library(json-helpers src/json-helpers.cpp)
target_include_directories (json-helpers PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/includes/rapidjson/include")

#delaunator-test
add_executable(delaunator-test src/delaunator-test.cpp)
target_link_libraries(delaunator-test delaunator)
target_include_directories (delaunator-test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/includes/catch/single_include/catch2")
target_link_libraries(delaunator-test json-helpers)
target_link_libraries(delaunator-test delaunator)


#triangulate
add_executable(triangulate src/triangulate.cpp)
target_include_directories (triangulate PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/includes/rapidjson/include")
target_include_directories (triangulate PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/includes/prettyprint")
target_link_libraries(triangulate delaunator)
target_link_libraries(triangulate json-helpers)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
