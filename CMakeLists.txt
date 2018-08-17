cmake_minimum_required(VERSION 3.0.0)
project(delaunator VERSION 0.1.0)

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/includes")
    execute_process(COMMAND bash "-c" "(cd ${PROJECT_SOURCE_DIR} && ./fetch-includes.sh)")
endif()

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/includes/rapidjson/include")
# add_subdirectory("includes/rapidjson")

add_executable(delaunator src/main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
