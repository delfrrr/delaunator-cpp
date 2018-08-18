cmake_minimum_required(VERSION 3.0.0)
project(delaunator VERSION 0.1.0)

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/includes")
    execute_process(COMMAND bash "-c" "(cd ${PROJECT_SOURCE_DIR} && ./fetch-includes.sh)")
endif()
# message("PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}")

add_executable(delaunator src/main.cpp)
target_include_directories (delaunator PRIVATE "${PROJECT_SOURCE_DIR}/includes/rapidjson/include")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
