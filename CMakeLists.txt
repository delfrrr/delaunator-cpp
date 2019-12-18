cmake_minimum_required(VERSION 3.8)
project(delaunator VERSION 0.4.0)
set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mason.cmake)

option(WERROR "Add -Werror flag to build (turns warnings into errors)" ON)
option(BENCHMARK_BIG_O "Calculate Big O in benchmark" OFF)
option(BENCHMARK_100M "Run against 100M points" OFF)
option(BENCHMARK_10M "Run against 100M points" OFF)
option(BUILD_BENCH_TEST "Use Mason" OFF)

# configure optimization
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(OPTIMIZATION_FLAGS "-O0 -DDEBUG")
    message("-- Configuring debug build")
else()
    set(OPTIMIZATION_FLAGS "-O3 -DNDEBUG")
    message("-- Configuring release build")
endif()

# Set relative install directories
set(INSTALL_INCLUDE_DIR "include")
set(INSTALL_CMAKE_DIR "lib${LIB_SUFFIX}/cmake/delaunator")

# Enable extra warnings to adhere to https://github.com/mapbox/cpp/issues/37
set(DESIRED_WARNINGS "-Wall -Wextra -Wconversion -Wunreachable-code -Wuninitialized -pedantic-errors -Wold-style-cast -Wno-error=unused-variable -Wshadow -Wfloat-equal -Weffc++")
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(DESIRED_WARNINGS "${DESIRED_WARNINGS} -Wmost")
endif()

# Note: -D_GLIBCXX_USE_CXX11_ABI=0 is needed to support mason packages that are precompiled libs
# Currently we only depend on a header only library, but this will help avoid issues when more libs are added via mason
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 ${DESIRED_WARNINGS}")

if (WERROR)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

if (BUILD_BENCH_TEST)
    # mason_use is a mason function within the mason.cmake file and provides ready-to-go vars, like "STATIC_LIBS" and "INCLUDE_DIRS"
    mason_use(catch VERSION 2.4.0 HEADER_ONLY)
    include_directories(SYSTEM ${MASON_PACKAGE_catch_INCLUDE_DIRS})

    mason_use(rapidjson VERSION 1.1.0 HEADER_ONLY)
    include_directories(SYSTEM ${MASON_PACKAGE_rapidjson_INCLUDE_DIRS})

    mason_use(benchmark VERSION 1.2.0)
    include_directories(SYSTEM ${MASON_PACKAGE_benchmark_INCLUDE_DIRS})
    include_directories("${PROJECT_SOURCE_DIR}/include")

    file(GLOB TEST_SOURCES test/*.cpp)
    add_executable(unit-tests ${TEST_SOURCES})

    # libbenchmark.a supports threads and therefore needs pthread support
    find_package(Threads REQUIRED)
    file(GLOB BENCH_SOURCES bench/*.cpp)
    add_executable(bench-tests ${BENCH_SOURCES})
    if(BENCHMARK_BIG_O)
        message("-- BENCHMARK_BIG_O=1")
        target_compile_definitions(bench-tests PUBLIC BENCHMARK_BIG_O=1)
    endif()
    if(BENCHMARK_100M)
        message("-- BENCHMARK_100M=1")
        target_compile_definitions(bench-tests PUBLIC BENCHMARK_100M=1)
    endif()
    if(BENCHMARK_10M)
        message("-- BENCHMARK_10M=1")
        target_compile_definitions(bench-tests PUBLIC BENCHMARK_10M=1)
    endif()
    
    #examples
    add_executable(triangulate-geojson examples/triangulate_geojson.cpp)
    add_executable(basic examples/basic.cpp)

    # link benchmark static library to the bench-tests binary so the bench tests know where to find the benchmark impl code
    target_link_libraries(bench-tests ${MASON_PACKAGE_benchmark_STATIC_LIBS} ${CMAKE_THREAD_LIBS_INIT})
endif()

add_library(delaunator INTERFACE)
target_include_directories(delaunator
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>)
    
install(TARGETS delaunator
        EXPORT delaunatorTargets)

add_library(delaunator::delaunator ALIAS delaunator)
 
 # Generate the cmake config and cmake config-version file:
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${delaunator_SOURCE_DIR}/delaunatorConfig.cmake.in"
    "${delaunator_BINARY_DIR}/delaunatorConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}
    PATH_VARS INSTALL_INCLUDE_DIR)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

install(EXPORT delaunatorTargets
        NAMESPACE delaunator::
        DESTINATION "${INSTALL_CMAKE_DIR}")

install(
    FILES "${delaunator_BINARY_DIR}/delaunatorConfig.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}" )
    
install(
	FILES "${delaunator_SOURCE_DIR}/include/delaunator.hpp"
	DESTINATION "${INSTALL_INCLUDE_DIR}" )
