cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU LANGUAGES C CXX CUDA
        DESCRIPTION "Social network Analysis algorithms for computing centrality metrics")

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif()

if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
    set(CMAKE_CUDA_ARCHITECTURES 61)
endif()
message(STATUS "CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}")

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

find_package(OpenMP REQUIRED)

add_subdirectory(lib/mmio)
add_subdirectory(lib/zf_log)
include_directories(lib/snap/lib)
include_directories(lib/snap/include)
include_directories(lib/doctest)
include_directories(include)

add_library(snap SHARED IMPORTED)
set_property(TARGET snap PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/lib/snap/lib/libsnap.so")

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(dataset)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CUDA_FLAGS_DEBUG "-g -lineinfo")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    message("-- CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
    message("-- CMAKE_CUDA_FLAGS_DEBUG is ${CMAKE_CUDA_FLAGS_DEBUG}")
endif()

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    message("-- CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
    message("-- CMAKE_CUDA_FLAGS_RELEASE is ${CMAKE_CUDA_FLAGS_RELEASE}")
endif ()
