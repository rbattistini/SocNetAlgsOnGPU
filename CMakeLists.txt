cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU C CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS regex graph)

add_library(mmio STATIC lib/mmio/mmio.c lib/mmio/mmio.h)

# from https://developer.nvidia.com/blog/building-cuda-applications-cmake/
set_target_properties(mmio PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

include_directories(lib/mmio)
include_directories(lib/snap-core)
include_directories(lib/glib-core)
include_directories(lib/doctest)
include_directories(include)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#include(cmake/code-coverage.cmake)
#include(cmake/sanitizers.cmake)
#include(cmake/tools.cmake)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(dataset)

#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")
#set(CMAKE_C_FLAGS_RELEASE "-O2 -DBENCHMARK")

# overwrote on purpose
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

#message("-- CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
#message("-- CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")

message("-- CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("-- CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
