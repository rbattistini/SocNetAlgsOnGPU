cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU LANGUAGES C CXX CUDA
        DESCRIPTION "Social network Analysis algorithms for computing centrality metrics")

include(FetchContent)

#FetchContent_Declare(Snap
#    GIT_REPOSITORY "https://github.com/snap-stanford/snap")
#FetchContent_MakeAvailable(Snap)
#
#add_custom_target(Snap COMMAND make)

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
find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem regex graph)

add_subdirectory(lib/mmio)
add_subdirectory(lib/zf_log)
include_directories(lib/snap/snap-core)
include_directories(lib/snap/glib-core)
include_directories(lib/doctest)
include_directories(include)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(dataset)

# overwrote on purpose
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    message("-- CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
    message("-- CMAKE_CUDA_FLAGS_DEBUG is ${CMAKE_CUDA_FLAGS_DEBUG}")
endif()

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    message("-- CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
    message("-- CMAKE_CUDA_FLAGS_RELEASE is ${CMAKE_CUDA_FLAGS_RELEASE}")
endif ()
