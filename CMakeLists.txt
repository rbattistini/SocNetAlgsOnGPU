cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU C CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

find_package(OpenMP REQUIRED)

# from https://developer.nvidia.com/blog/building-cuda-applications-cmake/

add_library(mmio STATIC libs/mmio/mmio.c libs/mmio/mmio.h)
target_compile_features(mmio PUBLIC cxx_std_11)
set_target_properties(mmio PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
include_directories(libs/mmio)

include_directories(libs/snap-core)
include_directories(libs/glib-core)
include_directories(libs/doctest)
include_directories(include)

file(GLOB samples samples/*.cu)

foreach(sample ${samples})
    GET_FILENAME_COMPONENT(samplename ${sample} NAME_WLE)
    add_executable(${samplename} ${sample})
endforeach()

#file(GLOB serial_srcs src/*.cpp src/*.h)

add_executable(parallel src/betweenness.cu)
set_target_properties(
        parallel PROPERTIES CUDA_SEPARABLE_COMPILATION ON
        COMPILE_FLAGS -DCUDA_DEBUG COMPILE_FLAGS -DDEBUG)
target_link_libraries(parallel PRIVATE mmio)

add_executable(serial src/main.cpp
        include/utils.h src/utils.cpp
        include/matio.h src/matio.cpp
        include/matstorage.h src/matstorage.cpp
        include/spvb.h src/spvb.cpp
        src/graphs.cpp include/graphs.h)

target_link_libraries(serial PRIVATE mmio)

if(OpenMP_CXX_FOUND)
    target_link_libraries(serial PUBLIC OpenMP::OpenMP_CXX libsnap.a)
endif()

add_subdirectory(tests)
