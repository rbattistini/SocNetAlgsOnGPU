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

# from https://developer.nvidia.com/blog/building-cuda-applications-cmake/

add_library(mmio STATIC libs/libmmio/mmio.c libs/libmmio/mmio.h)
target_compile_features(mmio PUBLIC cxx_std_11)
set_target_properties(mmio PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

include_directories(libs/libmmio)

add_executable(betweenness src/betweenness.cu)
set_target_properties(
        betweenness PROPERTIES CUDA_SEPARABLE_COMPILATION ON
        COMPILE_FLAGS -DCUDA_DEBUG COMPILE_FLAGS -DDEBUG)
target_link_libraries(betweenness PRIVATE mmio)

file(GLOB samples samples/*.cu)

foreach(sample ${samples})
    GET_FILENAME_COMPONENT(samplename ${sample} NAME_WLE)
    add_executable(${samplename} ${sample})
endforeach()

add_executable(brandes_serial src/brandes_serial.cpp)
target_link_libraries(brandes_serial PRIVATE mmio)
