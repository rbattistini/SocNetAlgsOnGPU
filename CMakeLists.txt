cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU CXX CUDA)

set(CMAKE_CUDA_FLAGS "-O3 -arch=sm_61")
# equal to -gencode arch=compute_60,code=sm_61

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

# from https://developer.nvidia.com/blog/building-cuda-applications-cmake/

add_library(libmmio STATIC mmio/mmio.cpp mmio/mmio.h)
target_compile_features(libmmio PUBLIC cxx_std_11)
set_target_properties(libmmio PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

include_directories(mmio)

add_executable(betweenness betweenness.cu)
set_property(TARGET betweenness PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(betweenness PRIVATE libmmio)

add_executable(bandwidthtest samples/bandwidthtest.cu)
add_executable(profile samples/profile.cu)
