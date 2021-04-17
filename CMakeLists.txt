cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(SocNetAlgsOnGPU C CXX CUDA)

# optimizations are disabled by default
set(CMAKE_CUDA_FLAGS "-arch=sm_61")
#set(CMAKE_CUDA_FLAGS "-O3 -arch=sm_61")
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

add_library(libmmio STATIC libs/mmio/mmio.c libs/mmio/mmio.h)
target_compile_features(libmmio PUBLIC cxx_std_11)
set_target_properties(libmmio PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

include_directories(libs/mmio)

# debugging is enabled by default
add_executable(betweenness src/betweenness.cu)
set_target_properties(
        betweenness PROPERTIES CUDA_SEPARABLE_COMPILATION ON
        COMPILE_FLAGS -DCUDA_DEBUG COMPILE_FLAGS -DDEBUG)
target_link_libraries(betweenness PRIVATE libmmio)

file(GLOB samples samples/*.cu)

foreach(sample ${samples})
    GET_FILENAME_COMPONENT(samplename ${sample} NAME_WLE)
    add_executable(${samplename} ${sample})
endforeach()

add_executable(brandes_serial src/brandes_serial.cpp)
target_link_libraries(brandes_serial PRIVATE libmmio)

add_executable(bsf_test src/bfs.cpp)
target_link_libraries(bsf_test PRIVATE libmmio)

add_executable(bsf src/bfs.c)
