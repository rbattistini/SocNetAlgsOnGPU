cmake_minimum_required(VERSION 3.9)
project(SocNetAlgsOnGPU  C CXX CUDA)

find_package(OpenMP REQUIRED)

set(CMAKE_CUDA_FLAGS "--gpu-architecture compute_60")
set(CMAKE_CUDA_STANDARD 11 REQUIRED)

file(GLOB cudafiles "*/cuda-*.cu")
file(GLOB cppfiles "*/*.cpp")
file(GLOB snap_examples "tutorials/*.cpp")
file(GLOB snap_core "snap-core/*.cpp" "snap-core/*.h")
file(GLOB glib_core "glib_core/*.cpp" "glib_core/*.h")

# Compile and link graph_gen
add_executable(graph_gen dataset/graph_gen.cpp)
target_include_directories(graph_gen PUBLIC snap-core)
target_include_directories(graph_gen PUBLIC glib-core)

if(OpenMP_CXX_FOUND)
    target_link_libraries(graph_gen PUBLIC OpenMP::OpenMP_CXX libsnap.a)
endif()

# Compile and link betweenness
add_executable(betweenness betweenness.cpp)
target_include_directories(betweenness PUBLIC snap-core)
target_include_directories(betweenness PUBLIC glib-core)

if(OpenMP_CXX_FOUND)
    target_link_libraries(betweenness PUBLIC OpenMP::OpenMP_CXX libsnap.a)
endif()

add_executable(graph_gen_boost dataset/graph_gen_boost.cpp)

foreach(cudafile ${cudafiles})
    GET_FILENAME_COMPONENT(cudafilename ${cudafile} NAME_WLE)
    add_executable(${cudafilename} ${cudafile})
endforeach()

#foreach(cppfile ${cppfiles})
#   GET_FILENAME_COMPONENT(cppfilename ${cppfile} NAME_WLE)
#    add_executable(${cppfilename} ${cppfile})
#endforeach()

foreach(snap_example ${snap_examples})
    GET_FILENAME_COMPONENT(snap_example_filename ${snap_example} NAME_WLE)
    add_executable(${snap_example_filename} ${snap_example})
    target_include_directories(${snap_example_filename} PRIVATE snap-core)
    target_include_directories(${snap_example_filename} PRIVATE glib-core)
endforeach()
