add_executable(parallel betweenness.cu
        gkernels.cu
        utils.cpp
        matio.cpp
        matstorage.cpp)

set_target_properties(
        parallel PROPERTIES CUDA_SEPARABLE_COMPILATION ON
        COMPILE_FLAGS -DCUDA_DEBUG COMPILE_FLAGS -DDEBUG)

target_link_libraries(parallel PRIVATE mmio)

add_executable(serial main.cpp
        utils.cpp
        matio.cpp
        matstorage.cpp
        spvb.cpp
        graphs.cpp)

target_link_libraries(serial PRIVATE mmio)

add_executable(graph_gen graph_gen.cpp
        utils.cpp
        matio.cpp
        matstorage.cpp)

target_link_libraries(graph_gen PRIVATE mmio)

if(OpenMP_CXX_FOUND)
    target_link_libraries(graph_gen PUBLIC OpenMP::OpenMP_CXX libsnap.a)
    target_link_libraries(serial PUBLIC OpenMP::OpenMP_CXX libsnap.a)
endif()
