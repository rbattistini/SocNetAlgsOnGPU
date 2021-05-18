add_executable(parallel betweenness.cu
        gkernels.cu
        common.cpp
        matio.cpp
        matstorage.cpp)

set_target_properties(parallel PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(parallel Boost::regex Boost::graph)
target_link_libraries(parallel mmio)
target_link_libraries(parallel zf_log)

add_executable(example_comp example_comp.cu
        common.cpp
        spmatops.cpp
        matio.cpp
        matstorage.cpp
        degree.cpp
        gkernels.cu
        spvb.cpp
        graphs.cpp)

target_link_libraries(example_comp PRIVATE mmio)
target_link_libraries(example_comp PRIVATE libsnap.a)
target_link_libraries(example_comp PRIVATE zf_log)

add_executable(example_bfs example_bfs.cu
        common.cpp
        spmatops.cpp
        matio.cpp
        matstorage.cpp
        degree.cpp
        gkernels.cu
        spvb.cpp
        graphs.cpp)

target_link_libraries(example_bfs PRIVATE mmio)
target_link_libraries(example_bfs PRIVATE libsnap.a)
target_link_libraries(example_bfs PRIVATE zf_log)

add_executable(example1 example1.cpp
        common.cpp
        matio.cpp
        degree.cpp
        graphs.cpp
        spmatops.cpp
        ecc.cpp
        matstorage.cpp)

target_link_libraries(example1 PRIVATE mmio)
target_link_libraries(example1 PRIVATE zf_log)

if(OpenMP_CXX_FOUND)
    target_link_libraries(example1 PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(example_comp PRIVATE OpenMP::OpenMP_CXX)
endif()

#define_file_basename_for_sources(try)
