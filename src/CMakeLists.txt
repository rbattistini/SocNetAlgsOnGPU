add_executable(parallel betweenness.cu
        gkernels.cu
        utils.cpp
        matio.cpp
        matstorage.cpp)

set_target_properties(parallel PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(parallel PRIVATE Boost::regex Boost::graph)
target_link_libraries(parallel PRIVATE mmio)

add_executable(serial main.cu
        utils.cpp
        spmatops.cpp
        matio.cpp
        matstorage.cpp
        spvb.cpp
        graphs.cpp)

target_link_libraries(serial PRIVATE mmio)
target_link_libraries(serial PRIVATE GSL)

add_executable(try try.cpp
        utils.cpp
        matio.cpp
        matstorage.cpp)

target_link_libraries(try PRIVATE mmio)

if(OpenMP_CXX_FOUND)
    target_link_libraries(serial PUBLIC OpenMP::OpenMP_CXX libsnap.a)
endif()
