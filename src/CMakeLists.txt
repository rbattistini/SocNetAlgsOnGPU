add_executable(sna_bc sna_bc.cu
        common.cpp
        spmatops.cpp
        matio.cpp
        matstorage.cpp
        degree.cpp
        profiler.cpp
        device_props.cu
        cli.cpp
        bc.cpp
        bc_statistics.cpp
        bc_kernels.cu
        bc_kernels_pitched.cu
        spvb.cpp
        graphs.cpp)

target_link_libraries(sna_bc PRIVATE mmio)
target_link_libraries(sna_bc PRIVATE libsnap.a)
target_link_libraries(sna_bc PRIVATE zf_log)

if(OpenMP_CXX_FOUND)
    target_link_libraries(sna_bc PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(sna_bc PUBLIC Boost::regex Boost::graph Boost::filesystem)
    target_link_libraries(sna_bc PRIVATE MPI::MPI_CXX)
endif()
