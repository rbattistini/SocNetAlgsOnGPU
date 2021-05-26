add_executable(sna_bc sna_bc.cu
        common.cpp
        spmatops.cpp
        matio.cpp
        matstorage.cpp
        degree.cpp
        ecc.cpp
        device_props.cu
        cli.cu
        bc_par.cpp
        bc_statistics.cpp
        bc_we_kernels.cu
        bc_kernels_pitched.cu
        bc_kernels_compressed.cu
        bc.cpp
        cl_kernels.cu
        cl.cpp
        graphs.cpp)

set_target_properties(sna_bc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(sna_bc PRIVATE mmio)
target_link_libraries(sna_bc PRIVATE zf_log)
target_link_libraries(sna_bc PRIVATE libsnap.a)

if(OpenMP_CXX_FOUND)
    target_link_libraries(sna_bc PRIVATE OpenMP::OpenMP_CXX)
endif()

if(MPI_FOUND)
    target_link_libraries(sna_bc PRIVATE MPI::MPI_CXX)
endif()

target_link_libraries(sna_bc PRIVATE Boost::regex Boost::graph Boost::filesystem)
