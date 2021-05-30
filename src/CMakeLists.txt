add_executable(sna_bc sna_bc.cu
        common.cpp
        spmatops.cpp
        matio.cpp
        matds.cpp
        degree.cpp
        ecc.cpp
        device_props.cu
        cli.cu
        bc_statistics.cpp
        bc_we_kernels.cu
        bc_we_kernels_pitched.cu
        bc_ep_kernels.cu
        bc.cpp
        bc_vp_kernels.cu
        cl_kernels.cu
        cl.cpp
        graphs.cpp)

set_target_properties(sna_bc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(sna_bc PRIVATE mmio)
target_link_libraries(sna_bc PRIVATE zf_log)

#add_executable(par_bc main_pbgl_bc.cpp
#        bc_pbgl.cpp
#        matio.cpp
#        graphs.cpp
#        degree.cpp
#        spmatops.cpp
#        ecc.cpp
#        matds.cpp
#        common.cpp)
#
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(par_bc PRIVATE OpenMP::OpenMP_CXX)
#endif()
#
#if(MPI_FOUND)
#    target_link_libraries(par_bc PRIVATE MPI::MPI_CXX)
#endif()
#
#target_link_libraries(par_bc PRIVATE mmio)
#target_link_libraries(par_bc PRIVATE zf_log)
#target_link_libraries(par_bc PRIVATE Boost::regex Boost::graph Boost::filesystem)
