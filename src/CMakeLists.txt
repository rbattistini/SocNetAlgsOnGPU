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
        bc_we_kernel_nopitch.cu
        bc_we_kernel.cu
        bc_ep_kernel.cu
        bc.cpp
        bc_vp_kernel.cu
        cl_kernels.cu
        cl.cpp
        graphs.cpp)

set_target_properties(sna_bc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(sna_bc PRIVATE mmio)
target_link_libraries(sna_bc PRIVATE zf_log)
