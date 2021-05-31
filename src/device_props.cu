/****************************************************************************
 * @file device_props.cu
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Utility functions for NVIDIA GPUs device properties querying
 *
 * Copyright 2021 (c) 2021 by Riccardo Battistini
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ****************************************************************************/

#include "device_props.cuh"

#ifndef DEVICE_ID
#define DEVICE_ID 0
#endif

void print_separator() {
    for (int i = 0; i < LINE_LENGTH; i++)
        printf("-");

    printf("\n");
}

int get_runtime_version() {
    int driverVersion;
    cudaSafeCall(cudaRuntimeGetVersion(&driverVersion));
    return driverVersion;
}

size_t get_global_mem_size() {
    size_t totalMem;
    cudaSafeCall(cudaMemGetInfo(nullptr, &totalMem));
    return totalMem;
}

int get_compute_capability_major() {
    int devId, computeCap;
    cudaSafeCall(cudaGetDevice(&devId));
    cudaSafeCall(cudaDeviceGetAttribute(&computeCap,
                                        cudaDevAttrComputeCapabilityMajor, devId));
    return computeCap;
}

int get_compute_capability_minor() {
    int devId, computeCap;
    cudaSafeCall(cudaGetDevice(&devId));
    cudaSafeCall(cudaDeviceGetAttribute(&computeCap,
                                        cudaDevAttrComputeCapabilityMinor, devId));
    return computeCap;
}

unsigned int get_max_threads_per_block() {
    int devId, threadsPerBlock;
    cudaSafeCall(cudaGetDevice(&devId));
    cudaSafeCall(cudaDeviceGetAttribute(&threadsPerBlock,
                                        cudaDevAttrMaxThreadsPerBlock, devId));
    return (unsigned int) threadsPerBlock;
}

unsigned int get_sm_count() {
    int devId, numProcs;
    cudaSafeCall(cudaGetDevice(&devId));
    cudaSafeCall(cudaDeviceGetAttribute(&numProcs,
                                        cudaDevAttrMultiProcessorCount, devId));
    return (unsigned int) numProcs;
}

void set_device(int device_id) {
    cudaSafeCall(cudaSetDevice(device_id));
}

int get_device_count() {
    int count = 0;
    cudaSafeCall(cudaGetDeviceCount(&count));
    return count;
}

void print_gpu_overview(int device_id) {

    cudaDeviceProp prop{};
    cudaSafeCall(cudaGetDeviceProperties(&prop, device_id));

    printf("GPU overview:\n\n");
    printf("\tDevice: \t\t%s\n", prop.name);
    printf("\tCUDA Runtime version: \t%d\n", get_runtime_version());
    printf("\tCompute Capability: \t%d.%d\n",
           get_compute_capability_major(),
           get_compute_capability_minor());
    printf("\tGlobal Memory size: \t%.2f GB\n",
           get_global_mem_size() / (double) (1 << 30));
    printf("\t# of SMs: \t\t%d\n", get_sm_count());

    print_separator();
}
