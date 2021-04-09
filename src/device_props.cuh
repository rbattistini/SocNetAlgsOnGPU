/****************************************************************************
 *
 * device_props.h - Utility functions for NVIDIA GPUs device properties querying
 *
 * Copyright 2021 (c) 2021 by Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
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
 *
 ****************************************************************************/

#ifndef DEVICE_PROPERTIES_CUH
#define DEVICE_PROPERTIES_CUH

#ifdef __CUDACC__

int get_driver_version() {
    int driverVersion;
    cudaSafeCall( cudaDriverGetVersion(&driverVersion) );
    return driverVersion;
}

/*
 * Fast property querying.
 *
 * from https://gist.github.com/teju85/9521e2224f0c31f71a93b593ff64e8da
 */
int get_compute_capability() {
    int devId, computeCap;
    cudaSafeCall( cudaGetDevice(&devId) );
    cudaSafeCall( cudaDeviceGetAttribute(&computeCap,
                                         cudaDevAttrComputeCapabilityMajor, devId) );
    return computeCap;
}

int get_max_shared_mem() {
    int devId, smemSize;
    cudaSafeCall( cudaGetDevice(&devId) );
    cudaSafeCall( cudaDeviceGetAttribute(&smemSize,
                                         cudaDevAttrMaxSharedMemoryPerBlock, devId) );
    return smemSize;
}

int get_sm_count() {
    int devId, numProcs;
    cudaSafeCall( cudaGetDevice(&devId) );
    cudaSafeCall( cudaDeviceGetAttribute(&numProcs,
                                         cudaDevAttrMultiProcessorCount, devId) );
    return numProcs;
}

#endif

#endif //DEVICE_PROPERTIES_CUH
