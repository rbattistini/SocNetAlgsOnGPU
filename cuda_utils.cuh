/****************************************************************************
 *
 * cuda_utils.h - Utility functions for NVIDIA GPUs error checking and device
 * management, from
 * https://developer.nvidia.com/blog/how-query-device-properties-and-handle-errors-cuda-cc/
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

#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#ifdef __CUDACC__
#include <cstdio>
#include <cstdlib>
#include <cassert>

/***************************************************************************
 *  Timing with CUDA Events API
 ***************************************************************************/

/*
 * from https://stackoverflow.com/questions/6959213/timing-a-cuda-application-using-events/6977536#6977536
 */
class EventTimer {

public:
    EventTimer() : mStarted(false), mStopped(false) {
        cudaEventCreate(&mStart);
        cudaEventCreate(&mStop);
    }

    ~EventTimer() {
        cudaEventDestroy(mStart);
        cudaEventDestroy(mStop);
    }

    void start(cudaStream_t s = 0) {
        cudaEventRecord(mStart, s);
        mStarted = true;
        mStopped = false;
    }

    void stop(cudaStream_t s = 0)  {
        assert(mStarted);
        cudaEventRecord(mStop, s);
        mStarted = false;
        mStopped = true;
    }

    float elapsed() {
        assert(mStopped);
        if (!mStopped)
            return 0;

        cudaEventSynchronize(mStop);
        float elapsed = 0;
        cudaEventElapsedTime(&elapsed, mStart, mStop);
        return elapsed;
    }

private:
    bool mStarted, mStopped;
    cudaEvent_t mStart, mStop;
};

/***************************************************************************
 *  Error checking
 ***************************************************************************/

/*
 * from https://gist.github.com/ashwin/2652488
 */
#define cudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define cudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

/*
 * from https://berenger.eu/blog/cusparse-cccuda-sparse-matrix-examples-csr-bcsr-spmv-and-conversions/
 */
#define cudaSparseCheck( test ) __cudaSparseCheckCore((test), __FILE__, __LINE__)

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifndef NO_CUDA_CHECK_ERROR
    if ( cudaSuccess != err ) {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        abort();
    }
#endif
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifndef NO_CUDA_CHECK_ERROR
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err ) {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        abort();
    }

    /* More careful checking. However, this will affect performance.
       Comment away if needed. */
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err ) {
        fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        abort();
    }
#endif
}

//inline void __cuSparseCheckErr(cusparseStatus_t code,
//                               const char *file, int line) {
//#ifndef NO_CUDA_CHECK_ERROR
//    if (code != CUSPARSE_STATUS_SUCCESS) {
//        fprintf(stderr,"Cuda Error %d : %s %s %d\n",
//                code, cusparseGetErrorString(code), file, line);
//        abort();
//    }
//#endif
//}

/***************************************************************************
 *  Device management
 ***************************************************************************/

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
    cudaDeviceGetAttribute(&computeCap,
                           cudaDevAttrComputeCapabilityMajor, devId);
    return computeCap;
}

int get_max_shared_mem() {
    int devId, smemSize;
    cudaSafeCall( cudaGetDevice(&devId) );
    cudaDeviceGetAttribute(&smemSize,
                           cudaDevAttrMaxSharedMemoryPerBlock, devId);
    return smemSize;
}

int get_sm_count() {
    int devId, numProcs;
    cudaSafeCall( cudaGetDevice(&devId) );
    cudaDeviceGetAttribute(&numProcs,
                           cudaDevAttrMultiProcessorCount, devId);
    return numProcs;
}

/*
 * Slow and exhaustive property querying.
 */
void get_device_info()
{
    struct cudaDeviceProp capabilities{};

    cudaGetDeviceProperties (&capabilities, 0);

    printf("->CUDA Platform & Capabilities\n");
    printf("Name: %s\n", capabilities.name);
    printf("totalGlobalMem: %.2f MB\n",
           capabilities.totalGlobalMem/1024.0f/1024.0f);
    printf("sharedMemPerBlock: %.2f KB\n",
           capabilities.sharedMemPerBlock/1024.0f);
    printf("regsPerBlock (32 bits): %d\n",
           capabilities.regsPerBlock);
    printf("warpSize: %d\n",
           capabilities.warpSize);
    printf("memPitch: %.2f KB\n",
           capabilities.memPitch/1024.0f);
    printf("maxThreadsPerBlock: %d\n",
           capabilities.maxThreadsPerBlock);
    printf("maxThreadsDim: %d x %d x %d\n",
           capabilities.maxThreadsDim[0],
           capabilities.maxThreadsDim[1],
           capabilities.maxThreadsDim[2]);
    printf("maxGridSize: %d x %d\n",
           capabilities.maxGridSize[0],
           capabilities.maxGridSize[1]);
    printf("totalConstMem: %.2f KB\n",
           capabilities.totalConstMem/1024.0f);
    printf("major.minor: %d.%d\n",
           capabilities.major,
           capabilities.minor);
    printf("clockRate: %.2f MHz\n",
           (float)capabilities.clockRate/1024.0f);
    printf("textureAlignment: %zu\n",
           capabilities.textureAlignment);
    printf("deviceOverlap: %d\n",
           capabilities.deviceOverlap);
    printf("multiProcessorCount: %d\n",
           capabilities.multiProcessorCount);
}

#endif

#endif //CUDA_UTILS_H
