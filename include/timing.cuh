/****************************************************************************
 *
 * timing.cuh - Utility functions for NVIDIA GPUs profiling using the CUDA
 * Event API
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

#include <cassert>
#include "errcheck.cuh"

// outputs size in Mb
#define size_mb(size) ((size) / (1024 * 1024))

// outputs bandwidth in GB/s
#define bandwidth(bytes, microseconds) ((bytes) * 1e-6 / (microseconds))

/*
 * from https://stackoverflow.com/questions/6959213/timing-a-cuda-application-using-events/6977536#6977536
 */
class EventTimer {

public:
    EventTimer() : mStarted(false), mStopped(false) {
        cudaSafeCall( cudaEventCreate(&mStart) );
        cudaSafeCall( cudaEventCreate(&mStop) );
    }

    ~EventTimer() {
        cudaSafeCall( cudaEventDestroy(mStart) );
        cudaSafeCall( cudaEventDestroy(mStop) );
    }

    void start(cudaStream_t s = 0) {
        cudaSafeCall( cudaEventRecord(mStart, s) );
        mStarted = true;
        mStopped = false;
    }

    void stop(cudaStream_t s = 0)  {
        assert(mStarted);
        cudaSafeCall( cudaEventRecord(mStop, s) );
        mStarted = false;
        mStopped = true;
    }

    float elapsed() {
        assert(mStopped);
        if (!mStopped)
            return 0;

        cudaSafeCall( cudaEventSynchronize(mStop) );
        float elapsed = 0;
        cudaSafeCall( cudaEventElapsedTime(&elapsed, mStart, mStop) );
        return elapsed;
    }

    void report() {
        printf("Time elapsed: %f\n", this->elapsed());
    }

private:
    bool mStarted, mStopped;
    cudaEvent_t mStart, mStop;
};

/*
 * The throughput is measured as ...
 */
class ThroughputComputer {

};

/*
 * Only for profiling and eventually debugging memory transfers from device to
 * host and from host to device. |nvperf| profiles them too.
 */
float cudaMemcpyProfiled(void *dst, const void *src, size_t count,
                         cudaMemcpyKind kind) {

    float time;
    cudaEvent_t startEvent, stopEvent;
    cudaSafeCall( cudaEventCreate(&startEvent) );
    cudaSafeCall( cudaEventCreate(&stopEvent) );

    cudaSafeCall( cudaEventRecord(startEvent, 0) );
    cudaSafeCall( cudaMemcpy(dst, src, count, kind) );
    cudaSafeCall( cudaEventRecord(stopEvent, 0) );
    cudaSafeCall( cudaEventSynchronize(stopEvent) );

    cudaSafeCall( cudaEventElapsedTime(&time, startEvent, stopEvent) );

#ifdef CUDA_DEBUG
    printf("Transfer size (MB): %lu\n", size_mb(count));
    printf("Effective Bandwidth %s (GB/s): %f\n",
           (kind == cudaMemcpyHostToDevice) ? "Host to Device"
                                            : "Device to Host", bandwidth(count, time));
#endif

    return time;
}

#endif

#endif //CUDA_UTILS_H
