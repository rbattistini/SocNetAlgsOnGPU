/****************************************************************************
 *
 * errcheck.cuh - Utility functions for NVIDIA GPUs error checking
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

#pragma once
#ifndef ERRCHECK_CUH
#define ERRCHECK_CUH

#ifdef __CUDACC__

/*
 * from https://gist.github.com/ashwin/2652488
 */
#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cudaCheckError() __cudaCheckError(__FILE__, __LINE__)

/*
 * from https://berenger.eu/blog/cusparse-cccuda-sparse-matrix-examples-csr-bcsr-spmv-and-conversions/
 */
#define cudaSparseCheck(test) __cudaSparseCheckCore((test), __FILE__, __LINE__)

inline void __cudaSafeCall(cudaError err, const char *file, const int line) {
#ifdef CUDA_DEBUG
    if (cudaSuccess != err) {
        fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(err));
        abort();
    }
#endif
}

inline void __cudaCheckError(const char *file, const int line) {
#ifdef CUDA_DEBUG
    cudaError err = cudaGetLastError();
    if (cudaSuccess != err) {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(err));
        abort();
    }

    /* More careful checking. However, this will affect performance.
       Comment away if needed. */
    err = cudaDeviceSynchronize();
    if (cudaSuccess != err) {
        fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(err));
        abort();
    }
#endif
}

//inline void __cuSparseCheckErr(cusparseStatus_t code,
//                               const char *file, int line) {
//#ifdef CUDA_DEBUG
//    if (code != CUSPARSE_STATUS_SUCCESS) {
//        fprintf(stderr,"Cuda Error %d : %s %s %d\n",
//                code, cusparseGetErrorString(code), file, line);
//        abort();
//    }
//#endif
//}
#endif

#endif//ERRCHECK_CUH
