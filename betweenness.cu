/****************************************************************************
 *
 * betweenness.cpp - Serial algorithm for computing betweenness centrality
 *
 * Based on ...
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
 * ---------------------------------------------------------------------------
 *
 * Compile with:
 * nvcc++ --gpu-architecture compute_60 -O3 betweenness.cpp -o betweenness
 * --linker-options "-lcusparse"
 *
 * Run with:
 * ./betweenness [input_filename]
 *
 ****************************************************************************/

#include <cstdlib>
#include <cstdio>
#include "utils.h"
#include "cuda_utils.cuh"

/*
 * Cuda Event API and Profiler API and nvperf tool for performance evaluation
 * from https://developer.nvidia.com/blog/how-implement-performance-metrics-cuda-cc/
 *
 * Cuda Graph for optimizations
 *
 * Both based on Cuda Streams, used for optimizing memory transfers
 * - https://developer.nvidia.com/blog/how-optimize-data-transfers-cuda-cc/
 * - https://developer.nvidia.com/blog/how-overlap-data-transfers-cuda-cc/
 * - https://developer.nvidia.com/blog/beyond-gpu-memory-limits-unified-memory-pascal/
 *
 * Memory management and Memory Access Patterns
 * - https://developer.nvidia.com/blog/unified-memory-cuda-beginners/
 * - https://developer.nvidia.com/blog/how-access-global-memory-efficiently-cuda-c-kernels/
 * - https://developer.nvidia.com/blog/using-shared-memory-cuda-cc/
 *
 * Warp-level and block-level synchronization
 * - https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/
 * - https://developer.nvidia.com/blog/cooperative-groups/
 *
 * Handling structures
 * - https://stackoverflow.com/questions/12778949/cuda-memory-alignment/12779757#12779757
 * - https://stackoverflow.com/questions/9309195/copying-a-struct-containing-pointers-to-cuda-device/9323898#9323898
 */
int main( int argc, char *argv[] ) {

    matrix_coo_t m_coo;
    matrix_csr_t m_csr;
    int *d_row_offsets, *d_cols, nedges;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s [input_filename] [nedges]", argv[0]);
        return EXIT_FAILURE;
    }

    if(get_compute_capability() != 6) {
        fprintf(stderr, "This program is meant to be executed only if compute"
                     " capability is 6.x\n");
        return EXIT_FAILURE;
    }

    /*
     * Load Matrix Market matrix stored in .mm format as a COO matrix.
     */
    nedges = (int) strtol(argv[2], nullptr, 10);
    loadMatrixMarketFile(argv[1], &m_coo, nedges);

    /*
     * Convert the internal storage representation of the matrix from COO to
     * the more efficient CSR.
     */
    coo2csr(&m_coo, &m_csr);

//    printf("\n%d %d %d\n", m_csr.nrows, m_csr.nrows, m_csr.nnz);
//    printArray(m_csr.cols, m_csr.nnz);
//    printArray(m_csr.row_offsets, m_csr.nrows + 1);

    /*
     * Allocate the matrix in CSR on the device.
     */
    size_t size_rows = m_coo.nnz;
    size_t size_cols = m_coo.nnz;
//    cudaSafeCall( cudaMalloc((void**)&d_row_offsets, size_rows) );
//    cudaSafeCall( cudaMalloc((void**)&d_cols, size_cols) );

    /*
     * TODO
     */

    free_matrix_coo(&m_coo);
    free_matrix_csr(&m_csr);
//    cudaSafeCall( cudaFree(d_row_offsets) );
//    cudaSafeCall( cudaFree(d_cols) );

    return EXIT_SUCCESS;
}
