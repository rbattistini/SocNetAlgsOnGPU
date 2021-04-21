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

#define CUDA_DEBUG
#include <cstdlib>
#include <cstdio>

#include "matio.h"
#include "matstorage.h"
//#include "timing.cuh"
//#include "errcheck.cuh"
#include "device_props.cuh"

typedef struct device_graph_t {
    int *row_offsets;
    int *cols;
    int *n;
    int *m;
} device_graph_t;

/*
 * TODO memory requirement
 * TODO diameter estimation
 * TODO sampling technique
 */
int main( int argc, char *argv[] ) {

    device_graph_t graph;
    matrix_coo_t m_coo;
    matrix_csr_t m_csr;
    int *d_row_offsets, *d_cols;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s [input_filename]", argv[0]);
        return EXIT_FAILURE;
    }

    if(get_compute_capability() != 6) {
        fprintf(stderr, "This program is meant to be executed only if compute"
                     " capability is 6.x\n");
        return EXIT_FAILURE;
    } else {
        print_gpu_overview();
    }

    /*
     * Coarse-grained parallelism.
     */
    const unsigned int threadsPerBlock = get_max_threads_per_block();
    const unsigned int blocksPerGrid = get_sm_count();
    dim3 block = {threadsPerBlock, 1, 1};
    dim3 grid = {blocksPerGrid, 1, 1};

    /*
     * Load Matrix Market matrix stored in .mm format as a COO matrix.
     */
    read_matrix(argv[1], &m_coo);

    /*
     * Convert the internal storage representation of the matrix from COO to
     * the more efficient CSR.
     */
    coo_to_csr(&m_coo, &m_csr);

    /*
     * Allocate the matrix in CSR on the device.
     */
//    size_t size_rows = m_coo.nnz;
//    size_t size_cols = m_coo.nnz;
//    cudaSafeCall( cudaMalloc((void**)&d_row_offsets, size_rows) );
//    cudaSafeCall( cudaMalloc((void**)&d_cols, size_cols) );

    free_matrix_coo(&m_coo);
    free_matrix_csr(&m_csr);
//    cudaSafeCall( cudaFree(d_row_offsets) );
//    cudaSafeCall( cudaFree(d_cols) );

    return EXIT_SUCCESS;
}
//
//__global__ void vertex_parallel_bfs(int nnz, int *sigma, int *d,
//                                       int source, const int *row_offsets,
//                                       const int *cols) {
//
//    int idx = (int) threadIdx.x;
//
//    /*
//     * Initialize d and sigma.
//     */
//    for(int v = idx; v < nnz; v += (int) blockDim.x) {
//        if(v == source) {
//            d[v] = 0;
//            sigma[v] = 1;
//        } else {
//            d[v] = -1;
//            sigma[v] = 0;
//        }
//    }
//    __shared__ int current_depth;
//    __shared__ bool done;
//
//    if(idx == 0) {
//        done = false;
//        current_depth = 0;
//    }
//    __syncthreads();
//
//    /*
//     * Calculate the number of shortest paths and the distance from s
//     * (the root) to each vertex.
//     */
//    while(!done)
//    {
//        __syncthreads();
//        done = true;
//        __syncthreads();
//
//        /*
//         * For each vertex, traverse its neighbours and update the
//         */
//        for(int v = idx; v < nnz; v += (int) blockDim.x) {
//
//            if(d[v] == current_depth) {
//
//                for(int r = row_offsets[v]; r < row_offsets[v+1]; r++) {
//
//                    int w = cols[r];
//
//                    if(d[w] == -1) {
//                        d[w] = d[v] + 1;
//                        done = false;
//                    }
//
//                    if(d[w] == (d[v] + 1))
//                        atomicAdd(&sigma[w], sigma[v]);
//                }
//            }
//        }
//        __syncthreads();
//
//        if(idx == 0)
//            current_depth++;
//    }
//}
//
//__global__ void edge_parallel_bfs(int nnz, int *sigma, int *d,
//                                     int source, const int *rows,
//                                     const int *cols, int nedges) {
//
//    int idx = (int) threadIdx.x;
//
//    /*
//     * Initialize d and sigma.
//     */
//    for(int k = idx; k < nnz; k += (int) blockDim.x)
//    {
//        if(k == source) {
//            d[k] = 0;
//            sigma[k] = 1;
//        } else {
//            d[k] = -1;
//            sigma[k] = 0;
//        }
//    }
//
//    __shared__ int current_depth;
//    __shared__ bool done;
//
//    if(idx == 0) {
//        done = false;
//        current_depth = 0;
//    }
//    __syncthreads();
//
//    /*
//     * Calculate the number of shortest paths and the distance from s
//     * (the root) to each vertex.
//     */
//    while(!done)
//    {
//        __syncthreads();
//        done = true;
//        __syncthreads();
//
//        for(int k = idx; k < nedges; k += (int) blockDim.x) {
//
//            int v = rows[k];
//
//            // If the head is in the vertex frontier, look at the tail
//            if(d[v] == current_depth) {
//
//                int w = cols[k];
//
//                if(d[w] == -1) {
//                    d[w] = d[v] + 1;
//                    done = false;
//                }
//
//                if(d[w] == (d[v] + 1))
//                    atomicAdd(&sigma[w], sigma[v]);
//            }
//        }
//
//        __syncthreads();
//        if (threadIdx.x == 0)
//            current_depth++;
//    }
//}
//
///*
// * REVIEW
// */
//__global__ void work_efficient_bfs(int nnz, int *sigma, int *d,
//                                      int source, const int *row_offsets,
//                                      const int *cols) {
//
//    int idx = (int) threadIdx.x;
//
//    /*
//     * Initialize d and sigma.
//     */
//    for(int k = idx; k < nnz; k += (int) blockDim.x)
//    {
//        if(k == source) {
//            d[k] = 0;
//            sigma[k] = 1;
//        } else {
//            d[k] = -1;
//            sigma[k] = 0;
//        }
//    }
//
//    __shared__ int Q_len;
//    __shared__ int Q2_len;
//
//    if(idx == 0)
//    {
//        Q[0] = source;
//        Q_len = 1;
//        Q2_len = 0;
//    }
//    __syncthreads();
//
//    while(true)
//    {
//        for(int k = idx; k < Q_len; k += (int) blockDim.x) {
//
//            int v = Q[k];
//
//            for(int r = row_offsets[v]; r < row_offsets[v+1]; r++) {
//
//                // Use atomicCAS to prevent duplicates
//                if(atomicCAS(&d[w], -1, d[v] + 1) == -1) {
//                    int t = atomicAdd(&Q2_len,1);
//                    Q2[t] = w;
//                }
//
//                if(d[w] == (d[v]+1))
//                    atomicAdd(&sigma[w],sigma[v]);
//            }
//        }
//        __syncthreads();
//
//        /*
//         * The next vertex frontier is empty, so we're done searching.
//         */
//        if(Q2_len == 0) {
//            break;
//        } else {
//
//            for(int k = idx; k < Q2_len; k += (int) blockDim.x) {
//                Q[k] = Q2[k];
//            }
//
//            __syncthreads();
//
//            if(idx == 0) {
//                Q_len = Q2_len;
//                Q2_len = 0;
//            }
//            __syncthreads();
//        }
//    }
//}
