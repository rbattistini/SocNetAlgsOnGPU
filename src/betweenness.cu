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

//#define BENCHMARK

#include <cstdlib>
#include <cstdio>

#include "../include/matio.h"
#include "../include/matstorage.h"
//#include "timing.cuh"
//#include "errcheck.cuh"
//#include "graphs_kernels.cuh"
#include "../include/device_props.cuh"

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
 * TODO handle output file
 *
 * TODO add measurement of throughput as edges traversed per second (TEPS)
 * refer to: https://gunrock.github.io/docs/#/gunrock/methodology
 */
int main( int argc, char *argv[] ) {

    device_graph_t graph;
    matrix_coo_t m_coo;
    matrix_csr_t m_csr;
    int *d_row_offsets, *d_cols;
    float *bc_cpu, *bc_gpu;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s [input_filename] [output_filename]", argv[0]);
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

#ifdef BENCHMARK
    /*
     * Compute BC with the algorithm that uses  multithreading of the BGL.
     */
    size_t nvertices = m_csr.nrows;
    bc_cpu = (float*) malloc(nvertices * sizeof(bc_cpu));

    /*
     * Check whether BC was computed correctly.
     */
    check_bc(m_csr, bc_cpu, bc_gpu);
    free(bc_cpu);
#endif

    free_matrix_coo(&m_coo);
    free_matrix_csr(&m_csr);
//    cudaSafeCall( cudaFree(d_row_offsets) );
//    cudaSafeCall( cudaFree(d_cols) );

    return EXIT_SUCCESS;
}
