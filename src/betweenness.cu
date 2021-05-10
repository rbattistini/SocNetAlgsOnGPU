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
 ****************************************************************************/

//#define BENCHMARK
// TODO bandwidth measure to be tested

#include <cstdio>
#include <cstdlib>
#include "matio.h"
#include "matstorage.h"
#include "device_props.cuh"
#include "timing.cuh"
#include "errcheck.cuh"
#include "gkernels.cuh"

int main(int argc, char *argv[]) {

    matrix_pcoo_t m_coo;
    matrix_pcsr_t m_csr;
    gprops_t gp;
    int *d_row_offsets, *d_cols, *d_dist, *d_sigma;
    float *bc_gpu, *d_bc, *d_delta;
    EventTimer chrono;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s [input_filename] [output_filename]", argv[0]);
        return EXIT_FAILURE;
    }

    if (get_compute_capability() != 6) {
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
//    dim3 block = {threadsPerBlock, 1, 1};
//    dim3 grid = {blocksPerGrid, 1, 1};

    /*
     * Load Matrix Market matrix stored in .mm format as a COO matrix.
     */
    read_matrix(argv[1], &m_coo, &gp);

    /*
     * Convert the internal storage representation of the matrix from COO to
     * the more efficient CSR.
     */
    pcoo_to_pcsr(&m_coo, &m_csr);

    print_matrix_csr(&m_csr);
    printf("\n");

    /*
     * Allocate the matrix in CSR on the device.
     */
    size_t nnz = m_csr.row_offsets[m_csr.nrows];
    cudaSafeCall( cudaMalloc((void**)&d_row_offsets,
                             (m_csr.nrows + 1) * sizeof(int)) );
    cudaSafeCall( cudaMalloc((void**)&d_cols, nnz * sizeof(int)) );
    cudaSafeCall( cudaMalloc((void**)&d_dist, m_csr.nrows * sizeof(int)) );
    cudaSafeCall( cudaMalloc((void**)&d_sigma, m_csr.nrows * sizeof(int)) );
    cudaSafeCall( cudaMalloc((void**)&d_delta, m_csr.nrows * sizeof(float)) );
    cudaSafeCall( cudaMalloc((void**)&d_bc, m_csr.nrows * sizeof(float)) );

    /*
     * Allocate memory for the bc scores on the host.
     */
    bc_gpu = (float*) malloc(m_csr.nrows * sizeof(float));

    /*
     * Compute bc.
     */
    chrono.start();
    for(int i = 0; i < m_csr.nrows; i++) {

        int s = i;
        vtx_par_bfs<<<1, threadsPerBlock>>>(s,
                                                        d_dist,
                                                        d_sigma,
                                                        m_csr.nrows,
                                                        nnz,
                                                        d_row_offsets,
                                                        d_cols);

        cudaCheckError();

        vtx_par_dep_acc<<<1, threadsPerBlock>>>(s,
                                                            d_dist,
                                                            d_sigma,
                                                            d_delta,
                                                            d_bc,
                                                            m_csr.nrows,
                                                            nnz,
                                                            d_row_offsets,
                                                            d_cols);

        cudaCheckError();
    }

    cudaSafeCall( cudaMemcpy(d_bc, bc_gpu, m_csr.nrows * sizeof(float),
                             cudaMemcpyDeviceToHost));
    chrono.stop();

    /*
     * Report time elapsed and throughput.
     */
//    printf("Time elapsed: %.2f\n", chrono.elapsed());
    print_array(bc_gpu, m_csr.nrows);

#ifdef BENCHMARK
    /*
     * Compute BC with the algorithm that uses  multithreading of the BGL.
     */
    size_t nvertices = m_csr.nrows;
    bc_cpu = (float *) malloc(nvertices * sizeof(bc_cpu));

    /*
     * Check whether BC was computed correctly.
     */
    check_bc(m_csr, bc_cpu, bc_gpu);
    free(bc_cpu);
#endif

    /*
     * Cleanup.
     */
    free_matrix_pcoo(&m_coo);
    free_matrix_pcsr(&m_csr);
    cudaSafeCall( cudaFree(d_row_offsets) );
    cudaSafeCall( cudaFree(d_cols) );
    cudaSafeCall( cudaFree(d_bc) );
    cudaSafeCall( cudaFree(d_delta) );
    cudaSafeCall( cudaFree(d_sigma) );
    cudaSafeCall( cudaFree(d_dist) );

    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);

    return EXIT_SUCCESS;
}
