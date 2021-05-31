/****************************************************************************
 * @file gkernels.cu
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * Kernels for computing Betweenness centrality on a Nvidia GPUs.
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
 *
 ****************************************************************************/

#include "cl_kernels.cuh"

__global__ void get_closeness_p(double *cl,
                                const int *rows,
                                const int *cols,
                                int nvertices,
                                int nnz,
                                int *d,
                                int *next_source,
                                size_t pitch_d) {

    int tid = (int) threadIdx.x;
    if(tid >= max(2 * nnz, nvertices)) {
        return;
    }

    __shared__ int depth;
    __shared__ bool done;
    __shared__ int s;

    if (tid == 0) {
        s = (int) blockIdx.x;
    }

    int *d_row = (int *) ((char *) d + blockIdx.x * pitch_d);

    __syncthreads();

    /*
     * For each vertex...
     */
    while (s < nvertices) {

        if(tid == 0) {
            done = false;
            depth = 0;
        }
        __syncthreads();

        /*
         * Initialization.
         */
        for (int v = tid; v < nvertices; v += (int) blockDim.x) {
            if (v == s) {
                d_row[v] = 0;
            } else {
                d_row[v] = INT_MAX;
            }
        }

        // wait for all threads to complete the initial configuration
        __syncthreads();

        /*
         * Graph traversal for shortest path discovery and counting.
         */
        while(!done) {
            __syncthreads();
            done = true;
            __syncthreads();

            /*
             * For each edge...
             */
            for(int i = tid; i < nnz; i += (int) blockDim.x) {
                int v = rows[i];

                /*
                 * If the edge is incident to a vertex in the current frontier.
                 */
                if(d_row[v] == depth) {
                    int w = cols[i];

                    if(d_row[w] == INT_MAX) {
                        d_row[w] = d_row[v] + 1;
                        done = false;
                    }
                }
            }
            __syncthreads();

            if(tid == 0)
                depth++;
            __syncthreads();
        }

        __syncthreads();

        /*
         * Compute closeness centrality.
         */
        for (int i = (int) threadIdx.x; i < nvertices; i += (int) blockDim.x) {
            atomicAdd(&cl[i], (double) d_row[i]);
        }

        if (tid == 0) {
            s = atomicAdd(next_source, 1);
        }
        __syncthreads();
    }
}

void compute_cl_gpu_p(matrix_pcsr_t *g, double *cl, stats_t *stats) {

    double tstart, tend, first_tstart, last_tend;

    first_tstart = get_time();
    const unsigned int sm_count = get_sm_count();
    int next_source = (int) sm_count;
    int nnz = (g->row_offsets[g->nrows]);

    double *d_cl;
    int *d_rows, *d_cols, *d_dist, *d_next_source;
    size_t pitch_d;

    auto rows = (int *) malloc(nnz * sizeof(int));
    expand_row_pointer(g->nrows, g->row_offsets, rows);

    /*
     * Setup block and grid dimensions.
     */
    const unsigned int threads_per_block = get_max_threads_per_block();
    const unsigned int blocks_per_grid = sm_count;
    dim3 block = {threads_per_block, 1, 1};
    dim3 grid = {blocks_per_grid, 1, 1};

    /*
    * Load the CSR matrix on the device.
    */
    cudaSafeCall(cudaMalloc((void **) &d_rows,
                            nnz * sizeof(int)));
    cudaSafeCall(cudaMalloc((void **) &d_cols,
                            nnz * sizeof(int)));

    cudaSafeCall(cudaMemcpy(d_rows, rows,
                            nnz * sizeof(int),
                            cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy(d_cols, g->cols,
                            nnz * sizeof(int),
                            cudaMemcpyHostToDevice));

    /*
     * Load cl.
     */
    cudaSafeCall(cudaMalloc((void **) &d_cl, g->nrows * sizeof(double)));
    cudaSafeCall(cudaMemset(d_cl, 0, g->nrows * sizeof(double)));

    /*
     * Load auxiliary arrays for cl.
     */
    cudaSafeCall(cudaMallocPitch((void **) &d_dist, &pitch_d,
                                 g->nrows * sizeof(int), grid.x));

    /*
     * Load single-variables.
     */
    cudaSafeCall(cudaMalloc((void **) &d_next_source, sizeof(int)));
    cudaSafeCall(cudaMemcpy(d_next_source, &next_source,
                            sizeof(int),
                            cudaMemcpyHostToDevice));

    tend = get_time();
    stats->load_time = tend - first_tstart;

    /*
     * Execute the cl computation.
     */
    tstart = get_time();
    get_closeness_p<<<grid, block>>>(d_cl,
                                     d_rows,
                                     d_cols,
                                     g->nrows,
                                     nnz,
                                     d_dist,
                                     d_next_source,
                                     pitch_d);

    cudaCheckError();

    cudaSafeCall(cudaMemcpy(cl, d_cl,
                            g->nrows * sizeof(double),
                            cudaMemcpyDeviceToHost));

    /*
     * Finish computation of the closeness on the CPU.
     */
    for (int i = 0; i < g->nrows; i++) {
        cl[i] = ((double) g->nrows - 1.0) / cl[i];
    }

    tend = get_time();
    stats->bc_comp_time = tend - tstart;

    tstart = get_time();

    /*
     * Device resource deallocation.
     */
    cudaSafeCall(cudaFree(d_rows));
    cudaSafeCall(cudaFree(d_cols));
    cudaSafeCall(cudaFree(d_next_source));
    cudaSafeCall(cudaFree(d_cl));
    cudaSafeCall(cudaFree(d_dist));

    last_tend = get_time();
    stats->unload_time = last_tend - tstart;
    stats->total_time = last_tend - first_tstart;
}
