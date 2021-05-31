/****************************************************************************
 * @file bc_ep_kernel.cu
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Kernel for computing Betweenness centrality on a Nvidia GPU using
 * the edge parallel technique.
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

#include "bc_ep_kernel.cuh"

__global__ void get_vertex_betweenness_epp(double *bc,
                                           const int *rows,
                                           const int *cols,
                                           int nnz,
                                           int nvertices,
                                           int *d,
                                           unsigned long long *sigma,
                                           double *delta,
                                           int *next_source,
                                           size_t pitch_d,
                                           size_t pitch_sigma,
                                           size_t pitch_delta) {

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
    auto *delta_row = (double *) ((char *) delta + blockIdx.x * pitch_delta);
    auto *sigma_row = (unsigned long long *) ((char *) sigma +
                                              blockIdx.x * pitch_sigma);
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
                sigma_row[v] = 1;
            } else {
                d_row[v] = INT_MAX;
                sigma_row[v] = 0;
            }
            delta_row[v] = 0.0;
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

                    if(d_row[w] == (d_row[v] + 1)) {
                        atomicAdd(&sigma_row[w], sigma_row[v]);
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
         * Dependency accumulation by back-propagation.
         */
        while(depth > 1) {

            if(tid == 0)
                depth--;
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

                    if(d_row[w] == (d_row[v] + 1)) {
                        if(sigma_row[w] != 0) {
                            atomicAdd(&delta_row[v],
                                      (1.0f + delta_row[w]) *
                                      ((double) sigma_row[v] / (double) sigma_row[w]));
                        }
                    }
                }
            }
            __syncthreads();
        }

        /*
         * Compute betweenness centrality.
         */
        for (int i = tid; i < nvertices; i += (int) blockDim.x) {
            if(i != s)
                atomicAdd(&bc[i], delta_row[i]);
        }
        __syncthreads();

        if(tid == 0)
            s = atomicAdd(next_source, 1);
        __syncthreads();
    }
}

void compute_bc_gpu_epp(matrix_pcsr_t *g, double *bc, stats_t *stats) {

    double tstart, tend, first_tstart, last_tend;

    first_tstart = get_time();
    int nnz = (g->row_offsets[g->nrows]);

    unsigned long long *d_sigma;
    double *d_bc, *d_delta;
    int *d_rows, *d_cols, *d_dist, *d_next_source;
    size_t pitch_d, pitch_sigma, pitch_delta;

    auto rows = (int *) malloc(nnz * sizeof(int));
    expand_row_pointer(g->nrows, g->row_offsets, rows);

    /*
     * Setup block and grid dimensions.
     */
    const unsigned int sm_count = get_sm_count();
    int next_source = (int) sm_count;

    const unsigned int threads_per_block = get_max_threads_per_block();
    const unsigned int blocks_per_grid = sm_count;
    dim3 block = {threads_per_block, 1, 1};
    dim3 grid = {blocks_per_grid, 1, 1};

    /*
    * Load the COO matrix on the device.
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

    free(rows);

    /*
    * Load bc.
    */
    cudaSafeCall(cudaMalloc((void **) &d_bc, g->nrows * sizeof(double)));
    cudaSafeCall(cudaMemset(d_bc, 0, g->nrows * sizeof(double)));

    /*
     * Load auxiliary arrays for bc.
     */
    cudaSafeCall(cudaMallocPitch((void **) &d_dist, &pitch_d,
                                 g->nrows * sizeof(int), grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_sigma, &pitch_sigma,
                                 g->nrows * sizeof(unsigned long long),
                                 grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_delta, &pitch_delta,
                                 g->nrows * sizeof(double), grid.x));

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
     * Execute the bc computation.
     */
    tstart = get_time();
    get_vertex_betweenness_epp<<<grid, block>>>(d_bc,
                                                d_rows,
                                                d_cols,
                                                nnz,
                                                g->nrows,
                                                d_dist,
                                                d_sigma,
                                                d_delta,
                                                d_next_source,
                                                pitch_d,
                                                pitch_sigma,
                                                pitch_delta);

    cudaCheckError();

    cudaSafeCall(cudaMemcpy(bc, d_bc,
                            g->nrows * sizeof(double),
                            cudaMemcpyDeviceToHost));

    /*
     * Count each edge only one time.
     */
    for (int k = 0; k < g->nrows; k++)
        bc[k] /= 2;

    cudaSafeCall(cudaDeviceSynchronize());
    tend = get_time();
    stats->bc_comp_time = tend - tstart;
    tstart = get_time();

    /*
     * Device resource deallocation.
     */
    cudaSafeCall(cudaFree(d_rows));
    cudaSafeCall(cudaFree(d_next_source));
    cudaSafeCall(cudaFree(d_cols));
    cudaSafeCall(cudaFree(d_bc));
    cudaSafeCall(cudaFree(d_sigma));
    cudaSafeCall(cudaFree(d_dist));
    cudaSafeCall(cudaFree(d_delta));

    last_tend = get_time();
    stats->unload_time = last_tend - tstart;
    stats->total_time = last_tend - first_tstart;
}
