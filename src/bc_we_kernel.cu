/****************************************************************************
 * @file bc_we_kernel.cu
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Kernel for computing Betweenness centrality on a Nvidia GPU using
 * the work efficient technique.
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

#include "bc_we_kernel.cuh"

/**
 * Updates auxiliary data structures used for both bfs and the dependency
 * accumulation step.
 *
 * @param[out] qcurr_len
 * @param[in, out] qnext_len
 * @param[out] qcurr
 * @param[in] qnext
 * @param[out] stack_len
 * @param[out] stack
 * @param[out] ends
 * @param[in, out] ends_len
 * @param[out] depth
 * @param[in] tid
 */
__device__ void bfs_update_ds_wpitched(int *qcurr_len,
                                       int *qcurr,
                                       int *qnext_len,
                                       const int *qnext,
                                       int *stack_len,
                                       int *stack,
                                       int *ends_len,
                                       int *ends,
                                       int *depth,
                                       int tid) {

    /*
     * Update the current frontier with elements from the next.
     * Update the stack with the element that will be visited in
     * the current frontier.
     */
    for (int i = tid; i < *qnext_len; i += (int) blockDim.x) {
        qcurr[i] = qnext[i];
        stack[i + *stack_len] = qnext[i];
    }
    __syncthreads();

    if (tid == 0) {

        /*
         * Update stack depth annotation and stack length.
         */
        ends[*ends_len] =
                ends[*ends_len - 1] + *qnext_len;
        *ends_len += 1;

        *stack_len += *qnext_len;

        /*
         * Update queues' length.
         */
        *qcurr_len = *qnext_len;
        *qnext_len = 0;

        /*
         * Keep track of the depth from which the dependency
         * accumulation will begin.
         */
        *depth += 1;
    }
}

__global__ void get_vertex_betweenness_wep(double *bc,
                                           const int *row_offsets,
                                           const int *cols,
                                           int nvertices,
                                           int *d,
                                           unsigned long long *sigma,
                                           double *delta,
                                           int *curr_queue,
                                           int *next_queue,
                                           int *stack,
                                           int *endpoints,
                                           int *next_source,
                                           size_t pitch_d,
                                           size_t pitch_sigma,
                                           size_t pitch_delta,
                                           size_t pitch_qcurr,
                                           size_t pitch_qnext,
                                           size_t pitch_stack,
                                           size_t pitch_endpoints) {
    __shared__ int s;

    int tid = (int) threadIdx.x;

    int *d_row = (int *) ((char *) d + blockIdx.x * pitch_d);
    auto *delta_row = (double *) ((char *) delta + blockIdx.x * pitch_delta);
    auto *sigma_row = (unsigned long long *) ((char *) sigma +
                                              blockIdx.x * pitch_sigma);

    __shared__ int *qcurr_row;
    __shared__ int *qnext_row;
    __shared__ int *stack_row;
    __shared__ int *ends_row;

    if (tid == 0) {
        s = (int) blockIdx.x;
        qcurr_row = (int *) ((char *) curr_queue + blockIdx.x * pitch_qcurr);
        qnext_row = (int *) ((char *) next_queue + blockIdx.x * pitch_qnext);
        stack_row = (int *) ((char *) stack + blockIdx.x * pitch_stack);
        ends_row = (int *) ((char *) endpoints +
                            blockIdx.x * pitch_endpoints);
    }

    __syncthreads();

    /*
     * For each vertex...
     */
    while (s < nvertices) {
        for (int k = tid; k < nvertices; k += (int) blockDim.x) {
            if (k == s) {
                d_row[k] = 0;
                sigma_row[k] = 1;
            } else {
                d_row[k] = INT_MAX;
                sigma_row[k] = 0;
            }
            delta_row[k] = 0.0;
        }

        __syncthreads();

        __shared__ int qcurr_len;
        __shared__ int qnext_len;
        __shared__ int stack_len;
        __shared__ int depth;
        __shared__ int ends_len;
        __shared__ bool done;

        if (tid == 0) {
            qcurr_row[0] = s;
            qcurr_len = 1;
            qnext_len = 0;
            stack_row[0] = s;
            stack_len = 1;
            ends_row[0] = 0;
            ends_row[1] = 1;
            ends_len = 2;
            depth = 0;
            done = false;
        }

        __syncthreads();

        /*
         * Graph traversal for shortest path discovery and counting.
         */
        while (!done) {

            __shared__ int next_index;
            if (tid == 0) {
                next_index = (int) blockDim.x;
            }
            __syncthreads();

            int k = tid;

            /*
             * For each vertex in the current frontier.
             *
             * Assigns a thread to each vertex such that edges from other
             * portions of the graph are not unnecessarily traversed.
             */
            while (k < qcurr_len) {
                int v = qcurr_row[k];

                /*
                 * Add the neighbours of the vertex of the current frontier
                 * to the queue of the vertices of the next frontier.
                 */
                for (int r = row_offsets[v]; r < row_offsets[v + 1]; r++) {
                    int w = cols[r];

                    if (atomicCAS(&d_row[w], INT_MAX, d_row[v] + 1) ==
                        INT_MAX) {
                        int t = atomicAdd(&qnext_len, 1);
                        qnext_row[t] = w;
                    }

                    if (d_row[w] == (d_row[v] + 1)) {
                        atomicAdd(&sigma_row[w], sigma_row[v]);
                    }
                }

                k = atomicAdd(&next_index, 1);
            }

            __syncthreads();

            /*
             * If the next frontier is empty the traversal is complete.
             */
            if (qnext_len == 0) {
                break;
            } else {
                bfs_update_ds_wpitched(&qcurr_len, qcurr_row,
                                       &qnext_len, qnext_row,
                                       &stack_len, stack_row,
                                       &ends_len, ends_row,
                                       &depth,
                                       tid);
            }
            __syncthreads();
        }

        if (tid == 0) {
            depth = d_row[stack_row[stack_len - 1]] - 1;
        }
        __syncthreads();

        /*
         * Dependency accumulation by back-propagation.
         */
        while (depth > 0) {

            int start = tid + ends_row[depth];
            int end = ends_row[depth + 1];
            for (int i = start; i < end; i += (int) blockDim.x) {
                int w = stack_row[i];
                float dsw = 0;
                auto sw = (float) sigma_row[w];
                for (int z = row_offsets[w]; z < row_offsets[w + 1]; z++) {
                    int v = cols[z];
                    if (d_row[v] == (d_row[w] + 1)) {
                        dsw += (sw / sigma_row[v]) * (1.0f + delta_row[v]);
                    }
                }
                delta_row[w] = dsw;
            }

            __syncthreads();
            if (tid == 0) {
                depth--;
            }
            __syncthreads();
        }

        /*
         * Compute betweenness centrality.
         */
        for (int i = tid; i < nvertices; i += (int) blockDim.x) {
            atomicAdd(&bc[i], delta_row[i]);
        }

        if (tid == 0) {
            s = atomicAdd(next_source, 1);
        }
        __syncthreads();
    }
}

void compute_bc_gpu_wep(matrix_pcsr_t *g, double *bc, stats_t *stats) {

    double tstart, tend, first_tstart, last_tend;

    first_tstart = get_time();
    const unsigned int sm_count = get_sm_count();
    int next_source = (int) sm_count;

    unsigned long long *d_sigma;
    double *d_bc, *d_delta;
    int *d_row_offsets,
            *d_cols,
            *d_dist,
            *d_qcurr,
            *d_qnext,
            *d_stack,
            *d_ends,
            *d_next_source;
    size_t pitch_d,
            pitch_sigma,
            pitch_delta,
            pitch_stack,
            pitch_qcurr,
            pitch_qnext,
            pitch_ends;

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
    int nnz = (g->row_offsets[g->nrows]);

    cudaSafeCall(cudaMalloc((void **) &d_row_offsets,
                            (g->nrows + 1) * sizeof(int)));
    cudaSafeCall(cudaMalloc((void **) &d_cols,
                            nnz * sizeof(int)));

    cudaSafeCall(cudaMemcpy(d_row_offsets, g->row_offsets,
                            (g->nrows + 1) * sizeof(int),
                            cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy(d_cols, g->cols,
                            nnz * sizeof(int),
                            cudaMemcpyHostToDevice));

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
    cudaSafeCall(cudaMallocPitch((void **) &d_stack, &pitch_stack,
                                 g->nrows * sizeof(int), grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_qcurr, &pitch_qcurr,
                                 g->nrows * sizeof(int), grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_qnext, &pitch_qnext,
                                 g->nrows * sizeof(int), grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_ends, &pitch_ends,
                                 (g->nrows + 1) * sizeof(int), grid.x));

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
    get_vertex_betweenness_wep<<<grid, block>>>(d_bc,
                                                d_row_offsets,
                                                d_cols,
                                                g->nrows,
                                                d_dist,
                                                d_sigma,
                                                d_delta,
                                                d_qcurr,
                                                d_qnext,
                                                d_stack,
                                                d_ends,
                                                d_next_source,
                                                pitch_d,
                                                pitch_sigma,
                                                pitch_delta,
                                                pitch_qcurr,
                                                pitch_qnext,
                                                pitch_stack,
                                                pitch_ends);

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
    cudaSafeCall(cudaFree(d_row_offsets));
    cudaSafeCall(cudaFree(d_cols));
    cudaSafeCall(cudaFree(d_next_source));
    cudaSafeCall(cudaFree(d_bc));
    cudaSafeCall(cudaFree(d_sigma));
    cudaSafeCall(cudaFree(d_dist));
    cudaSafeCall(cudaFree(d_delta));
    cudaSafeCall(cudaFree(d_stack));
    cudaSafeCall(cudaFree(d_qcurr));
    cudaSafeCall(cudaFree(d_qnext));
    cudaSafeCall(cudaFree(d_ends));

    last_tend = get_time();
    stats->unload_time = last_tend - tstart;
    stats->total_time = last_tend - first_tstart;
}
