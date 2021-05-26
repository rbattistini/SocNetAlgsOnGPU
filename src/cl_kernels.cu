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

__global__ void get_closeness_p(float *cl,
                                const int *row_offsets,
                                const int *cols,
                                int nvertices,
                                int *d,
                                int *curr_queue,
                                int *next_queue,
                                int *next_source,
                                size_t pitch_d,
                                size_t pitch_qcurr,
                                size_t pitch_qnext) {
    __shared__ int ind;
    __shared__ int s;

    int tid = (int) threadIdx.x;

    int *d_row = (int *) ((char *) d + blockIdx.x * pitch_d);

    __shared__ int *qcurr_row;
    __shared__ int *qnext_row;

    if (tid == 0) {
        ind = (int) blockIdx.x;
        s = ind;
        qcurr_row = (int *) ((char *) curr_queue + blockIdx.x * pitch_qcurr);
        qnext_row = (int *) ((char *) next_queue + blockIdx.x * pitch_qnext);
    }

    __syncthreads();

    /*
     * For each vertex...
     */
    while (ind < nvertices) {

        for (int k = (int) threadIdx.x; k < nvertices; k += (int) blockDim.x) {
            if (k == s) {
                d_row[k] = 0;
            } else {
                d_row[k] = INT_MAX;
            }
        }

        __syncthreads();

        __shared__ int qcurr_len;
        __shared__ int qnext_len;
        __shared__ int depth;
        __shared__ bool done;

        if (tid == 0) {
            qcurr_row[0] = s;
            qcurr_len = 1;
            qnext_len = 0;
            depth = 0;
            done = false;
        }

        __syncthreads();

        //Do first iteration separately since we already know the edges to traverse
        int start = (int) threadIdx.x + row_offsets[s];
        int end = row_offsets[s + 1];
        for (int r = start; r < end; r += (int) blockDim.x) {
            int w = cols[r];

            if (d_row[w] == INT_MAX) {
                d_row[w] = 1;
                int t = atomicAdd(&qnext_len, 1);
                qnext_row[t] = w;
            }
        }
        __syncthreads();

        /*
         * If the next frontier is empty the traversal is complete.
         */
        if (qnext_len == 0) {
            done = true;
        } else {
            /*
             * Update the current frontier with elements from the next.
             * Update the stack with the element that will be visited in
             * the current frontier.
             */
            for (int i = (int) threadIdx.x; i < qnext_len; i += (int) blockDim.x) {
                qcurr_row[i] = qnext_row[i];
            }
            __syncthreads();

            if (tid == 0) {

                /*
                 * Update queues' length.
                 */
                qcurr_len = qnext_len;
                qnext_len = 0;

                /*
                 * Keep track of the depth from which the dependency
                 * accumulation will begin.
                 */
                depth += 1;
            }
        }
        __syncthreads();

        /*
         * Bfs. One iteration for each depth of the graph.
         */
        while (!done) {

            __shared__ int next_index;
            if (tid == 0) {
                next_index = blockDim.x;
            }
            __syncthreads();

            int k = threadIdx.x;

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
                }

                k = atomicAdd(&next_index, 1);
            }

            __syncthreads();

            /*
             * If the next frontier is empty the traversal is complete.
             */
            if (qnext_len == 0) {
                done = true;
            } else {
                /*
                 * Update the current frontier with elements from the next.
                 * Update the stack with the element that will be visited in
                 * the current frontier.
                 */
                for (int i = (int) threadIdx.x; i < qnext_len;
                     i += (int) blockDim.x) {
                    qcurr_row[i] = qnext_row[i];
                }
                __syncthreads();

                if (tid == 0) {

                    /*
                     * Update queues' length.
                     */
                    qcurr_len = qnext_len;
                    qnext_len = 0;

                    /*
                     * Keep track of the depth from which the dependency
                     * accumulation will begin.
                     */
                    depth += 1;
                }
            }
            __syncthreads();
        }

        /*
         * Compute closeness centrality.
         */
        for (int i = (int) threadIdx.x; i < nvertices; i += (int) blockDim.x) {
            atomicAdd(&cl[i], (float) d_row[i]);
        }

        if (tid == 0) {
            ind = atomicAdd(next_source, 1);
            s = ind;
        }
        __syncthreads();
    }
}

void compute_cl_gpu_p(matrix_pcsr_t *g, float *cl, stats_t *stats) {

    double tstart, tend, first_tstart, last_tend;
    int n = stats->nrun;

    first_tstart = get_time();
    const unsigned int sm_count = get_sm_count();
    int next_source = (int) sm_count;

    float *d_cl;
    int *d_row_offsets,
            *d_cols,
            *d_dist,
            *d_qcurr,
            *d_qnext,
            *d_next_source;
    size_t pitch_d,
            pitch_qcurr,
            pitch_qnext;

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
     * Load cl.
     */
    cudaSafeCall(cudaMalloc((void **) &d_cl, g->nrows * sizeof(float)));
    cudaSafeCall(cudaMemset(d_cl, 0, g->nrows * sizeof(float)));

    /*
     * Load auxiliary arrays for cl.
     */
    cudaSafeCall(cudaMallocPitch((void **) &d_dist, &pitch_d,
                                 g->nrows * sizeof(int), grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_qcurr, &pitch_qcurr,
                                 g->nrows * sizeof(int), grid.x));
    cudaSafeCall(cudaMallocPitch((void **) &d_qnext, &pitch_qnext,
                                 g->nrows * sizeof(int), grid.x));

    /*
     * Load single-variables.
     */
    cudaSafeCall(cudaMalloc((void **) &d_next_source, sizeof(int)));
    cudaSafeCall(cudaMemcpy(d_next_source, &next_source,
                            sizeof(int),
                            cudaMemcpyHostToDevice));

    tend = get_time();
    stats->load_time[n] = tend - first_tstart;

    /*
     * Execute the cl computation.
     */
    tstart = get_time();
    get_closeness_p<<<grid, block>>>(d_cl,
                                     d_row_offsets,
                                     d_cols,
                                     g->nrows,
                                     d_dist,
                                     d_qcurr,
                                     d_qnext,
                                     d_next_source,
                                     pitch_d,
                                     pitch_qcurr,
                                     pitch_qnext);

    cudaCheckError();

    cudaSafeCall(cudaMemcpy(cl, d_cl,
                            g->nrows * sizeof(float),
                            cudaMemcpyDeviceToHost));

    /*
     * Finish computation of the closeness on the CPU.
     */
    for (int i = 0; i < g->nrows; i++) {
        cl[i] = ((float) g->nrows - 1.0f) / cl[i];
    }

    tend = get_time();
    stats->bc_comp_time[n] = tend - tstart;

    tstart = get_time();

    /*
     * Device resource deallocation.
     */
    cudaSafeCall(cudaFree(d_row_offsets));
    cudaSafeCall(cudaFree(d_cols));
    cudaSafeCall(cudaFree(d_next_source));
    cudaSafeCall(cudaFree(d_cl));
    cudaSafeCall(cudaFree(d_dist));
    cudaSafeCall(cudaFree(d_qcurr));
    cudaSafeCall(cudaFree(d_qnext));

    last_tend = get_time();
    stats->unload_time[n] = last_tend - tstart;
    stats->total_time[n] = last_tend - first_tstart;
}
