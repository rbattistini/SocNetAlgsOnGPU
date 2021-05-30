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

#include "bc_we_kernels.cuh"

__device__ void bfs_update_ds(int *qcurr_len,
                              int *qcurr,
                              int *qnext_len,
                              const int *qnext,
                              int *stack_len,
                              int *stack,
                              int *ends_len,
                              int *ends,
                              int *depth,
                              int tid,
                              int nvertices) {

    /*
     * Update the current frontier with elements from the next.
     * Update the stack with the element that will be visited in
     * the current frontier.
     */
    for (int i = (int) threadIdx.x; i < *qnext_len; i += (int) blockDim.x) {
        qcurr[blockIdx.x * nvertices + i] = qnext[blockIdx.x * nvertices + i];
        stack[blockIdx.x * nvertices + i + *stack_len] = qnext[blockIdx.x * nvertices + i];
    }
    __syncthreads();

    if (tid == 0) {

        /*
         * Update stack depth annotation and stack length.
         */
        ends[blockIdx.x * nvertices + *ends_len] =
                ends[blockIdx.x * nvertices + *ends_len - 1] + *qnext_len;
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

__global__ void get_vertex_betweenness_we(float *bc,
                                          const int *row_offsets,
                                          const int *cols,
                                          int nvertices,
                                          int *d,
                                          unsigned long long *sigma,
                                          float *delta,
                                          int *qcurr,
                                          int *qnext,
                                          int *stack,
                                          int *ends,
                                          int *next_source,
                                          bool normalized) {

    __shared__ int ind;
    __shared__ int s;

    int tid = (int) threadIdx.x;

    if (tid == 0) {
        ind = (int) blockIdx.x;
        s = ind;
    }

    __syncthreads();

    /*
     * For each vertex...
     */
    while (ind < nvertices) {

        for (int k = (int) threadIdx.x; k < nvertices; k += (int) blockDim.x) {
            if (k == s) {
                d[blockIdx.x * nvertices + k] = 0;
                sigma[blockIdx.x * nvertices + k] = 1;
            } else {
                d[blockIdx.x * nvertices + k] = INT_MAX;
                sigma[blockIdx.x * nvertices + k] = 0;
            }
            delta[blockIdx.x * nvertices + k] = 0;
        }

        __syncthreads();

        __shared__ int qcurr_len;
        __shared__ int qnext_len;
        __shared__ int stack_len;
        __shared__ int ends_len;
        __shared__ int depth;
        __shared__ bool done;

        if (tid == 0) {
            qcurr[blockIdx.x * nvertices] = s;
            qcurr_len = 1;
            qnext_len = 0;
            stack[blockIdx.x * nvertices] = s;
            stack_len = 1;
            ends[blockIdx.x * nvertices] = 0;
            ends[blockIdx.x * nvertices + 1] = 1;
            ends_len = 2;
            depth = 0;
            done = false;
        }

        __syncthreads();

        //Do first iteration separately since we already know the edges to traverse
        int start = (int) threadIdx.x + row_offsets[s];
        int end = row_offsets[s + 1];
        for (int r = start; r < end; r += (int) blockDim.x) {
            int w = cols[r];

            if (d[blockIdx.x * nvertices + w] == INT_MAX) {
                d[blockIdx.x * nvertices + w] = 1;
                int t = atomicAdd(&qnext_len, 1);
                qnext[blockIdx.x * nvertices + t] = w;
            }

            if (d[blockIdx.x * nvertices + w] ==
                (d[blockIdx.x * nvertices + s] + 1)) {
                atomicAdd(&sigma[blockIdx.x * nvertices + w], 1);
            }
        }
        __syncthreads();

        /*
         * If the next frontier is empty the traversal is complete.
         */
        if (qnext_len == 0) {
            done = true;
        } else {
            bfs_update_ds(&qcurr_len, qcurr,
                          &qnext_len, qnext,
                          &stack_len, stack,
                          &ends_len, ends,
                          &depth,
                          tid,
                          nvertices);
        }
        __syncthreads();

        /*
         * Bfs. One iteration for each depth of the graph.
         */
        while (!done) {

            __shared__ int next_index;
            if (tid == 0) {
                next_index = (int) blockDim.x;
            }
            __syncthreads();

            int k = (int) threadIdx.x;

            /*
             * For each vertex in the current frontier.
             *
             * Assigns a thread to each vertex such that edges from other
             * portions of the graph are not unnecessarily traversed.
             */
            while (k < qcurr_len) {
                int v = qcurr[blockIdx.x * nvertices + k];

                /*
                 * Add the neighbours of the vertex of the current frontier
                 * to the queue of the vertices of the next frontier.
                 */
                for (int r = row_offsets[v]; r < row_offsets[v + 1]; r++) {
                    int w = cols[r];

                    if (atomicCAS(&d[blockIdx.x * nvertices + w], INT_MAX,
                                  d[blockIdx.x * nvertices + v] + 1) ==
                        INT_MAX) {
                        int t = atomicAdd(&qnext_len, 1);
                        qnext[blockIdx.x * nvertices + t] = w;
                    }

                    if (d[blockIdx.x * nvertices + w] ==
                        (d[blockIdx.x * nvertices + v] + 1)) {
                        atomicAdd(&sigma[blockIdx.x * nvertices + w],
                                  sigma[blockIdx.x * nvertices + v]);
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
                bfs_update_ds(&qcurr_len, qcurr,
                              &qnext_len, qnext,
                              &stack_len, stack,
                              &ends_len, ends,
                              &depth,
                              tid,
                              nvertices);
            }
            __syncthreads();
        }

        //The elements at the end of the stack will have the largest distance from the source
        //Using the successor method, we can start from one depth earlier
        if (tid == 0) {
            depth = d[blockIdx.x * nvertices +
                      stack[blockIdx.x * nvertices + stack_len - 1]] -
                    1;
        }
        __syncthreads();

        //Dependency Accumulation (Madduri/Ediger successor method)
        while (depth > 0) {

            int start =
                    (int) threadIdx.x + ends[blockIdx.x * nvertices + depth];
            int end = ends[blockIdx.x * nvertices + depth + 1];
            for (int i = start; i < end; i += (int) blockDim.x) {
                int w = stack[blockIdx.x * nvertices + i];
                float dsw = 0;
                auto sw = (float) sigma[blockIdx.x * nvertices + w];
                for (int z = row_offsets[w]; z < row_offsets[w + 1]; z++) {
                    int v = cols[z];
                    if (d[blockIdx.x * nvertices + v] ==
                        (d[blockIdx.x * nvertices + w] + 1)) {
                        dsw += (sw /
                                (float) sigma[blockIdx.x * nvertices + v]) *
                               (1.0f + delta[blockIdx.x * nvertices + v]);
                    }
                }
                delta[blockIdx.x * nvertices + w] = dsw;
            }

            __syncthreads();
            if (tid == 0) {
                depth--;
            }
            __syncthreads();
        }

        if(normalized) {
            for (int i = (int) threadIdx.x; i < nvertices; i += (int) blockDim.x) {
                atomicAdd(&bc[i],
                          (2 * delta[blockIdx.x * nvertices + i]) /
                          (float) (nvertices * nvertices - 3 * nvertices + 2));
            }
        } else {
            for (int i = (int) threadIdx.x; i < nvertices; i += (int) blockDim.x) {
                atomicAdd(&bc[i], delta[blockIdx.x * nvertices + i]);
            }
        }


        if (tid == 0) {
            ind = atomicAdd(next_source, 1);
            s = ind;
        }
        __syncthreads();
    }
}

void compute_bc_gpu_we(matrix_pcsr_t *g, float *bc, bool normalized) {

    const unsigned int sm_count = get_sm_count();
    int next_source = (int) sm_count;

    unsigned long long *d_sigma;
    float *d_bc, *d_delta;
    int *d_row_offsets,
            *d_cols,
            *d_dist,
            *d_qcurr,
            *d_qnext,
            *d_stack,
            *d_endpoints,
            *d_next_source;

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
    cudaSafeCall(cudaMalloc((void **) &d_bc, g->nrows * sizeof(float)));
    cudaSafeCall(cudaMemset(d_bc, 0, g->nrows * sizeof(float)));

    /*
     * Load auxiliary arrays for bc.
     */
    cudaSafeCall(cudaMalloc((void **) &d_dist,
                            g->nrows * grid.x * sizeof(int)));
    cudaSafeCall(cudaMalloc((void **) &d_sigma,
                            g->nrows * grid.x * sizeof(unsigned long long)));
    cudaSafeCall(cudaMalloc((void **) &d_delta,
                            g->nrows * grid.x * sizeof(float)));
    cudaSafeCall(cudaMalloc((void **) &d_stack,
                            g->nrows * grid.x * sizeof(int)));
    cudaSafeCall(cudaMalloc((void **) &d_qcurr,
                            g->nrows * grid.x * sizeof(int)));
    cudaSafeCall(cudaMalloc((void **) &d_qnext,
                            g->nrows * grid.x * sizeof(int)));
    cudaSafeCall(cudaMalloc((void **) &d_endpoints,
                            (g->nrows + 1) * grid.x * sizeof(int)));

    /*
     * Load single-variables.
     */
    //    cudaCalloc((void **) &d_nedges, 1, sizeof(unsigned int));
    cudaSafeCall(cudaMalloc((void **) &d_next_source, sizeof(int)));
    cudaSafeCall(cudaMemcpy(d_next_source, &next_source,
                            sizeof(int),
                            cudaMemcpyHostToDevice));

    //    statistics->load_time[ntry] = chrono.elapsed();

    /*
     * Execute the bc computation.
     */
    get_vertex_betweenness_we<<<grid, block>>>(d_bc,
                                               d_row_offsets,
                                               d_cols,
                                               g->nrows,
                                               d_dist,
                                               d_sigma,
                                               d_delta,
                                               d_qcurr,
                                               d_qnext,
                                               d_stack,
                                               d_endpoints,
                                               d_next_source,
                                               normalized);

    cudaCheckError();

    //    statistics->bc_comp_time[ntry] = chrono.elapsed();

    //    unsigned int nedges;

    //    cudaSafeCall(cudaMemcpy(&nedges, d_nedges,
    //                            sizeof(unsigned int),
    //                            cudaMemcpyDeviceToHost));
    //    cudaSafeCall(cudaMemcpy(&next_source, d_next_source,
    //                            sizeof(int),
    //                            cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy(bc, d_bc,
                            g->nrows * sizeof(float),
                            cudaMemcpyDeviceToHost));

    /*
     * Count each edge only one time.
     */
    for (int k = 0; k < g->nrows; k++)
        bc[k] /= 2;

    /*
     * Device resource deallocation.
     */
    cudaSafeCall(cudaFree(d_row_offsets));
    cudaSafeCall(cudaFree(d_cols));
    //    cudaSafeCall(cudaFree(d_nedges));
    cudaSafeCall(cudaFree(d_next_source));
    cudaSafeCall(cudaFree(d_bc));
    cudaSafeCall(cudaFree(d_sigma));
    cudaSafeCall(cudaFree(d_dist));
    cudaSafeCall(cudaFree(d_delta));
    cudaSafeCall(cudaFree(d_stack));
    cudaSafeCall(cudaFree(d_qcurr));
    cudaSafeCall(cudaFree(d_qnext));
    cudaSafeCall(cudaFree(d_endpoints));

    //    statistics->unload_time[ntry] = chrono.elapsed();
    //    statistics.nedges_traversed = nedges;
}
