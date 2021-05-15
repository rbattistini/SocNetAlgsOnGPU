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

#include "gkernels.cuh"

#define NULL nullptr

/*
 * Each CUDA thread block processes its own set of roots and the threads
 * within the block traverse the graph cooperatively, starting at the
 * roots assigned to their block
 */

__global__ void del1deg(int *deg,
                        int nrows,
                        const int *row_offsets,
                        int *cols,
                        float *bc,
                        float *p,
                        bool *keep_on) {

    int v = (int) (blockIdx.x * blockDim.x + threadIdx.x);

    /*
     * One thread is assigned to each vertex. Only threads that are assigned
     * to a degree 1 vertex are active.
     */
    if (v < nrows && deg[v] == 1) {

        *keep_on = true;
        int w = row_offsets[v];

        /*
         * Mark edges to be removed.
         */
        cols[w] = -1;
        cols[row_offsets[w]] = -1;
        atomicAdd(deg + w, -1);

        /*
         * Update bc score and partial betweenness of the neighbour.
         */
        bc[w] += 2 * ((float) nrows - p[v] - p[w] - 2) * (p[v] + 1);
        p[w] += p[v] + 1;
    }
}

__global__ void vtx_par_bfs(int s,
                            int *d,
                            int *sigma,
                            int nrows,
                            int nnz,
                            const int *row_offsets,
                            const int *cols) {

    __shared__ int current_depth;
    __shared__ bool done;

    int tid = (int) threadIdx.x;

    /*
     * Initialize d and sigma.
     */
    for (int v = tid; v < nrows; v += (int) blockDim.x) {
        if (v == s) {
            d[v] = 0;
            sigma[v] = 1;
        } else {
            d[v] = -1;
            sigma[v] = 0;
        }
    }

    if (tid == 0) {
        done = false;
        current_depth = 0;
    }
    // wait for all threads to complete the initial configuration
    __syncthreads();

    /*
     * Compute the number of shortest paths (sigma) and the distance from s
     * (the root) to each vertex.
     */
    while (!done) {
        done = true;

        // wait for all threads to see if all vertices have been discovered
        __syncthreads();

        /*
         *  For each edge...
         */
        for (int v = tid; v < nnz; v += (int) blockDim.x) {

            /*
             * Only threads assigned to the vertices in the current forward
             * propagation frontier are active. Thus there is a work imbalance.
             */
            if (d[v] == current_depth) {

                /*
                 * If a thread is assigned to a vertex with a low degree and
                 * there is another thread with a high degree in the same warp
                 * then there is work imbalance.
                 */
                for (int i = row_offsets[v]; i < row_offsets[v + 1]; i++) {

                    int w = cols[i];

                    if (d[w] == -1) {
                        d[w] = d[v] + 1;
                        done = false;
                    }

                    if (d[w] == (d[v] + 1))
                        atomicAdd(&sigma[w], sigma[v]);
                }
            }
        }
        // wait for all threads to compute vertices of the current frontier
        __syncthreads();

        if (tid == 0)
            current_depth++;

        // wait for all threads to see the updated current depth
        __syncthreads();
    }
}

__global__ void vtx_par_dep_acc(int s,
                                const int *d,
                                const int *sigma,
                                float *delta,
                                float *bc,
                                int nrows,
                                int nnz,
                                const int *row_offsets,
                                const int *cols) {

    int tid = (int) threadIdx.x;

    /*
     * Initialize delta and bc scores.
     */
    for (int v = tid; v < nrows; v += (int) blockDim.x) {
        bc[v] = 0.0;
        delta[v] = 0.0;
    }

    __shared__ int current_depth;
    __shared__ bool done;

    if (tid == 0) {
        done = false;
        current_depth = 0;
    }
    __syncthreads();

    /*
     * Compute the dependency accumulation (delta) and the bc score of each
     * vertex.
     */
    while (!done) {
        __syncthreads();
        done = true;
        __syncthreads();

        /*
         *  For each edge...
         */
        for (int v = tid; v < nnz; v += (int) blockDim.x) {

            /*
             * Only threads assigned to the vertices in the current backward
             * propagation frontier are active.
             */
            if (d[v] == current_depth) {

                /*
                 * If a thread is assigned to a vertex with a low degree and
                 * there is another thread with a high degree in the same warp
                 * then there is work imbalance.
                 */
                for (int i = row_offsets[v]; i < row_offsets[v + 1]; i++) {

                    int w = cols[i];

                    if (d[v] == (d[w] - 1)) {
                        delta[v] +=
                                ((float) sigma[v] / (float) sigma[w]) *
                                (1.0f + delta[w]);
                    }

                    if (w != s) {
                        bc[w] += delta[w];
                    }
                }
            }
        }
        __syncthreads();

        if (tid == 0)
            current_depth++;
    }
}
