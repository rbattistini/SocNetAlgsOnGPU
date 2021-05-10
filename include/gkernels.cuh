/****************************************************************************
 *
 * gkernels.cuh - Kernels for computing Betweenness centrality on a Nvidia GPU
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

#pragma one
#ifndef SOCNETALGSONGPU_GKERNELS_CUH
#define SOCNETALGSONGPU_GKERNELS_CUH

/**
 *
 * @param deg
 * @param nrows
 * @param row_offsets
 * @param cols
 * @param bc
 * @param p
 * @param keep_on
 */
__global__ void del1deg(int *deg,
                        int nrows,
                        const int *row_offsets,
                        int *cols,
                        float *bc,
                        float *p,
                        bool *keep_on);

/**
 *
 * @param s
 * @param d
 * @param sigma
 * @param nrows
 * @param nnz
 * @param row_offsets
 * @param cols
 */
__global__ void vtx_par_bfs(int s,
                            int *d,
                            int *sigma,
                            int nrows,
                            int nnz,
                            const int *row_offsets,
                            const int *cols);

/**
 *
 * @param s
 * @param d
 * @param sigma
 * @param delta
 * @param bc
 * @param nrows
 * @param nnz
 * @param row_offsets
 * @param cols
 */
__global__ void vtx_par_dep_acc(int s,
                                const int *d,
                                const int *sigma,
                                float *delta,
                                float *bc,
                                int nrows,
                                int nnz,
                                const int *row_offsets,
                                const int *cols);

#endif //SOCNETALGSONGPU_GKERNELS_CUH
