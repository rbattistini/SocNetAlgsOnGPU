/****************************************************************************
 * @file cl_kernels.cuh
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Kernel for computing Closeness centrality.
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

#pragma once
#ifndef CL_KERNELS_CUH
#define CL_KERNELS_CUH

#ifdef __CUDACC__

#include "device_props.cuh"
#include "matds.h"
#include <bc_statistics.h>
#include <common.h>

/**
 * @brief Computes Closeness Centrality using pitched memory and exploiting
 * two levels of parallelism:
 *
 * - coarse grained because a parallel BFS is assigned to each Streaming
 * Multiprocessor;
 * - fine grained because each in each BFS work is distributed among threads.
 *
 * @param[out] cl
 * @param[in] rows
 * @param[in] cols
 * @param[in] nvertices
 * @param[in] nnz
 * @param[in, out] d
 * @param[in] next_source
 * @param[in] pitch_d
 */
__global__ void get_closeness_p(double *cl,
                                const int *rows,
                                const int *cols,
                                int nvertices,
                                int nnz,
                                int *d,
                                int *next_source,
                                size_t pitch_d);

void compute_cl_gpu_p(matrix_pcsr_t *g, double *cl, stats_t *stats);

#endif

#endif//CL_KERNELS_CUH
