/****************************************************************************
 * @file bc_we_kernel.cuh
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

#pragma once
#ifndef BC_KERNELS_PITCHED_CUH
#define BC_KERNELS_PITCHED_CUH

#ifdef __CUDACC__

#include "device_props.cuh"
#include "matds.h"
#include <bc_statistics.h>
#include <common.h>

/**
 * @brief Computes Betweenness Centrality exploiting two levels of
 * parallelism:
 *
 * - coarse grained because a parallel BFS or dependency accumulation is
 * assigned to each Streaming Multiprocessor;
 * - fine grained because each in each BFS or dependency accumulation work
 * is distributed among threads according to a certain strategy.
 *
 * Computes dependency accumulation and bc scores using Brandes' recursive
 * formula. For the parallel dependency accumulation it is used the method
 * proposed by Madduri and Ediger.
 *
 * The algorithm uses the work efficient technique and the optimizations
 * proposed by McLaughlin and Bader.
 *
 * @note Uses pitched memory.
 *
 * @cite brandes_faster_2001
 * @cite madduri_faster_2009
 * @cite mclaughlin_scalable_2014
 *
 * @param[out] bc
 * @param[in] row_offsets
 * @param[in] cols
 * @param[in] nvertices
 * @param[in] d
 * @param[in] sigma
 * @param[out] delta
 * @param[in] curr_queue
 * @param[in] next_queue
 * @param[in] stack
 * @param[in] endpoints
 * @param[in] next_source
 * @param[in] pitch_d
 * @param[in] pitch_sigma
 * @param[in] pitch_delta
 * @param[in] pitch_qcurr
 * @param[in] pitch_qnext
 * @param[in] pitch_stack
 * @param[in] pitch_endpoints
 */
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
                                           size_t pitch_endpoints);

void compute_bc_gpu_wep(matrix_pcsr_t *g, double *bc, stats_t *stats);

#endif

#endif//BC_KERNELS_PITCHED_CUH
