/****************************************************************************
 * @file bc_vp_kernels.cu
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * Kernels for computing Betweenness centrality on a Nvidia GPU.
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
#ifndef BC_VP_KERNELS_CUH
#define BC_VP_KERNELS_CUH

#ifdef __CUDACC__

#include "device_props.cuh"
#include "matds.h"
#include <bc_statistics.h>
#include <common.h>

__global__ void get_vertex_betweenness_vpp(float *bc,
                                           const int *row_offsets,
                                           const int *cols,
                                           int nvertices,
                                           int *d,
                                           unsigned long long *sigma,
                                           float *delta,
                                           int *next_source,
                                           size_t pitch_d,
                                           size_t pitch_sigma,
                                           size_t pitch_delta);

void compute_bc_gpu_vpp(matrix_pcsr_t *g, float *bc);

#endif

#endif//BC_VP_KERNELS_CUH
