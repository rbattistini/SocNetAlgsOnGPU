/****************************************************************************
 * @file bc_statistics.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Compute and store statistics such as runtime and teps.
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
#ifndef BC_STATISTICS_H
#define BC_STATISTICS_H

#include "matds.h"
#include <cstring>

/**
 * @brief Structure of array that stores parameters used to measure runtime
 * and throughput.
 */
typedef struct stats_t {
    double load_time = 0;
    double unload_time = 0;
    double bc_comp_time = 0;
    double total_time = 0;
    double cpu_time = 0;
    unsigned long long nedges_traversed = 0;
} stats_t;

/**
 * @brief Computes the RMSE (Root-Mean-Square Error) of the bc scores obtained
 * by the CPU and the GPU.
 *
 * @note The RMSE describes the standard deviation of the differences
 * between the predicted bc scores of the CPU and the observed bc scores on
 * the GPU.
 *
 * @param nrows number of vertices
 * @param score_cpu
 * @param score_gpu
 * @return the maximum value of the RMSE
 */
double check_score(int nrows, const double *score_cpu, const double *score_gpu);

void print_stats(stats_t *s);

/**
 * @brief Dump runtime and throughput of the GPU algorithm to a file for each
 * run requested.
 *
 * @param stats structure with measured parameters
 * @param f file where the dump happens
 * @return 0 if successful, -1 if the stream was not closed correctly,
 * 1 if another error occurred
 */
int append_stats(stats_t *stats, char *fname, int technique_id);

/**
 * @brief Compute Traversed Edges Per Second for the BC algorithm on the GPU
 * as a measure of the throughput of the GPU.
 *
 * @note TEPS is defined as nedges / time, where nedges is computed as the
 * number of edges traversed during the forward propagation phase multiplied by
 * two.
 *
 * @note The Bfs performed during the forward propagation phase are
 * non-idempotent. This means that each vertex is visited only one time.
 *
 * @param nedges number of edges traversed during both the forward and backward
 * propagation phases of the Brandes algorithm on the GPU
 * @param time_elapsed includes both memory transfers and the time used by the
 * BC computation kernel
 * @return TEPS value
 */
double inline get_bc_teps(unsigned int nedges, double time_elapsed) {
    return (nedges * 2) / time_elapsed;
}

/**
 * @brief Dump scores to a file.
 *
 * @param nvertices
 * @param scores array that stores the scores to be dumped
 * @param fname file where the dump happens
 * @param score_name the name of the score to be dumped
 * @return 0 if successful, -1 if the stream was not closed correctly,
 * 1 if another error occurred
 */
int dump_scores(int nvertices,
                const int *degree_scores,
                const double *bc_scores,
                const double *cl_scores,
                char *fname);

#endif//BC_STATISTICS_H
