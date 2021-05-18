/****************************************************************************
 * @file bc.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions to compute the betweenness centrality of the vertices of an
 * undirected and unweighted graph stored as a sparse pattern matrix in CSR
 * format.
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
#ifndef SOCNETALGSONGPU_BC_H
#define SOCNETALGSONGPU_BC_H

#include "matstorage.h"
#include "common.h"
#include <climits>
#include <queue>
#include <stack>
#include <vector>

using std::queue;
using std::stack;

void compute_dep_acc(matrix_pcsr_t *g,
                     stack<int> S,
                     const unsigned long *sigma,
                     const int *d,
                     float *delta,
                     float *bc,
                     int s);

void compute_bfs(matrix_pcsr_t *g,
                 queue<int> Q,
                 stack<int> S,
                 unsigned long *sigma,
                 int *d);

void BC_dec_comp(matrix_pcsr_t *g, float *bc_scores, bool directed);

void get_vertex_betweenness(matrix_pcsr_t *g, float *bc_scores, bool directed);

void print_bc_scores(matrix_pcsr_t *g, const float *bc_scores, FILE *f);

void get_closeness() {
    // TODO
}

#endif //SOCNETALGSONGPU_BC_H