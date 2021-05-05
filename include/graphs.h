/****************************************************************************
 *
 * graphs.h - Algorithms for graph manipulation
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

#ifndef GRAPHS_H
#define GRAPHS_H

#include <vector>
#include <queue>
#include <climits>
#include <cstdio>
#include "matstorage.h"

typedef int vertex;

typedef struct gprops_t {
    bool is_directed;
    bool is_weighted;
    bool is_connected;
} gprops_t;

typedef struct components_t {
    vertex *ccs_array;  // vertices ids of each cc
    int *ccs_size;      // size of each cc
} components_t;

void BC_computation(matrix_pcsr_t *g, float *bc_scores, bool directed);

void print_bc_scores(matrix_pcsr_t *g, const float *bc_scores, FILE* fout);

void extract_und_subgraph(const vertex *vertices, int nvertices, matrix_pcsr_t *g,
                          matrix_pcsr_t *m);

int* DFS_visit(matrix_pcsr_t *g, bool *visited, int s, int *cc_size);

int get_cc(matrix_pcsr_t *g, components_t *ccs);

void compute_degrees_undirected(matrix_pcoo_t* g, int *degree);

void compute_degrees_directed(matrix_pcoo_t*  g, int *in_degree,
                              int *out_degree);

#endif //GRAPHS_H
