/****************************************************************************
 * @file graphs.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * Algorithms for graphs manipulation.
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
 * --------------------------------------------------------------------------
 *
 * TODO: Compute bc score of the network
 *
 ****************************************************************************/

#pragma once
#ifndef GRAPHS_H
#define GRAPHS_H

#include "matstorage.h"
#include "spmatops.h"
#include "utils.h"
#include <climits>
#include <queue>
#include <stack>

typedef struct gprops_t {
    bool is_directed;
    bool is_weighted;
    bool is_connected;
    bool has_self_loops;
} gprops_t;

typedef struct components_t {
    int *array;     // vertices ids of each cc
    int *cc_size;   // size of the i-th cc at index i
    int cc_count;
} components_t;

void print_gprops(gprops_t *gp);

void BFS_visit(matrix_pcsr_t *g, int *d, int s);

int* DFS_visit(matrix_pcsr_t *g, bool *visited, int s, int *cc_size);

/**
 * Get the largest cc and extract a subgraph from it.
 *
 * @param A
 * @param C
 * @param ccs
 */
void get_largest_cc(matrix_pcsr_t *A,
                    matrix_pcsr_t *C,
                    components_t *ccs);

void get_cc(matrix_pcsr_t *g, components_t *ccs);

int get_diameter(matrix_pcsr_t *g);

void free_ccs(components_t *ccs);

void extract_subgraph(const int *vertices,
                      const int nvertices,
                      matrix_pcsr_t *A,
                      matrix_pcsr_t *C);

#endif//GRAPHS_H
