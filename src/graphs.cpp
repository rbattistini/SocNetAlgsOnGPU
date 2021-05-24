/****************************************************************************
 * @file graphs.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Algorithms for graphs manipulation.
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

#include "graphs.h"
#include <algorithm>
#include <ecc.h>

using std::queue;
using std::stack;

void print_graph_overview(matrix_pcsr_t *g, int *degree) {

    int nvertices = g->nrows;
    int nedges = g->row_offsets[g->nrows];

    printf("Graph overview:\n\n");
    printf("\tVertices:\t\t%d\n", nvertices);
    printf("\tEdges:\t\t\t%d\n", nedges);
    printf("\tDensity:\t\t%.2f %%\n", get_density(nvertices, nedges) * 100);
    printf("\tMax degree: \t\t%d\n", degree[get_max_idx(degree, nvertices)]);
//    printf("\tDiameter: \t\t%d\n", get_diameter(g));

    print_separator();
}

void print_graph_properties(gprops_t *gp) {

    printf("Graph properties:\n\n");
    printf("\tDirected: \t\t%s\n", (gp->is_directed == 1) ? "yes" : "no");
    printf("\tWeighted: \t\t%s\n", (gp->is_weighted == 1) ? "yes" : "no");
    printf("\tConnected: \t\t%s\n", (gp->is_connected == 1) ? "yes" : "no");
    printf("\tHas self loops: \t%s\n", (gp->has_self_loops == 1) ? "yes" : "no");

    print_separator();
}

int *DFS_visit(matrix_pcsr_t *g, bool *visited, int s, int *cc_size) {

    /*
     * Let the STL vector manage memory allocation.
     */
    std::stack<int> S;
    std::vector<int> subgraph_vertices;
    S.push(s);

    while (!S.empty()) {
        s = S.top();
        S.pop();

        if (!visited[s]) {
            visited[s] = true;
            subgraph_vertices.push_back(s);
        }

        /*
         * Get all adjacent vertices of the vertex s.
         * If a adjacent has not been visited, then push it to the S.
         */
        for (int i = g->row_offsets[s]; i < g->row_offsets[s + 1]; i++) {
            int v = g->cols[i];
            if (!visited[v]) {
                S.push(v);
            }
        }
    }

    /*
     * Copy the vector into an array of known size.
     */
    *cc_size = subgraph_vertices.size();
    int *cc_array = stlvector_to_array_int(subgraph_vertices, *cc_size);

    return cc_array;
}

void get_cc(matrix_pcsr_t *g, components_t *ccs) {

    auto visited = (bool *) malloc(g->nrows * sizeof(bool));
    assert(visited);
    int cc_count = 0;
    std::vector<int> ccs_array, ccs_size;

    for (int i = 0; i < g->nrows; i++)
        visited[i] = false;

    for (int i = 0; i < g->nrows; i++) {
        if (!visited[i]) {
            int cc_size = 0;
            int *cc_array = DFS_visit(g, visited, i, &cc_size);

            for (int j = 0; j < cc_size; j++) {
                ccs_array.push_back(cc_array[j]);
            }

            ccs_size.push_back(cc_size);
            cc_count++;
        }
    }

    int *tmp_ccs_array = stlvector_to_array_int(ccs_array, ccs_array.size());
    ccs->array = tmp_ccs_array;
    int *tmp_ccs_size = stlvector_to_array_int(ccs_size, ccs_size.size());
    ccs->cc_size = tmp_ccs_size;
    ccs->cc_count = cc_count;
}

void BFS_visit(matrix_pcsr_t *g, int *d, int s) {

    queue<int> Q;
    d[s] = 0;
    Q.push(s);

    while (!Q.empty()) {

        int v = Q.front();
        Q.pop();

        for (int k = g->row_offsets[v]; k < g->row_offsets[v + 1]; k++) {

            int w = g->cols[k];

            /*
             * If the vertex was not discovered, discover it and add to
             * the queue of new vertices to visit. Update its distance from
             * v.
             */
            if (d[w] == INT_MAX) {
                Q.push(w);
                d[w] = d[v] + 1;
            }
        }
    }
}

void extract_subgraph(const int *vertices,
                      int nvertices,
                      matrix_pcsr_t *A,
                      matrix_pcsr_t *C) {

    matrix_pcsr_t R, Q;

    get_R_matrix(&R, vertices, nvertices, A->nrows);
    transpose(&R, &Q);
    spref(&R, A, &Q, C);

    free_matrix_pcsr(&R);
    free_matrix_pcsr(&Q);
}

void get_largest_cc(matrix_pcsr_t *A, matrix_pcsr_t *C, components_t *ccs) {

    int max_idx = get_max_idx(ccs->cc_size, ccs->cc_count);
    int largest_cc_size = ccs->cc_size[max_idx];

    int *largest_cc_vertices =
            (int *) malloc(largest_cc_size * sizeof(*largest_cc_vertices));
    assert(largest_cc_vertices);

    int start, end;

    if (max_idx == 0) {
        start = 0;
        end = largest_cc_size;
    } else {
        int csum = 0;
        for (int i = 0; i < max_idx; i++)
            csum += ccs->cc_size[i];
        start = csum;
        end = start + ccs->cc_size[max_idx];
    }

    for (int i = start, j = 0; i < end; j++, i++)
        largest_cc_vertices[j] = ccs->array[i];

    std::sort(largest_cc_vertices, largest_cc_vertices + largest_cc_size);
    extract_subgraph(largest_cc_vertices, largest_cc_size, A, C);
    free(largest_cc_vertices);
}

void free_ccs(components_t *ccs) {
    free(ccs->cc_size);
    free(ccs->array);
}
