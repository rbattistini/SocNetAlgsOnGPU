/****************************************************************************
 *
 * graphs.cpp - Algorithms for graph manipulation
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
 * --------------------------------------------------------------------------
 *
 * TODO compute_network_bc_score
 *
 ****************************************************************************/

#include "graphs.h"
#include "utils.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <stack>
#include <vector>

using std::queue;
using std::stack;

void print_gprops(gprops_t *gp) {
    printf("Directed: %s\n", (gp->is_directed == 1) ? "yes" : "no");
    printf("Weighted: %s\n", (gp->is_weighted == 1) ? "yes" : "no");
    printf("Connected: %s\n", (gp->is_connected == 1) ? "yes" : "no");
}

void extract_und_subgraph(const int *vertices, int nvertices, matrix_pcsr_t *g,
                          matrix_pcsr_t *m) {

    int c = 0, *rows, *cols;
    std::vector<int> vcols, vrows;
    matrix_pcoo_t subgraph;

    /*
     * Store neighbours of each vertex in the cc in vrows and vcols.
     */
    for (int i = 0; i < nvertices; i++) {
        int vi = vertices[i];
        for (int j = g->row_offsets[vi]; j < g->row_offsets[vi + 1]; j++) {
            vrows.push_back(vi);
            vcols.push_back(g->cols[j]);
            c++;
        }
    }

    //    print_array()
    //    printf("c: %d\n", c);
    rows = (int *) malloc(c * sizeof(*rows));
    cols = (int *) malloc(c * sizeof(*cols));

    /*
     * Compute row and column indexes for the new submatrix.
     */
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < nvertices; j++) {
            if (vertices[j] == vcols[i]) {
                cols[i] = j;
            }
            if (vertices[j] == vrows[i]) {
                rows[i] = j;
            }
        }
    }

    subgraph.nnz = c;
    subgraph.nrows = nvertices;
    subgraph.rows = rows;
    subgraph.cols = cols;

    //    print_matrix_coo(&subgraph);

    pcoo_to_pcsr(&subgraph, m);
    free_matrix_pcoo(&subgraph);
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

int get_cc(matrix_pcsr_t *g, components_t *ccs) {

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

    return cc_count;
}

void compute_bfs(matrix_pcsr_t *g, queue<int> Q, stack<int> S,
                 unsigned long *sigma, int *d) {

    while (!Q.empty()) {

        int v = Q.front();
        Q.pop();
        S.push(v);

        for (int k = g->row_offsets[v]; k < g->row_offsets[v + 1]; k++) {

            int w = g->cols[k];

            /*
             * If the vertex was not discovered, discover it and add to
             * the queue of new vertices to visit. Update its d from
             * v.
             */
            if (d[w] == INT_MAX) {
                Q.push(w);
                d[w] = d[v] + 1;
            }

            /*
             * If the vertex is "safe" give him all the power v has
             * in terms of shortest paths crossing it.
             */
            if (d[w] == (d[v] + 1)) {
                sigma[w] += sigma[v];
            }
        }
    }
}

void compute_dep_acc(matrix_pcsr_t *g, stack<int> S,
                     const unsigned long *sigma, const int *d, float *delta,
                     float *bc, int s) {
    while (!S.empty()) {

        int w = S.top();
        S.pop();

        for (int i = g->row_offsets[w]; i < g->row_offsets[w + 1]; i++) {
            int v = g->cols[i];
            if (d[v] == (d[w] - 1)) {
                delta[v] +=
                        ((float) sigma[v] / (float) sigma[w]) *
                        (1.0f + delta[w]);
            }
        }

        if (w != s) {
            bc[w] += delta[w];
        }
    }
}

void BC_dec_comp(matrix_pcsr_t *g, float *bc_scores, bool directed) {

    for (int i = 0; i < g->nrows; i++)
        bc_scores[i] = 0.0;

    for (int j = 0; j < g->nrows; j++) {

        /*
         * Workspace setup.
         */
        int s = j;
        auto *sigma = (unsigned long *) malloc(
                g->nrows * sizeof(unsigned long));
        auto *d = (int *) malloc(g->nrows * sizeof(int));
        auto *delta = (float *) malloc(g->nrows * sizeof(float));
        assert(sigma);
        assert(d);
        assert(delta);

        for (int i = 0; i < g->nrows; i++) {
            sigma[i] = 0;
            delta[i] = 0.0f;
            d[i] = INT_MAX;
        }

        queue<int> Q;
        stack<int> S;

        sigma[s] = 1;
        d[s] = 0;
        Q.push(s);

        /*
         *  First forward propagation phase.
         */
        compute_bfs(g, Q, S, sigma, d);

        /*
         * Second backward propagation phase.
         */
        compute_dep_acc(g, S, sigma, d, delta, bc_scores, j);

        /*
         * Cleanup.
         */
        free(sigma);
        free(d);
        free(delta);
    }

    /*
     * Scores are duplicated if the graph is undirected because each edge is
     * counted two times.
     */
    if (!directed) {
        for (int k = 0; k < g->nrows; k++)
            bc_scores[k] /= 2;
    }
}

void BC_computation(matrix_pcsr_t *g, float *bc_scores, bool directed) {

    for (int i = 0; i < g->nrows; i++)
        bc_scores[i] = 0;

    for (int j = 0; j < g->nrows; j++) {

        int s = j;
        auto *sigma = (unsigned long *) malloc(
                g->nrows * sizeof(unsigned long));
        auto *d = (int *) malloc(g->nrows * sizeof(int));
        auto *delta = (float *) malloc(g->nrows * sizeof(float));
        assert(sigma);
        assert(d);
        assert(delta);

        for (int i = 0; i < g->nrows; i++) {
            sigma[i] = 0;
            delta[i] = 0.0f;
            d[i] = INT_MAX;
        }

        queue<int> Q;
        stack<int> S;

        sigma[s] = 1;
        d[s] = 0;
        Q.push(s);

        while (!Q.empty()) {

            int v = Q.front();
            Q.pop();
            // update for the backward propagation phase
            S.push(v);

            //            printf("%d | ", v);
            for (int k = g->row_offsets[v]; k < g->row_offsets[v + 1]; k++) {

                int w = g->cols[k];
                //                printf("%d ", w);

                /*
                 * If the vertex was not discovered, discover it and add to
                 * the queue of new vertices to visit. Update its d from
                 * v.
                 */
                if (d[w] == INT_MAX) {
                    Q.push(w);
                    d[w] = d[v] + 1;
                }

                /*
                 * If the vertex is "safe" give him all the power v has
                 * in terms of shortest paths crossing it.
                 */
                if (d[w] == (d[v] + 1)) {
                    sigma[w] += sigma[v];
                }
            }

            //            printf("\n");
        }

        while (!S.empty()) {

            int w = S.top();
            S.pop();

            for (int i = g->row_offsets[w]; i < g->row_offsets[w + 1]; i++) {
                int v = g->cols[i];
                if (d[v] == (d[w] - 1)) {
                    delta[v] +=
                            (sigma[v] / (float) sigma[w]) * (1.0f + delta[w]);
                }
            }

            if (w != s) {
                bc_scores[w] += delta[w];
            }
        }

        free(sigma);
        free(d);
        free(delta);
    }

    /*
     * Scores are duplicated if the graph is undirected because each edge is
     * counted two times.
     */
    if (!directed) {
        for (int k = 0; k < g->nrows; k++)
            bc_scores[k] /= 2;
    }
}

void print_bc_scores(matrix_pcsr_t *g, const float *bc_scores, FILE *fout) {

    unsigned int nvertices = g->nrows;
    unsigned int nedges = g->row_offsets[nvertices];

    if (fout != nullptr) {
        fprintf(fout, "Number of vertices: %d\n", nvertices);
        fprintf(fout, "Number of edges: %d\n", nedges);

        for (size_t i = 0; i < nvertices; i++) {
            fprintf(fout, "%.2f\n", bc_scores[i]);
        }
    } else {
        fprintf(stderr, "Failed to create output file\n");
    }
}

void compute_degrees_undirected(matrix_pcoo_t *g, int *degree) {

    if (!check_matrix_pcoo_init(g)) {
        fprintf(stderr, "The graph is not initialized");
        return;
    }

    int *rows = g->rows; // row indices of A
    int nnz = g->nnz;    // number of nnz in A
    int nrows = g->nrows;// number of rows in A

    fill(degree, nrows, 0);

    /*
     * Compute number of non-zero entries per row of A.
     */
    for (int n = 0; n < nnz; n++) {
        degree[rows[n]]++;
    }
}

void compute_degrees_directed(matrix_pcoo_t *g, int *in_degree,
                              int *out_degree) {

    if (g->rows == nullptr) {
        fprintf(stderr, "The graph is not initialized");
        return;
    }

    int *rows = g->rows;  // row indices of A
    int nnz = g->nnz;     // number of nnz in A
    int length = g->nrows;// number of rows and columns in A
    int *cols = g->cols;  // column indices of A

    //    offsets = (int *) malloc((length + 1) * sizeof(*offsets));
    fill(in_degree, length, 0);
    fill(out_degree, length, 0);

    /*
     * Compute number of non-zero entries per column of A.
     */
    for (int n = 0; n < nnz; n++) {
        out_degree[rows[n]]++;
    }

    /*
     * Compute number of non-zero entries per row of A.
     */
    for (int n = 0; n < nnz; n++) {
        in_degree[cols[n]]++;
    }
}

void free_ccs(components_t *ccs) {
    free(ccs->cc_size);
    free(ccs->array);
}
