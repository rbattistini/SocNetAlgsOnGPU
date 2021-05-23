/****************************************************************************
 * @file bc.cpp
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

#include "bc.h"

void compute_ser_bc_cpu(matrix_pcsr_t *g, float *bc_scores, bool directed) {

    for (int i = 0; i < g->nrows; i++)
        bc_scores[i] = 0;
    //    int cnt = 0;
    for (int j = 0; j < g->nrows; j++) {

        int s = j;
        auto *sigma = (int *) malloc(
                g->nrows * sizeof(int));
        auto *d = (int *) malloc(g->nrows * sizeof(int));
        auto *delta = (double *) malloc(g->nrows * sizeof(double));
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

        //        printf("s: %d\n", s);
        //        print_int_array(sigma, g->nrows - 1);
        //        print_int_array(d, g->nrows - 1);

        while (!S.empty()) {

            int w = S.top();
            S.pop();

            for (int i = g->row_offsets[w]; i < g->row_offsets[w + 1]; i++) {
                int v = g->cols[i];
                if (d[v] == (d[w] - 1)) {
                    printf("v: %d, delta[v]: %.1f\n", v, delta[v]);
                    delta[v] +=
                            (sigma[v] / (double) sigma[w]) * (1.0f + delta[w]);
                }
            }

            if (w != s) {
                bc_scores[w] += delta[w];
            }
        }

        //        print_double_array(delta, g->nrows - 1);
        //        print_float_array(bc_scores, g->nrows - 1);

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

void compute_par_bc_cpu(matrix_pcsr_t *g_tmp, float *bc_cpu) {
    typedef boost::adjacency_list<> graph;

    Stopwatch chrono;
    int nnz = g_tmp->row_offsets[g_tmp->nrows];

    graph g((unsigned long) g_tmp->nrows);

    auto rows = (int *) malloc(nnz * sizeof(int));
    expand_row_pointer(g_tmp->nrows, g_tmp->row_offsets, rows);

    for (int i = 0; i < nnz; i++) {
        boost::add_edge((unsigned long) rows[i], (unsigned long) g_tmp->cols[i],
                        g);
    }

    /*
     * Compute BC with the algorithm that uses multithreading of the BGL.
     */
    boost::shared_array_property_map<double, boost::property_map<graph,
            boost::vertex_index_t>::const_type>
            centrality_map(num_vertices(g), get(boost::vertex_index, g));

    chrono.start();
    boost::brandes_betweenness_centrality(g, centrality_map);
    chrono.stop();

    for (int i = 0; i < g_tmp->nrows; i++) {
        bc_cpu[i] = (float) (centrality_map[i]);
    }

    /*
     * Count each edge only one time.
     */
    for (int k = 0; k < g_tmp->nrows; k++)
        bc_cpu[k] /= 2;
}
