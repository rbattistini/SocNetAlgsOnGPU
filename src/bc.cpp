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

        while (!S.empty()) {

            int w = S.top();
            S.pop();

            for (int i = g->row_offsets[w]; i < g->row_offsets[w + 1]; i++) {
                int v = g->cols[i];
                if (d[v] == (d[w] - 1)) {
                    delta[v] +=
                            (sigma[v] / (double) sigma[w]) * (1.0f + delta[w]);
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
