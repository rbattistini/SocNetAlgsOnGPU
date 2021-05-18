/****************************************************************************
 * @file spvb.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Modified version of serial algorithm of Brandes for computing
 * betweenness centrality optimized for social networks.
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

#include "spvb.h"

using std::queue;
using std::stack;

int spvb(const PUNGraph &g_i, float *bc_scores, float *p, bool directed) {

    int nit = 0, tree_nodes = 0;
    int nvertices = g_i->GetNodes();
    stack<int> S;

    for (int i = 0; i < nvertices; i++) {
        p[i] = 0.0f;
        bc_scores[i] = 0.0f;
    }

    /*
     * Initialize first set of 1-degree vertices.
     */
    for (TUNGraph::TNodeI NI = g_i->BegNI(); NI < g_i->EndNI(); NI++) {
        if (NI.GetDeg() == 1) {
            S.push(NI.GetId());
            tree_nodes++;
        }
    }

    do {

        int v = S.top();
        S.pop();

        const TUNGraph::TNodeI NI = g_i->GetNI(v);
        int w = NI.GetOutNId(0);

        bc_scores[w] += 2 * ((float) nvertices - p[v] - p[w] - 2) * (p[v] + 1);
        p[w] += p[v] + 1;

        nit++;
        g_i->DelEdge(v, w);
        g_i->DelNode(v);

        const TUNGraph::TNodeI NIw = g_i->GetNI(w);
        if (NIw.GetDeg() == 1) {
            S.push(w);
            tree_nodes++;
        }


    } while (!S.empty());

    //    printf("Nit: %d\n", nit);

    print_float_array(bc_scores, nvertices);
    print_float_array(p, nvertices);

    /*
     * Call the BC of Brandes procedure on the new graph g_i.
     */
    if (g_i->GetNodes() > 1)
        BC_mod_computation(g_i, p, bc_scores);

    if (!directed) {
        for (int k = 0; k < nvertices; k++)
            bc_scores[k] /= 2;
    }

    return tree_nodes;
}

void BC_mod_computation(const PUNGraph &g, const float *p, float *bc_scores) {

    int nvertices = g->GetNodes();
    for (TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {

        int s = NI.GetId();
        auto *sigma = (unsigned long *) malloc(nvertices * sizeof(unsigned long));
        auto *distance = (int *) malloc(nvertices * sizeof(int));
        auto *delta = (float *) malloc(nvertices * sizeof(float));

        if(sigma == 0 || distance == 0 || delta == 0) {
            ZF_LOGF("Memory allocation failed!");
            return;
        }

        for (int i = 0; i < nvertices; i++) {
            sigma[i] = 0;
            delta[i] = 0.0f;
            distance[i] = INT_MAX;
        }

        queue<int> Q;
        stack<int> S;

        sigma[s] = 1;
        distance[s] = 0;
        Q.push(s);

        while (!Q.empty()) {

            int v = Q.front();
            Q.pop();
            // update for the backward propagation phase
            S.push(v);

            const TUNGraph::TNodeI NI = g->GetNI(v);
            for (int e = 0; e < NI.GetOutDeg(); e++) {

                int w = NI.GetOutNId(e);

                /*
                 * If the vertex was not discovered, discover it and add to
                 * the queue of new vertices to visit. Update its distance from
                 * v.
                 */
                if (distance[w] == INT_MAX) {
                    Q.push(w);
                    distance[w] = distance[v] + 1;
                }

                /*
                 * If the vertex is "safe" give him all the power v has
                 * in terms of shortest paths crossing it.
                 */
                if (distance[w] == (distance[v] + 1)) {
                    sigma[w] += sigma[v];
                }
            }
        }

        while (!S.empty()) {

            int w = S.top();
            S.pop();

            const TUNGraph::TNodeI &NI = g->GetNI(w);
            for (int e = 0; e < NI.GetOutDeg(); e++) {

                int v = NI.GetOutNId(e);

                if (distance[v] == (distance[w] - 1)) {
                    delta[v] += (sigma[v] / (float) sigma[w]) *
                                (1.0f + delta[w] + p[w]);
                }
            }

            if (w != s) {
                bc_scores[w] += delta[w] * (1 + p[s]);
            }
        }

        free(sigma);
        free(distance);
        free(delta);
    }
}
