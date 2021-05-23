/****************************************************************************
 * @file graph_gen.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Program to generate synthetic undirected and unweighted graphs
 * dataset using the SNAP Library. Based on examples available in the SNAP
 * documentation.
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

#include "matio.h"
#include "mmio.h"
#include <Snap.h>
#include <cstdlib>
#include <iostream>

enum GraphType {
    Random      = 1,    // Erdős-Rényi random graph
    SmallWorld  = 2,    // Watts-Strogatz Small World graph
    ScaleFree   = 3     // Barabási-Albert Scale Free graph
};

using namespace TSnap;

/*
 * Basic configuration parameters.
 */
const int min_vertices = 20;
const int max_vertices = 5000;

/**
 * @brief Parsing with error handling.
 *
 * @param p pointer to the null-terminated byte string to be interpreted
 * @param endp pointer to a pointer to character
 * @param base base of the interpreted integer value
 * @return the parsed number or -1 if unsuccessful
 */
static long strtol_wcheck(const char *p, char* endp, int base) {
    int done = 0;
    long i = -1;

    while(!done) {
        /*
         * errno can be set to any non-zero value by a library function call
         * regardless of whether there was an error, so it needs to be cleared
         * in order to check the error set by strtol.
         */
        errno = 0;
        i = strtol(p, &endp, base);

        if (p == endp)
            done = 1;

        const bool range_error = errno == ERANGE;
        ZF_LOGI("Extracted '%.*s', strtol returned %ld.",
                (int)(endp - p), p, i);
        p = endp;

        if (range_error) {
            ZF_LOGE("Range error occurred");
            return -1;
        }

        putchar('\n');
    }

    return i;
}

static int bin_coeff(int n, int k) {

    int C[k + 1];
    memset(C, 0, sizeof(C));

    C[0] = 1;

    for (int i = 1; i <= n; i++) {
        // Compute next row of pascal triangle using the previous row
        for (int j = min(i, k); j > 0; j--)
            C[j] = C[j] + C[j - 1];
    }
    return C[k];
}

static int write_snap_to_mtx(FILE *f, const PUNGraph &g) {

    /*
     * Write the banner.
     */
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_pattern(&matcode);
    mm_set_symmetric(&matcode);

    if (mm_write_banner(f, matcode)) {
        ZF_LOGF("Error writing mm banner");
        return EXIT_FAILURE;
    }

    /*
     * Write the header and the values.
     */
    if (mm_write_mtx_crd_size(f, g->GetNodes(), g->GetNodes(),
                              g->GetEdges())) {
        ZF_LOGF("Error writing header");
        return EXIT_FAILURE;
    }

    for (TUNGraph::TEdgeI EI = g->BegEI(); EI < g->EndEI(); EI++) {
        if (fprintf(f, "%d %d\n", EI.GetSrcNId() + 1,
                    EI.GetDstNId() + 1) < 0) {
            ZF_LOGF("Error writing values");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

/**
 * @brief Print graph's basic statistics, namely the name of the graph, the
 * number of edges and of vertices and whether the graph is empty and
 * connected.
 *
 * @tparam PGraph type of graph, directed or undirected
 * @param s pointer to the name of the graph
 * @param G pointer to the object representing the graph
 */
template<class PGraph>
void PrintGStats(const char *s, PGraph G) {
    printf("statistics of undirected, unweighted graph in %s: \n"
           "nodes:\t\t %d\n"
           "edges:\t\t %d\n"
           "isEmpty:\t %s\n"
           "isConnected:\t %s\n",
           s, G->GetNodes(), G->GetEdges(),
           G->Empty() ? "yes" : "no",
           IsConnected(G) ? "yes" : "no");
}

int main(int argc, char *argv[]) {

    int v, gtype;
    char *fname;
    TInt::Rnd.PutSeed(0);   // random number generator
    PUNGraph g;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s [filename] [graph_type] [number_vertices]\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    fname = argv[1];
    gtype = (int) strtol_wcheck(argv[2], nullptr, 10);
    v = (int) strtol_wcheck(argv[3], nullptr, 10);

    if (v < min_vertices || v > max_vertices) {
        ZF_LOGF("A graph must have at least 20 nodes and no more "
                     "than 1000 nodes");
        return EXIT_FAILURE;
    }

    switch (gtype) {
        case Random: {
            /*
             * Number of edges, computed so that there is a high probability of
             * generating a connected random graph.
             */
            const int e = (int) (bin_coeff(v, 2) * (log(v) / v)) + 2;
            g = GenRndGnm<PUNGraph>(v, e, false, TInt::Rnd);
            break;
        }
        case SmallWorld: {
            // number of neighbours to which each node is connected (Out-Degree)
            const int k = (int) (0.8 * v);
            // rewiring probability
            const float p = 0.03;
            g = GenSmallWorld(v, k, p, TInt::Rnd);
            break;
        }
        case ScaleFree: {
            // number of neighbours to which each node is connected (Out-Degree)
            const int k = (int) (0.8 * v);
            g = GenPrefAttach(v, k);
            break;
        }
        default:
            ZF_LOGF("Graph types accepted are: \n"
                         "[1]: Random graph \n"
                         "[2]: Small World graph \n"
                         "[3]: Scale Free graph");
            return EXIT_FAILURE;
    }

    /*
     * Print graph basic statistics.
     */
    PrintGStats(fname, g);

    /*
     * Write SNAP graph to a Matrix Market file.
     */
    FILE *f = fopen(fname, "w");

    if (write_snap_to_mtx(f, g))
        return EXIT_FAILURE;

    close_stream(f);

    return EXIT_SUCCESS;
}
