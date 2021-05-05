/****************************************************************************
 *
 * betweenness.cpp - Serial algorithm for computing betweenness centrality
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

#include <cassert>
#include "Snap.h"
#include "utils.h"
#include "matio.h"
#include "matstorage.h"
#include "spvb.h"
#include "graphs.h"

int main( int argc, char *argv[] ) {

    if (argc != 3) {
        fprintf(stderr, "Usage: %s [input filename] [output_filename]\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    float *bc_scores, *bc_scores_mod, *partial_bc_scores;
    int *degree;
    const char *fname = argv[1];

    matrix_pcsr_t m_csr;
    matrix_pcoo_t m_coo;
    gprops_t gp;
    query_gprops(fname, &gp);

    /*
     * Load matrix in COO format.
     */
    if(read_matrix(argv[1], &m_coo, &gp)) {
        fprintf(stderr, "Error reading matrix\n");
        return EXIT_FAILURE;
    }

    /*
     * Load the matrix in COO as an undirected graph in SNAP.
     * This way enhanced algorithm can run efficiently.
     */
//    matrix_pcoo_t subgraph_coo;
//    matrix_pcsr_t subgraph_csr;
    degree = (int*) malloc(m_coo.nrows * sizeof(int));
    assert(degree);
    compute_degrees_undirected(&m_coo, degree);

    /*
     * Get connected components of undirected graph.
     */
    pcoo_to_pcsr(&m_coo, &m_csr);
//    print_matrix_csr(&m_csr);
//
//    components_t cc_array;
//    int cc_count = get_cc(&m_csr, &cc_array);
//    gp.is_connected = (cc_count == 1);
//    printf("cc: %d\n", cc_count);


//    PUNGraph g = TUNGraph::New();
//    for(int i = 0; i < m_coo.nrows; i++) {
//        g->AddNode(i);
//    }
//
//    for(int i = 0; i < m_coo.nnz; i++) {
//        g->AddEdge(m_coo.rows[i], m_coo.cols[i]);
//    }
//
    /*
     * Extract ccs and induced subgraphs.
     */
//    TCnComV cc_array;
//    TSnap::GetWccs(g, cc_array);
//
//    PUNGraph subg = TSnap::GetSubGraph(g, TIntV::GetV(cc_array[0]));

    /*
     * Convert the subgraphs obtained in COO, then in CSR.
     */

    // traverse the nodes
//    int i = 0;
//    degree = (int*) malloc(m_coo.nrows * sizeof(*degree));
//    assert(degree);
//
//    for(TUNGraph::TNodeI NI = g->BegNI(); NI < g->EndNI(); NI++) {
//        degree[i] = NI.GetDeg();
//        i++;
//    }
//
//    // traverse the edges
//    i = 0;
//    subgraph_coo.nnz = m_coo.nrows;
//    subgraph_coo.nrows = m_coo.nrows;
//
//    for(TUNGraph::TEdgeI EI = g->BegEI(); EI < g->EndEI(); EI++) {
//        subgraph_coo.rows[i] = EI.GetSrcNId();
//        subgraph_coo.cols[i] = EI.GetDstNId();
//        printf("edge (%d, %d)\n", EI.GetSrcNId(), EI.GetDstNId());
//        i++;
//    }

//    pcoo_to_pcsr(&subgraph_coo, &subgraph_csr);
//    print_matrix_coo(&m_coo);
//    print_matrix_csr(&m_csr);
//    print_edge_list(m_csr.row_offsets, m_csr.cols, m_csr.nrows);

    /*
     * Compute BC using CSR.
     */
    bc_scores = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
    assert(bc_scores);
    BC_computation(&m_csr, bc_scores, gp.is_directed);
    FILE *fout = fopen(argv[2], "w");
    print_bc_scores(&m_csr, bc_scores, stdout);
    fclose(fout);

    /*
     * Compute BC with enhanced algorithm, using TUNGraph.
     */
//    bc_scores_mod = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
//    partial_bc_scores = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
//    assert(bc_scores_mod);
//    assert(partial_bc_scores);
//    spvb(g, degree, bc_scores_mod, partial_bc_scores, gp.is_directed);
//    FILE *fout = fopen(argv[2], "w");
//    print_bc_scores(m_csr, bc_scores, stdout);
//    fclose(fout);

    /*
     * Cleanup.
     */
    free_matrix_pcoo(&m_coo);
    free_matrix_pcsr(&m_csr);
//    free(bc_scores);
    free(degree);

    /*
     * Closing standard streams with error checking.
     */
    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);

    return EXIT_SUCCESS;
}
