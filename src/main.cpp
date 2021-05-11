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

#include "Snap.h"
#include "graphs.h"
#include "matio.h"
#include "matstorage.h"
#include "spvb.h"
#include "utils.h"
#include <cassert>

int main(int argc, char *argv[]) {

    if (argc != 3) {
        fprintf(stderr, "Usage: %s [input filename] [output_filename]\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    float *bc_scores, *bc_scores_mod, *partial_bc_scores;
    int *degree;

    matrix_pcsr_t m_csr;
    matrix_pcoo_t m_coo;
    gprops_t gp;

    /*
     * Load matrix in COO format.
     */
    if (query_gprops(argv[1], &gp) || read_matrix(argv[1], &m_coo, &gp)) {
        return EXIT_FAILURE;
    }

    /*
     * Get largest connected component of the undirected graph.
     */
    pcoo_to_pcsr(&m_coo, &m_csr);

    matrix_pcsr_t g;
    components_t cc_array;
    int cc_count = get_cc(&m_csr, &cc_array);
    gp.is_connected = (cc_count == 1);

    print_gprops(&gp);

    if (!gp.is_connected) {
        int max_idx = 0, cmax = 0;
        for (int i = 0; i < cc_count; i++) {
            if (cc_array.cc_size[i] > cmax) {
                max_idx = i;
                cmax = cc_array.cc_size[i];
            }
        }

        extract_und_subgraph(cc_array.array,
                             cc_array.cc_size[max_idx],
                             &m_csr,
                             &g);
    } else {
        g = m_csr;
    }

    /*
     * Load the matrix of the subgraph in CSR as an undirected graph in SNAP.
     *
     * TODO handle conversion from CSR to Snap graph
     */
    //    PUNGraph snap_g = TUNGraph::New();
    //    for(int i = 0; i < m_csr.nrows; i++) {
    //        snap_g->AddNode(i);
    //    }
    //
    //    for(int i = 0; i < m_csr.row_offsets[m_csr.nrows]; i++) {
    //        snap_g->AddEdge(m_csr.rows[i], m_csr.cols[i]);
    //    }

    /*
     * Compute BC using CSR.
     */
    bc_scores = (float *) malloc(m_coo.nrows * sizeof(*bc_scores));
    assert(bc_scores);

    BC_computation(&m_csr, bc_scores, gp.is_directed);

    print_bc_scores(&m_csr, bc_scores, stdout);

    /*
     * Compute BC with SPVB, using TUNGraph.
     */
    //    bc_scores_mod = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
    //    partial_bc_scores = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
    //    assert(bc_scores_mod);
    //    assert(partial_bc_scores);
    //
    //    spvb(snap_g, bc_scores_mod, partial_bc_scores, gp.is_directed);
    //
    //    print_bc_scores(&m_csr, bc_scores_mod, stdout);

    /*
     * Cleanup.
     */
    free_matrix_pcoo(&m_coo);
    free_matrix_pcsr(&m_csr);
    free(bc_scores);
    //    free(bc_scores_mod);
    //    free(partial_bc_scores);
    free_matrix_pcsr(&g);
    free_ccs(&cc_array);

    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);

    return EXIT_SUCCESS;
}
