/****************************************************************************
 * @file bc_par.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Function to compute the betweenness centrality using the parallel
 * algorithm of the Parallel Boost Graph Library.
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

#include "bc_par.h"

void compute_par_bc_cpu(matrix_pcsr_t *g_tmp, float *bc_cpu) {

    /*
     * Build the Boost graph from the pattern CSR matrix given in input.
     */
    int nvertices = g_tmp->nrows;
    int nnz = g_tmp->row_offsets[nvertices];
    typedef boost::adjacency_list<> graph;
    graph g((unsigned long) nvertices);

    auto rows = (int *) malloc(nnz * sizeof(int));
    expand_row_pointer(nvertices, g_tmp->row_offsets, rows);

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

    boost::brandes_betweenness_centrality(g, centrality_map);

    for (int i = 0; i < nvertices; i++) {
        bc_cpu[i] = (float) (centrality_map[i]);
    }

    /*
     * Count each edge only one time.
     */
    for (int k = 0; k < nvertices; k++)
        bc_cpu[k] /= 2;

}
