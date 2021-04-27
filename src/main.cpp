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
 * ---------------------------------------------------------------------------
 *
 * Compile with:
 * g++ -O3 betweenness.cpp -o betweenness
 *
 * Run with:
 * ./betweenness [input_filename] [output_filename]
 *
 * ---------------------------------------------------------------------------
 *
 * TODO Test BC of directed graphs
 *
 ****************************************************************************/

#include <cassert>
#include "utils.h"
#include "matio.h"
#include "matstorage.h"
#include "bc.h"

void extract_subgraph(int* ids, int nids, matrix_csr_t* g, matrix_coo_t* m) {
    const float *p;
    for(int i = 0; i < nids; i++) {
        for(int j = 0; j < g->nrows; j++) {
            spvb(g, p);
        }
    }
}

void DFS_visit(matrix_csr_t* g, bool* visited, int s, int* cc_array, int nids) {

    /*
     * Let the STL vector manage memory allocation.
     */
    std::stack<int> S;
    std::vector<int> subgraph_vertices;
    S.push(s);
    int vertices_cnt = 0;

    while (!S.empty()) {
        s = S.top();
        S.pop();

        if (!visited[s]) {
            visited[s] = true;
            subgraph_vertices.push_back(s);
            vertices_cnt++;
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
    nids = subgraph_vertices.size();
    cc_array = (int*) malloc(nids * sizeof(int));
    assert(cc_array);

    for(int i = 0; i < nids; i++) {
        cc_array[i] = cc_array[i];
    }
}

int get_cc(matrix_csr_t* g, matrix_csr_t** ccs_array) {

    auto visited = (bool*) malloc(g->nrows * sizeof(bool));
    assert(visited);
    int *ids = nullptr, nids = 0, cc_count = 0;

    for(int i = 0; i < g->nrows; i++)
        visited[i] = false;

    for(int i = 0; i < g->nrows; i++) {
        if(!visited[i]) {
            DFS_visit(g, visited, i, ids, nids);
            matrix_coo_t subgraph_coo;
            matrix_csr_t subgraph;

            /*
             * WARNING: too much memory required here. Extract in csr directly.
             */
            extract_subgraph(ids, nids, g, &subgraph_coo);
            coo_to_csr(&subgraph_coo, &subgraph);
            ccs_array[cc_count] = &subgraph;
            print_matrix_csr(&subgraph);

            cc_count++;
            printf(" |\n");
        }
    }

    return cc_count;
}

int main( int argc, char *argv[] ) {

    if (argc != 3) {
        fprintf(stderr, "Usage: %s [input filename] [output_filename]\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    matrix_csr_t m_csr;
    float* bc_scores;
    gprops_t graph_props;
    int *in_degree, *out_degree, *degree;

    /*
     * State graph type.
     */
    graph_props.is_directed = false;
    graph_props.has_self_loops = false;
    graph_props.is_connected = false;
    graph_props.is_weighted = false;

   /*
    * Load matrix in COO format.
    */
    matrix_coo_t m_coo;
    if(read_matrix(argv[1], &m_coo)) {
        fprintf(stderr, "Error reading matrix\n");
        return EXIT_FAILURE;
    }

    /*
     * Convert the internal storage representation of the matrix from COO to
     * CSR.
     *
     * WARNING: Too much memory used, parse into CSR directly.
     */
    coo_to_csr(&m_coo, &m_csr);
    print_matrix_csr(&m_csr);
    print_edge_list(m_csr.row_offsets, m_csr.cols, m_csr.nrows);

    /*
     * Compute degrees.
     */
    if(graph_props.is_directed) {
        in_degree = (int*) malloc(m_coo.nrows * sizeof(int));
        out_degree = (int*) malloc(m_coo.nrows * sizeof(int));
        assert(in_degree);
        assert(out_degree);
        compute_degrees_directed(&m_coo, in_degree, out_degree);
//        print_array(in_degree, m_coo.nrows - 1);
//        print_array(out_degree, m_coo.nrows - 1);

    } else {
        degree = (int*) malloc(m_coo.nrows * sizeof(int));
        assert(degree);
        compute_degrees_undirected(&m_coo, degree);
//        print_array(degree, m_coo.nrows - 1);
    }

    /*
     * Get connected components of undirected graph.
     */
//    matrix_csr_t **cc_list = nullptr;
//    int cc_count = get_cc(&m_csr, cc_list);
//    graph_props.is_connected = (cc_count > 1) ? false : true;
//    printf("cc: %d\n", cc_count);

    /*
     * Compute BC.
     */
    bc_scores = (float*) malloc(m_csr.nrows * sizeof(*bc_scores));
    assert(bc_scores);
    BC_computation(&m_csr, bc_scores, graph_props.is_directed);
//    FILE *fout = fopen(argv[2], "w");
//    print_bc_scores(m_csr, bc_scores, stdout);
//    fclose(fout);

    /*
     * Cleanup.
     */
    free_matrix_coo(&m_coo);
    free_matrix_csr(&m_csr);
    free(bc_scores);
    if(graph_props.is_directed) {
        free(in_degree);
        free(out_degree);
    } else {
        free(degree);
    }

    /*
     * Closing standard streams with error checking.
     */
    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);

    return EXIT_SUCCESS;
}
