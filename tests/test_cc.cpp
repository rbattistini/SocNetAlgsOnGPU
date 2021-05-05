/****************************************************************************
 *
 * test_connectivity.cpp
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

#include <cstdio>
#include <cstdlib>
#include "utils.h"
#include "graphs.h"
#include "tests.h"

static matrix_pcsr_t csr;
static matrix_pcoo_t coo;
static gprops_t gprops;

typedef matrix_pcsr_t graph;

static int clear_workspace() {

    /*
     * Closing standard streams with error checking.
     */
    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);

    return EXIT_SUCCESS;
}

TEST_CASE("Test connected component algorithm with disconnected undirected graph") {

    /*
     * Workspace setup for this test.
     */
    int source_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
    int source_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

    csr.nrows = 7;
    csr.cols = source_cols;
    csr.row_offsets = source_row_offsets;

    components_t cc_list;
    int cc_count = get_cc(&csr, &cc_list);
    gprops.is_connected = (cc_count == 1);

    /*
     * Ensure the graph is not connected.
     */
    CHECK_FALSE(gprops.is_connected);

    /*
     * Ensure the graph has two connected components.
     */
    CHECK(cc_count == 2);

    /*
     * Ensure that the first cc is composed by vertices: 0 1 4 2
     *
     * in CSR this means:
     *
     * nrows = 4
     * R = [ 0 1 4 5 6 ]
     * C = [ 1 0 2 3 1 1 ]
     */
    graph subgraph;
    size_t size = cc_list.ccs_size[0];
    auto *vertices = (vertex*) malloc(sizeof(vertex) * size);

    for(int i = 0; i < size; i++) {
        vertices[i] = cc_list.ccs_array[i];
    }

    extract_und_subgraph(cc_list.ccs_array, cc_list.ccs_size[0],
                         &csr, &subgraph);

    int nrows = subgraph.nrows;
    int row_offsets[] = {0, 1, 4, 5, 6};
    int ncols = 6;
    int cols[] = {1, 0, 2, 3, 1, 1};

    CHECK(nrows == 4);

    for(int i = 0; i < ncols; i++) {
        CHECK(subgraph.cols[i] == cols[i]);
    }

    for(int i = 0; i < nrows; i++) {
        CHECK(subgraph.row_offsets[i] == row_offsets[i]);
    }

    /*
     * Ensure that the second cc is composed by vertices:  3 5 6
     *
     * in CSR this means:
     *
     * nrows = 3
     * R = [ 0 2 4 6 ]
     * C = [ 1 2 0 2 0 1 ]
     */
    free_matrix_csr(&subgraph);
    free(vertices);
    size = cc_list.ccs_size[1];
    vertices = (vertex*) malloc(sizeof(vertex) * size);

    for(int i = cc_list.ccs_size[1] - cc_list.ccs_size[0], j = 0; i < size; j++, i++) {
        vertices[j] = cc_list.ccs_array[i];
    }

    extract_und_subgraph(vertices, size, &csr, &subgraph);

    nrows = subgraph.nrows;
    int row_offsets2[] = {0, 2, 4, 6};
    ncols = 6;
    int cols2[] = {1, 2, 0, 2, 0, 1};

    CHECK(nrows == 3);

    for(int i = 0; i < ncols; i++) {
        CHECK(subgraph.cols[i] == cols2[i]);
    }

    for(int i = 0; i < nrows; i++) {
        CHECK(subgraph.row_offsets[i] == row_offsets2[i]);
    }

    free(vertices);
    clear_workspace();
}

//
//TEST_CASE("Test connected component algorithm with connected undirected graph") {
//
//    /*
//     * Workspace setup for this test.
//     */
//    int source_row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
//    int source_cols[] = { 1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};
//
//    csr->nrows = 9;
//    csr->cols = source_cols;
//    csr->row_offsets = source_row_offsets;
//
//    components_t *cc_list = nullptr;
//    int cc_count = get_cc(csr, cc_list);
//    gprops.is_connected = (cc_count == 1);
//
//    /*
//     * The graph is connected.
//     */
//    CHECK(gprops.is_connected);
//
//    /*
//     * The graph has one connected component.
//     */
//    REQUIRE(cc_list);
//    CHECK(cc_count == 1);
//
//    /*
//     * The cc is composed by all vertices.
//     *
//     * in CSR this means:
//     *
//     * nrows = 9
//     * R = [ 0 4 6 9 10 12 14 15 17 18 ]
//     * C = [ 1 3 4 5 0 2 1 6 7 0 0 5 0 4 2 2 8 7 ]
//     */
//    graph *subgraph;
//    extract_und_subgraph(cc_list->ccs_array, cc_list->ccs_size[0],
//                         csr, &subgraph);
//
//    int nrows = subgraph->nrows;
//    int row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
//    int ncols = 18;
//    int cols[] = {1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};
//
//    CHECK(nrows == 9);
//
//    for(int i = 0; i < ncols; i++) {
//        CHECK(subgraph->cols[i] == cols[i]);
//    }
//
//    for(int i = 0; i < nrows; i++) {
//        CHECK(subgraph->row_offsets[i] == row_offsets[i]);
//    }
//
//    clear_workspace();
//}
