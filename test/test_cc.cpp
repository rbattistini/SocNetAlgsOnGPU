/****************************************************************************
 *
 * test_cc.cpp
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

#include <cstdlib>
#include "graphs.h"
#include "tests.h"

static matrix_pcsr_t csr;
static gprops_t gprops;

TEST_CASE(
        "Test connected component algorithm with disconnected undirected graph") {
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
    CHECK_UNARY_FALSE(gprops.is_connected);

    /*
     * Ensure the graph has two connected components.
     */
    CHECK_EQ(cc_count, 2);

    /*
     * Ensure that the first cc is composed by vertices: 0 1 4 2.
     */
    int cc[] = {0, 1, 4, 2};
    for (int i = 0; i < cc_list.cc_size[0]; i++) {
        CHECK_EQ(cc_list.array[i], cc[i]);
    }

    /*
     * Ensure that the second cc is composed by vertices:  3 5 6.
     */
    int cc2[] = {3, 5, 6};
    for (int i = cc_list.cc_size[0] + 1; i < cc_list.cc_size[1]; i++) {
        CHECK_EQ(cc_list.array[i], cc2[i]);
    }
}

TEST_CASE(
        "Test connected component algorithm with connected undirected graph") {

    /*
     * Workspace setup for this test.
     */
    int source_row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
    int source_cols[] = {1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};

    csr.nrows = 9;
    csr.cols = source_cols;
    csr.row_offsets = source_row_offsets;

    components_t cc_list;
    int cc_count = get_cc(&csr, &cc_list);
    gprops.is_connected = (cc_count == 1);

    /*
     * The graph is connected.
     */
    CHECK_UNARY(gprops.is_connected);

    /*
     * The graph has one connected component.
     */
    CHECK_EQ(cc_count, 1);

    std::sort(cc_list.array, cc_list.array + csr.nrows);

    /*
     * The cc is composed by all vertices [0-8].
     */
    for (int i = 0; i < cc_list.cc_size[0]; i++) {
        CHECK_EQ(cc_list.array[i], i);
    }
}

TEST_CASE("Test subgraph extraction from undirected graph given vertices ids") {

    /*
     * Workspace setup for this test.
     */
    matrix_pcsr_t subgraph;
    components_t ccs;
    int nccs = 2, nvertices = 7, size;
    int *vertices;
    int ccs_array[] = {0, 1, 4, 2, 3, 6, 5};

    ccs.array = (int *) malloc(sizeof(*ccs.array) * nvertices);
    ccs.cc_size = (int *) malloc(sizeof(*ccs.cc_size) * nccs);

    for (int i = 0; i < nvertices; i++)
        ccs.array[i] = ccs_array[i];

    ccs.cc_size[0] = 4;
    ccs.cc_size[1] = 3;

    int source_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
    int source_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

    csr.nrows = nvertices;
    csr.cols = source_cols;
    csr.row_offsets = source_row_offsets;

    /*
     * Ensure that the first cc is composed by vertices: 0 1 4 2
     *
     * in CSR this means:
     *
     * nrows = 4
     * R = [ 0 1 4 5 6 ]
     * C = [ 1 0 2 3 1 1 ]
     */
    size = ccs.cc_size[0];
    vertices = (int *) malloc(sizeof(*vertices) * size);

    for (int i = 0; i < size; i++) {
        vertices[i] = ccs.array[i];
    }

    extract_und_subgraph(vertices, size, &csr, &subgraph);

    int nrows = subgraph.nrows;
    int expected_row_offsets[] = {0, 1, 4, 5, 6};
    int ncols = 6;
    int expected_cols[] = {1, 0, 3, 2, 1, 1};

    REQUIRE_EQ(nrows, 4);

    for (int i = 0; i < ncols; i++) {
        CHECK_EQ(subgraph.cols[i], expected_cols[i]);
    }

    for (int i = 0; i < nrows; i++) {
        CHECK_EQ(subgraph.row_offsets[i], expected_row_offsets[i]);
    }

    free(vertices);
    free_matrix_pcsr(&subgraph);

    /*
     * Ensure that the second cc is composed by vertices:  3 5 6
     *
     * in CSR this means:
     *
     * nrows = 3
     * R = [ 0 2 4 6 ]
     * C = [ 1 2 0 2 0 1 ]
     */
    size = ccs.cc_size[1];
    vertices = (int *) malloc(sizeof(*vertices) * size);
    int start = ccs.cc_size[0];
    int end = ccs.cc_size[1] + ccs.cc_size[0];

    for (int i = start, j = 0; i < end; j++, i++) {
        vertices[j] = ccs.array[i];
    }

    extract_und_subgraph(vertices, size, &csr, &subgraph);

    nrows = subgraph.nrows;
    int expected_row_offsets2[] = {0, 2, 4, 6};
    ncols = 6;
    int expected_cols2[] = {2, 1, 0, 2, 0, 1};

    REQUIRE_EQ(nrows, size);

    for (int i = 0; i < ncols; i++) {
        CHECK_EQ(subgraph.cols[i], expected_cols2[i]);
    }

    for (int i = 0; i < nrows; i++) {
        CHECK_EQ(subgraph.row_offsets[i], expected_row_offsets2[i]);
    }

    free(vertices);
    free_matrix_pcsr(&subgraph);
}
