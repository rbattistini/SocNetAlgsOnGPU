/****************************************************************************
 * @file test_cc.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
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

#include "graphs.h"
#include "tests.h"

static matrix_pcsr_t A;
static gprops_t gprops;

TEST_CASE(
        "Test connected component algorithm with disconnected undirected graph") {
    /*
     * Workspace setup for this test.
     */
    int source_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
    int source_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

    A.nrows = 7;
    A.cols = source_cols;
    A.row_offsets = source_row_offsets;

    components_t ccs;
    get_cc(&A, &ccs);
    gprops.is_connected = (ccs.cc_count == 1);

    /*
     * Ensure the graph is not connected.
     */
    CHECK_UNARY_FALSE(gprops.is_connected);

    /*
     * Ensure the graph has two connected components.
     */
    CHECK_EQ(ccs.cc_count, 2);

    /*
     * Ensure that the first cc is composed by vertices: 0 1 4 2.
     */
    int cc[] = {0, 1, 4, 2};
    for (int i = 0; i < ccs.cc_size[0]; i++) {
        CHECK_EQ(ccs.array[i], cc[i]);
    }

    /*
     * Ensure that the second cc is composed by vertices:  3 5 6.
     */
    int cc2[] = {3, 5, 6};
    for (int i = ccs.cc_size[0] + 1; i < ccs.cc_size[1]; i++) {
        CHECK_EQ(ccs.array[i], cc2[i]);
    }
}

TEST_CASE(
        "Test connected component algorithm with connected undirected graph") {

    /*
     * Workspace setup for this test.
     */
    int source_row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
    int source_cols[] = {1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};

    A.nrows = 9;
    A.cols = source_cols;
    A.row_offsets = source_row_offsets;

    components_t ccs;
    get_cc(&A, &ccs);
    gprops.is_connected = (ccs.cc_count == 1);

    /*
     * The graph is connected.
     */
    CHECK_UNARY(gprops.is_connected);

    /*
     * The graph has one connected component.
     */
    CHECK_EQ(ccs.cc_count, 1);

    std::sort(ccs.array, ccs.array + A.nrows);

    /*
     * The cc is composed by all vertices [0-8].
     */
    for (int i = 0; i < ccs.cc_size[0]; i++) {
        CHECK_EQ(ccs.array[i], i);
    }
}

TEST_CASE("Test subgraph extraction from undirected graph given vertices ids") {

    /*
     * Workspace setup for this test.
     */
    matrix_pcsr_t R, Q, C;
    components_t ccs;
    int nvertices = 7, size;
    int *vertices;
    int ccs_array[] = {0, 1, 4, 2, 3, 6, 5};

    ccs.cc_count = 2;
    ccs.array = (int *) malloc(sizeof(*ccs.array) * nvertices);
    ccs.cc_size = (int *) malloc(sizeof(*ccs.cc_size) * ccs.cc_count);

    for (int i = 0; i < nvertices; i++)
        ccs.array[i] = ccs_array[i];

    ccs.cc_size[0] = 4;
    ccs.cc_size[1] = 3;

    int source_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
    int source_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

    A.nrows = nvertices;
    A.ncols = nvertices;
    A.cols = source_cols;
    A.row_offsets = source_row_offsets;

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

    std::sort(vertices, vertices + size);
    extract_subgraph(vertices, size, &A, &C);

    int nrows = C.nrows;
    int expected_row_offsets[] = {0, 1, 4, 5, 6};
    int nnz = C.row_offsets[C.nrows];
    int expected_cols[] = {1, 0, 2, 3, 1, 1};

    REQUIRE_EQ(nrows, 4);

    for (int i = 0; i < nnz; i++) {
        CHECK_EQ(C.cols[i], expected_cols[i]);
    }

    for (int i = 0; i < nrows; i++) {
        CHECK_EQ(C.row_offsets[i], expected_row_offsets[i]);
    }

    free(vertices);
    free_matrix_pcsr(&C);

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

    std::sort(vertices, vertices + size);
    extract_subgraph(vertices, size, &A, &C);

    nrows = C.nrows;
    int expected_row_offsets2[] = {0, 2, 4, 6};
    nnz = C.row_offsets[C.nrows];
    int expected_cols2[] = {1, 2, 0, 2, 0, 1};

    print_matrix_pcsr(&C);

    REQUIRE_EQ(nrows, size);

    for (int i = 0; i < nnz; i++) {
        CHECK_EQ(C.cols[i], expected_cols2[i]);
    }

    for (int i = 0; i < nrows; i++) {
        CHECK_EQ(C.row_offsets[i], expected_row_offsets2[i]);
    }

    free(vertices);
    free_matrix_pcsr(&C);
}

TEST_CASE("Test largest connected component extraction of undirected unweighted graph") {

    SUBCASE("Test extraction of the largest cc when it is the first") {
        /*
         * Workspace setup for this test.
         */
        matrix_pcsr_t subgraph;
        components_t ccs;
        int nvertices = 7;
        int ccs_array[] = {0, 1, 4, 2, 3, 6, 5};

        ccs.cc_count = 2;
        ccs.array = (int *) malloc(sizeof(*ccs.array) * nvertices);
        ccs.cc_size = (int *) malloc(sizeof(*ccs.cc_size) * ccs.cc_count);

        ccs.cc_size[0] = 4;
        ccs.cc_size[1] = 3;

        for (int i = 0; i < nvertices; i++)
            ccs.array[i] = ccs_array[i];

        int source_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
        int source_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

        A.nrows = nvertices;
        A.ncols = nvertices;
        A.row_offsets = source_row_offsets;
        A.cols = source_cols;

        /*
         * Ensure the cc extracted is the largest and is composed by
         * vertices 0 1 4 2.
         *
         * in CSR this means:
         *
         * nrows = 4
         * R = [ 0 1 4 5 6 ]
         * C = [ 1 0 2 3 1 1 ]
         */
        get_largest_cc(&A, &subgraph, &ccs);

        print_matrix_pcsr(&subgraph);

        int nrows = subgraph.nrows;
        int expected_row_offsets[] = {0, 1, 4, 5, 6};
        int ncols = subgraph.ncols;
        int expected_cols[] = {1, 0, 2, 3, 1, 1};

        REQUIRE_EQ(nrows, 4);
        REQUIRE_EQ(ncols, 4);

        for (int i = 0; i < ncols; i++) {
            CHECK_EQ(subgraph.cols[i], expected_cols[i]);
        }

        for (int i = 0; i < nrows; i++) {
            CHECK_EQ(subgraph.row_offsets[i], expected_row_offsets[i]);
        }

        free_matrix_pcsr(&subgraph);
    }

    SUBCASE("Test connection of the largest cc when it is not the first") {

        /*
         * Workspace setup for this test.
         */
        matrix_pcsr_t subgraph;
        components_t ccs;
        int nvertices = 15;
        int ccs_array[] =
                {0, 1, 4, 2, 3, 6, 5, 7, 8, 9, 11, 10, 12, 13, 14};

        ccs.cc_count = 4;
        ccs.array = (int *) malloc(sizeof(*ccs.array) * nvertices);
        ccs.cc_size = (int *) malloc(sizeof(*ccs.cc_size) * ccs.cc_count);

        ccs.cc_size[0] = 4;
        ccs.cc_size[1] = 3;
        ccs.cc_size[2] = 5;
        ccs.cc_size[3] = 3;

        for (int i = 0; i < nvertices; i++)
            ccs.array[i] = ccs_array[i];

        int source_row_offsets[] =
                {0, 1, 4, 5, 7, 8, 10, 12, 13, 15, 18, 19, 20, 21, 23, 24 };
        int source_cols[] =
                {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5, 8, 7, 9, 8, 10, 11, 9, 9, 13, 12, 14, 13};

        A.nrows = nvertices;
        A.ncols = nvertices;
        A.cols = source_cols;
        A.row_offsets = source_row_offsets;

        /*
         * Ensure the cc extracted is the largest and is composed by
         * vertices 7 8 9 10 11.
         *
         * in CSR this means:
         *
         * nrows = 5
         * R = [ 0 1 3 6 7 8 ]
         * C = [ 1 0 2 1 4 3 2 2 ]
         */
        get_largest_cc(&A, &subgraph, &ccs);

        int nrows = subgraph.nrows;
        int expected_row_offsets[] = {0, 1, 3, 6, 7, 8};
        int ncols = subgraph.ncols;
        int expected_cols[] = {1, 0, 2, 1, 3, 4, 2, 2};

        REQUIRE_EQ(nrows, 5);

        for (int i = 0; i < ncols; i++) {
            CHECK_EQ(subgraph.cols[i], expected_cols[i]);
        }

        for (int i = 0; i < nrows; i++) {
            CHECK_EQ(subgraph.row_offsets[i], expected_row_offsets[i]);
        }

        free_matrix_pcsr(&subgraph);
    }
}
