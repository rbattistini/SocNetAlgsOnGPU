/****************************************************************************
 *
 * test_matrix_io.cpp
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
#include "utils.h"
#include "graphs.h"
#include "matio.h"
#include "matstorage.h"
#include "tests.h"

static matrix_pcoo_t coo;
static matrix_pcsr_t csr;
static gprops_t gprops;

static void clear_workspace() {

    /*
     * Closing standard streams with error checking.
     */
    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);
}

TEST_CASE("Test graph loading with properties querying") {

    SUBCASE("undirected unweighted graph") {
        const char *gname = "und_unw_test.mtx";

        query_gprops(gname, &gprops);

        CHECK(gprops.is_directed == false);
        CHECK(gprops.is_weighted == false);
    }

    SUBCASE("undirected weighted graph") {
        const char *gname = "und_wgh_test.mtx";

        query_gprops(gname, &gprops);

        CHECK(gprops.is_directed == false);
        CHECK(gprops.is_weighted == true);
    }

    SUBCASE("directed unweighted graph") {
        const char *gname = "dir_unw_test.mtx";

        query_gprops(gname, &gprops);

        CHECK(gprops.is_directed == true);
        CHECK(gprops.is_weighted == false);
    }

    SUBCASE("directed weighted graph") {
        const char *gname = "dir_wgh_test.mtx";

        query_gprops(gname, &gprops);

        CHECK(gprops.is_directed == true);
        CHECK(gprops.is_weighted == true);
    }
}

TEST_CASE("Test graph loading using mmio with some corner cases") {

    SUBCASE("edges missing for undirected graph without symmetric property") {
        const char *digraph = "error1.mtx";

        CHECK(read_matrix(digraph, &coo, &gprops));
        CHECK(coo.nnz != 18);
    }

    SUBCASE("check out of bound indices") {
        const char *gname = "error2.mtx";

        CHECK_FALSE(read_matrix(gname, &coo, &gprops));
    }

    SUBCASE("handling numbers near infinite") {
        const char *gname = "error3.mtx";

        CHECK_FALSE(read_matrix(gname, &coo, &gprops));
    }

    SUBCASE("handling non-numeric values") {
        const char *gname = "error4.mtx";

        CHECK_FALSE(read_matrix(gname, &coo, &gprops));
    }

    clear_workspace();
}

TEST_CASE("Test sparse matrix storage format conversion from COO to CSR") {

    /*
     * Workspace setup for this test.
     */
    gprops.is_directed = false;
    int rows[] = {1, 0, 3, 0, 4, 0, 5, 0, 2, 1, 6, 2, 7, 2, 5, 4, 8, 7};
    int cols[] = {0, 1, 0, 3, 0, 4, 0, 5, 1, 2, 2, 6, 2, 7, 4, 5, 7, 8};

    coo.nnz = 18;
    coo.nrows = 9;
    coo.rows = rows;
    coo.cols = cols;

    pcoo_to_pcsr(&coo, &csr);

    int expected_row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
    int expected_cols[] =
            { 1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};

    CHECK(csr.nrows == 9);

    for(int i = 0; i < coo.nrows; i++) {
        CHECK(csr.row_offsets[i] == expected_row_offsets[i]);
    }

    for(int i = 0; i < coo.nnz; i++) {
        CHECK(csr.cols[i] == expected_cols[i]);
    }
}
