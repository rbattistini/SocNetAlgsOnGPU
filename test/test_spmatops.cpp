/****************************************************************************
 * @file test_spmatops.cpp
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

#include "spmatops.h"
#include "tests.h"
#include <matds.h>

static matrix_pcsr_t A;
static matrix_pcsr_t B;
static matrix_pcsr_t C;
static matrix_pcsr_t R;
static matrix_pcsr_t Q;

static matrix_pcoo_t coo;
static matrix_pcsr_t csr;

TEST_CASE("Test sparse matrix storage format conversion from COO to CSR") {

    /*
     * Workspace setup for this test.
     */
    int rows[] = {1, 0, 3, 0, 4, 0, 5, 0, 2, 1, 6, 2, 7, 2, 5, 4, 8, 7};
    int cols[] = {0, 1, 0, 3, 0, 4, 0, 5, 1, 2, 2, 6, 2, 7, 4, 5, 7, 8};

    coo.nnz = 18;
    coo.nrows = 9;
    coo.rows = rows;
    coo.cols = cols;

    coo_to_csr(&coo, &csr);

    int expected_row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
    int expected_cols[] =
            {1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};

    CHECK_EQ(csr.nrows, 9);

    for (int i = 0; i < coo.nrows; i++) {
        CHECK_EQ(csr.row_offsets[i], expected_row_offsets[i]);
    }

    for (int i = 0; i < coo.nnz; i++) {
        CHECK_EQ(csr.cols[i], expected_cols[i]);
    }
}

TEST_CASE("Test sparse matrix in CSR format transposition") {

    /*
     * Matrix R (4 x 7)
     */
    int r_row_offsets[] = {0, 1, 2, 3, 4};
    int r_cols[] = {0, 1, 2, 4};

    R.nrows = 4;
    R.ncols = 7;
    R.cols = r_cols;
    R.row_offsets = r_row_offsets;

    /*
     * Matrix Q, transpose of R (7 x 4)
     */
    int q_row_offsets[] = {0, 1, 2, 3, 3, 4, 4, 4};
    int q_cols[] = {0, 1, 2, 3};

    transpose(&R, &Q);

    /*
     * Ensure that transposition is computed correctly.
     */
    REQUIRE_EQ(Q.nrows, R.ncols);
    REQUIRE_EQ(Q.ncols, R.nrows);

    for (int i = 0; i < Q.nrows; i++) {
        CHECK_EQ(Q.row_offsets[i], q_row_offsets[i]);
    }

    for (int i = 0; i < Q.row_offsets[Q.nrows]; i++) {
        CHECK_EQ(Q.cols[i], q_cols[i]);
    }
}

TEST_CASE("Test spgemm on two pattern matrices") {

    /*******************
     * Workspace setup.
     *******************/

    /*
     * Matrix R (4 x 7)
     */
    int r_row_offsets[] = /* {0, 1, 2, 3, 4} */ {0, 1, 2, 3};
    int r_cols[] = /* {0, 1, 2, 4} */ {3, 5, 6};

    R.nrows = 3;
    R.ncols = 7;
    R.cols = r_cols;
    R.row_offsets = r_row_offsets;

    /*
     * Matrix A (7 x 7)
     */
    int a_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
    int a_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

    A.nrows = 7;
    A.ncols = 7;
    A.cols = a_cols;
    A.row_offsets = a_row_offsets;

    /*
     * B = RA
     */
    spgemm(&R, &A, &B);

    /*
     * Ensure the product is correct.
     */
    REQUIRE_EQ(B.nrows, R.nrows);
    REQUIRE_EQ(B.ncols, A.ncols);

    int expected_row_offsets[] = /* {0, 1, 4, 5, 6} */ {0, 2, 4, 6};
    int expected_cols[] = /* {1, 0, 2, 4, 1, 1} */ {5, 6, 3, 6, 3, 5};

    for (int i = 0; i < B.nrows; i++) {
        CHECK_EQ(B.row_offsets[i], expected_row_offsets[i]);
    }

    for (int i = 0; i < B.row_offsets[B.nrows]; i++) {
        CHECK_EQ(B.cols[i], expected_cols[i]);
    }
}

TEST_CASE("Test spref on three pattern matrices") {

    /*******************
     * Workspace setup.
     *******************/

    /*
     * Matrix R (4 x 7)
     */
    int r_row_offsets[] = {0, 1, 2, 3, 4};
    int r_cols[] = {0, 1, 2, 4};

    R.nrows = 4;
    R.ncols = 7;
    R.cols = r_cols;
    R.row_offsets = r_row_offsets;

    /*
     * Matrix A (7 x 7)
     */
    int a_row_offsets[] = {0, 1, 4, 5, 7, 8, 10, 12};
    int a_cols[] = {1, 0, 2, 4, 1, 5, 6, 1, 3, 6, 3, 5};

    A.nrows = 7;
    A.ncols = 7;
    A.cols = a_cols;
    A.row_offsets = a_row_offsets;

    /*
     * Matrix Q, transpose of R (7 x 4)
     */
    int q_row_offsets[] = {0, 1, 2, 3, 3, 4, 4, 4};
    int q_cols[] = {0, 1, 2, 3};

    Q.nrows = 7;
    Q.ncols = 4;
    Q.cols = q_cols;
    Q.row_offsets = q_row_offsets;

    spref(&R, &A, &Q, &C);

    /*
     * Ensure the triple product is correct.
     */
    REQUIRE_EQ(C.nrows, R.nrows);
    REQUIRE_EQ(C.ncols, Q.ncols);

    int expected_row_offsets[] = {0, 1, 4, 5, 6};
    int expected_cols[] = {1, 0, 2, 3, 1, 1};

    for (int i = 0; i < C.nrows; i++) {
        CHECK_EQ(C.row_offsets[i], expected_row_offsets[i]);
    }

    for (int i = 0; i < C.row_offsets[C.nrows]; i++) {
        CHECK_EQ(C.cols[i], expected_cols[i]);
    }
}
