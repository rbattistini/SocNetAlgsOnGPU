/****************************************************************************
 * @file spmatops.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions for handling sparse matrix-matrix multiplication and
 * generalized indexing into a sparse matrix.
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

#include "spmatops.h"

int spgemm(matrix_pcsr_t *A, matrix_pcsr_t *B, matrix_pcsr_t *C) {

    /*
     * Check inputs.
     */
    if (!check_matrix_pcsr(A) &&
        !check_matrix_pcsr(B) &&
        A->ncols != B->nrows) {
        ZF_LOGF("Input matrices not initialized or with mismatched shapes");
        return EXIT_FAILURE;
    }

    /*
     * Load A matrix.
     */
    int a_nrows = A->nrows;
    int *a_row_offsets = A->row_offsets;
    int *a_cols = A->cols;

    /*
     * Load B matrix.
     */
    int b_ncols = B->ncols;
    int *b_row_offsets = B->row_offsets;
    int *b_cols = B->cols;

    /*
     * Get size of C.
     */
    int c_nrows = a_nrows;
    int c_ncols = b_ncols;

    /*
     * Count nnz of C.
     */
    int c_nnz = 0;
    int *flag = (int *) malloc(c_ncols * sizeof(*flag));

    /*
     * Clear the flag array at each iteration.
     */
    fill(flag, c_ncols, -1);

    int *c_row_offsets = (int *) malloc((c_nrows + 1) * sizeof(int));

    if (c_row_offsets == 0) {
        ZF_LOGF("Memory allocation failed!");
        return EXIT_FAILURE;
    }

    for (int ic = 0; ic < c_nrows; ic++) {

        c_row_offsets[ic] = c_nnz;

        /*
         * For each row of A...
         */
        for (int ia = a_row_offsets[ic]; ia < a_row_offsets[ic + 1]; ia++) {
            /*
             * For each row of B...
             */
            for (int ib = b_row_offsets[a_cols[ia]]; ib < b_row_offsets[a_cols[ia] + 1]; ib++) {
                int j_b = b_cols[ib];

                if (flag[j_b] != ic) {
                    c_nnz++;
                    flag[j_b] = ic;
                }
            }
        }

        /*
         * Check for integer overflow.
         */
        if (c_nnz < 0) {
            ZF_LOGF("Integer overflow occurred!");
            return EXIT_FAILURE;
        }
    }

    /*
     * Allocate C matrix.
     */
    int *c_cols = (int *) malloc(c_nnz * sizeof(int));
    if (c_cols == 0) {
        ZF_LOGF("Memory allocation failed!");
        return EXIT_FAILURE;
    }

    /*******************
     * Perform C = AB.
     *******************/

    c_nnz = 0;

    /*
     * Clear the flag array
     */
    fill(flag, c_ncols, -1);

    /*
     * For each row of C...
     */
    int row_start;
    for (int ic = 0; ic < c_nrows; ic++) {

        row_start = c_row_offsets[ic];

        /*
         * For each row of A...
         */
        for (int k = a_row_offsets[ic]; k < a_row_offsets[ic + 1]; k++) {
            int t = a_cols[k];

            /*
             * For each row of B...
             */
            for (int v = b_row_offsets[t]; v < b_row_offsets[t + 1]; v++) {
                int j = b_cols[v];

                if (flag[j] < row_start) {
                    flag[j] = c_nnz;
                    c_cols[c_nnz] = j;
                    c_nnz++;
                }
            }
        }
    }
    c_row_offsets[c_nrows] = c_nnz;

    C->nrows = c_nrows;
    C->ncols = c_ncols;
    C->row_offsets = c_row_offsets;
    C->cols = c_cols;

    return EXIT_SUCCESS;
}

int spref(matrix_pcsr_t *R,
          matrix_pcsr_t *A,
          matrix_pcsr_t *Q,
          matrix_pcsr_t *C) {

    if (!check_matrix_pcsr(R) ||
        !check_matrix_pcsr(A) ||
        !check_matrix_pcsr(Q)) {
        ZF_LOGF("Input matrices not initialized or mismatched");
        return EXIT_FAILURE;
    }

    matrix_pcsr_t B;

    spgemm(R, A, &B);

    spgemm(&B, Q, C);

    free_matrix_pcsr(&B);

    return EXIT_SUCCESS;
}

int get_R_matrix(matrix_pcsr_t *R,
                 const int *vertices,
                 int nvertices,
                 int nrows) {

    if (vertices == 0 || nvertices <= 0 || nrows <= 0) {
        ZF_LOGF("Input values not valid");
        return EXIT_FAILURE;
    }

    int *rows =
            (int *) malloc(nrows * sizeof(*rows));
    int *cols =
            (int *) malloc(nvertices * sizeof(*cols));

    if (rows == 0 || cols == 0) {
        ZF_LOGF("Memory allocation failed!");
        return EXIT_FAILURE;
    }

    for (int k = 0; k < nvertices; k++) {
        rows[k] = k;
        cols[k] = vertices[k];
    }

    int *row_offsets =
            (int *) calloc((nvertices + 1), sizeof(*row_offsets));

    if (row_offsets == 0) {
        ZF_LOGF("Memory allocation failed!");
        return EXIT_FAILURE;
    }

    /*
     * Compute number of non-zero entries per column of R.
     */
    for (int n = 0; n < nvertices; n++) {
        row_offsets[rows[n]]++;
    }

    /*
     * Compute row offsets
     */
    for (int i = 0, psum = 0; i < nvertices; i++) {
        int temp = row_offsets[i];
        row_offsets[i] = psum;
        psum += temp;
    }
    row_offsets[nvertices] = nvertices;

    free(rows);

    R->nrows = nvertices;
    R->ncols = nrows;
    R->row_offsets = row_offsets;
    R->cols = cols;

    return EXIT_SUCCESS;
}
