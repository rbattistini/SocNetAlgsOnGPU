/****************************************************************************
 * @file matstorage.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions to convert from the COOrdinate format to the Compressed
 * Sparse Rows. In addition it provides a way of representing adjacency
 * matrices with structures of arrays.
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

#include "matstorage.h"

int check_matrix_pcoo(matrix_pcoo_t *matrix) {
    return (matrix->cols != 0) &&
           (matrix->rows != 0) &&
           (matrix->nrows > -1) &&
           (matrix->nnz > -1);
}

int check_matrix_rcoo(matrix_rcoo_t *matrix) {
    return (matrix->cols != 0) &&
           (matrix->rows != 0) &&
           (matrix->weights != 0) &&
           (matrix->nrows > -1) &&
           (matrix->nnz > -1);
}

int check_matrix_pcsr(matrix_pcsr_t *matrix) {
    return (matrix->cols != 0) &&
           (matrix->row_offsets != 0) &&
           (matrix->nrows > -1);
}

int check_matrix_rcsr(matrix_rcsr_t *matrix) {
    return (matrix->cols != 0) &&
           (matrix->row_offsets != 0) &&
           (matrix->weights != 0) &&
           (matrix->nrows > -1);
}

void free_matrix_pcoo(matrix_pcoo_t *matrix) {
    free(matrix->rows);
    free(matrix->cols);
    matrix->rows = 0;
    matrix->cols = 0;
    matrix->nnz = -1;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void free_matrix_pcsr(matrix_pcsr_t *matrix) {
    free(matrix->row_offsets);
    free(matrix->cols);
    matrix->row_offsets = 0;
    matrix->cols = 0;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void free_matrix_rcoo(matrix_rcoo_t *matrix) {
    free(matrix->rows);
    free(matrix->cols);
    free(matrix->weights);
    matrix->rows = 0;
    matrix->cols = 0;
    matrix->nnz = -1;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void free_matrix_rcsr(matrix_rcsr_t *matrix) {
    free(matrix->row_offsets);
    free(matrix->cols);
    free(matrix->weights);
    matrix->row_offsets = 0;
    matrix->cols = 0;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void print_matrix_pcoo(matrix_pcoo_t *matrix) {

    if (!check_matrix_pcoo(matrix)) {
        ZF_LOGE("The matrix is not initialized");
        return;
    }

    printf("nrows = %d\n", matrix->nrows);
    printf("ncols = %d\n", matrix->ncols);
    printf("nnz = %d\n", matrix->nnz);
    printf("rows = \n");
    print_int_array(matrix->rows, matrix->nnz - 1);
    printf("cols = \n");
    print_int_array(matrix->cols, matrix->nnz - 1);
}

void print_matrix_pcsr(matrix_pcsr_t *matrix) {

    if (!check_matrix_pcsr(matrix)) {
        ZF_LOGE("The matrix is not initialized");
        return;
    }

    int nnz = matrix->row_offsets[matrix->nrows];
    printf("nrows = %d\n", matrix->nrows);
    printf("ncols = %d\n", matrix->ncols);
    printf("offsets = \n");
    print_int_array(matrix->row_offsets, matrix->nrows);
    printf("cols = \n");
    print_int_array(matrix->cols, nnz - 1);
}

void expand_row_pointer(int nrows, const int *row_offsets, int *rows) {

    for (int i = 0; i < nrows; i++) {
        for (int j = row_offsets[i]; j < row_offsets[i + 1]; j++) {
            rows[j] = i;
        }
    }
}

int coo_to_csr(matrix_pcoo_t *A, matrix_pcsr_t *B) {

    if (!check_matrix_pcoo(A)) {
        ZF_LOGE("The matrix is not initialized");
        return EXIT_FAILURE;
    }

    int *row_offsets, *cols;
    int *rows = A->rows; // row indices of A
    int nnz = A->nnz;    // number of nnz in A
    int nrows = A->nrows;// number of rows in A
    int ncols = A->ncols;

    row_offsets = (int *) calloc((nrows + 1), sizeof(*row_offsets));
    cols = (int *) malloc(nnz * sizeof(*cols));

    if (row_offsets == 0 || cols == 0) {
        ZF_LOGF("Memory allocation failed!");
        return EXIT_FAILURE;
    }

    /*
     * Compute number of non-zero entries per row of A.
     */
    for (int n = 0; n < nnz; n++) {
        row_offsets[rows[n]]++;
    }

    /*
     * Compute row offsets
     */
    for (int i = 0, psum = 0; i < nrows; i++) {
        int temp = row_offsets[i];
        row_offsets[i] = psum;
        psum += temp;
    }
    row_offsets[nrows] = nnz;

    /*
     * Copy cols array of A in cols of B
     */
    for (int n = 0; n < nnz; n++) {
        int row = rows[n];
        int dest = row_offsets[row];
        cols[dest] = A->cols[n];
        row_offsets[row]++;
    }

    for (int i = 0, last = 0; i <= nrows; i++) {
        int temp = row_offsets[i];
        row_offsets[i] = last;
        last = temp;
    }

    B->nrows = nrows;
    B->ncols = ncols;
    B->cols = cols;
    B->row_offsets = row_offsets;

    return EXIT_SUCCESS;
}

int transpose(matrix_pcsr_t *A, matrix_pcsr_t *B) {

    if (!check_matrix_pcsr(A)) {
        ZF_LOGE("The matrix is not initialized");
        return EXIT_FAILURE;
    }

    int *row_offsets, *cols;
    int nrows = A->nrows;
    int ncols = A->ncols;
    int nnz = A->row_offsets[A->nrows];

    row_offsets = (int *) calloc((ncols + 1), sizeof(*row_offsets));
    cols = (int *) malloc(nnz * sizeof(*cols));

    if (row_offsets == 0 || cols == 0) {
        ZF_LOGF("Memory allocation failed!");
        return EXIT_FAILURE;
    }

    /*
     * Compute number of non-zero entries per column of A.
     */
    for (int n = 0; n < nnz; n++) {
        row_offsets[A->cols[n]]++;
    }

    /*
     * Compute row offsets (column offsets of A)
     */
    for (int i = 0, psum = 0; i < ncols; i++) {
        int temp = row_offsets[i];
        row_offsets[i] = psum;
        psum += temp;
    }
    row_offsets[ncols] = nnz;

    /*
     * Copy cols array of A in cols of B
     */
    for (int n = 0; n < nrows; n++) {
        for (int j_a = A->row_offsets[n]; j_a < A->row_offsets[n + 1]; j_a++) {
            int col = A->cols[j_a];
            int dest = row_offsets[col];
            cols[dest] = n;
            row_offsets[col]++;
        }
    }

    for (int i = 0, last = 0; i <= ncols; i++) {
        int temp = row_offsets[i];
        row_offsets[i] = last;
        last = temp;
    }

    B->nrows = ncols;
    B->ncols = nrows;
    B->row_offsets = row_offsets;
    B->cols = cols;

    return EXIT_SUCCESS;
}
