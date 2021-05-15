/****************************************************************************
 * @file matstorage.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * Functions to convert from the COOrdinate format to the Compressed Sparse
 * Rows. In addition it provides a way of representing adjacency matrices
 * with structures of arrays.
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

void check_bc(matrix_pcsr_t g, const float *bc_cpu, const float *bc_gpu) {
    for (int i = 0; i < g.row_offsets[g.nrows]; i++) {
        assert(bc_cpu[i] == bc_gpu[i]);
    }
}

int check_matrix_init(matrix_pcoo_t *matrix) {
    return (matrix->cols != nullptr) &&
           (matrix->rows != nullptr) &&
           (matrix->nrows > -1) &&
           (matrix->nnz > -1);
}

int check_matrix_init(matrix_rcoo_t *matrix) {
    return (matrix->cols != nullptr) &&
           (matrix->rows != nullptr) &&
           (matrix->weights != nullptr) &&
           (matrix->nrows > -1) &&
           (matrix->nnz > -1);
}

int check_matrix_init(matrix_pcsr_t *matrix) {
    return (matrix->cols != nullptr) &&
           (matrix->row_offsets != nullptr) &&
           (matrix->nrows > -1);
}

int check_matrix_init(matrix_rcsr_t *matrix) {
    return (matrix->cols != nullptr) &&
           (matrix->row_offsets != nullptr) &&
           (matrix->weights != nullptr) &&
           (matrix->nrows > -1);
}

int transpose(matrix_pcsr_t *A, matrix_pcsr_t *B) {

    if(!check_matrix_init(A)) {
        fprintf(stderr, "Input matrix not initialized\n");
        return EXIT_FAILURE;
    }

    int *row_offsets, *cols;
    int nrows = A->nrows;
    int ncols = A->ncols;
    int nnz = A->row_offsets[A->nrows];

    row_offsets = (int*) calloc((ncols + 1), sizeof(*row_offsets));
    assert(row_offsets);
    cols = (int*) malloc(nnz * sizeof(*cols));
    assert(cols);

    /*
     * Compute number of non-zero entries per column of A.
     */
    for (int n = 0; n < nnz; n++){
        row_offsets[A->cols[n]]++;
    }

    /*
     * Compute row offsets (column offsets of A)
     */
    for(int i = 0, psum = 0; i < ncols; i++){
        int temp  = row_offsets[i];
        row_offsets[i] = psum;
        psum += temp;
    }
    row_offsets[ncols] = nnz;

    print_array(row_offsets, ncols - 1);

    /*
     * Copy cols array of A in cols of B
     */
    for(int n = 0; n < nrows; n++){
        for(int j_a = A->row_offsets[n]; j_a < A->row_offsets[n + 1]; j_a++){
            int col  = A->cols[j_a];
            int dest = row_offsets[col];
            cols[dest] = n;
            row_offsets[col]++;
        }
    }

    for(int i = 0, last = 0; i <= ncols; i++){
        int temp  = row_offsets[i];
        row_offsets[i] = last;
        last    = temp;
    }

    B->nrows = ncols;
    B->ncols = nrows;
    B->row_offsets = row_offsets;
    B->cols = cols;

    return EXIT_SUCCESS;
}

int coo_to_csr(matrix_pcoo_t *A, matrix_pcsr_t *B) {

    if(!check_matrix_init(A)) {
        fprintf(stderr, "Input matrix not initialized\n");
        return EXIT_FAILURE;
    }

    int *row_offsets, *cols;
    int *rows = A->rows;    // row indices of A
    int nnz = A->nnz;       // number of nnz in A
    int nrows = A->nrows;   // number of rows in A
    int ncols = A->ncols;

    row_offsets = (int *) calloc((nrows + 1), sizeof(*row_offsets));
    assert(row_offsets);
    cols = (int *) malloc(nnz * sizeof(*cols));
    assert(cols);

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

void print_matrix(matrix_pcoo_t *matrix) {

    if (!check_matrix_init(matrix)) {
        fprintf(stderr, "The matrix is not initialized\n");
        return;
    }

    printf("nrows = %d\n", matrix->nrows);
    printf("ncols = %d\n", matrix->ncols);
    printf("nnz = %d\n", matrix->nnz);
    printf("rows = \n");
    print_array(matrix->rows, matrix->nnz - 1);
    printf("cols = \n");
    print_array(matrix->cols, matrix->nnz - 1);
}

void print_matrix(matrix_pcsr_t *matrix) {

    if (!check_matrix_init(matrix)) {
        fprintf(stderr, "The matrix is not initialized\n");
        return;
    }

    int nnz = matrix->row_offsets[matrix->nrows];
    printf("nrows = %d\n", matrix->nrows);
    printf("ncols = %d\n", matrix->ncols);
    printf("offsets = \n");
    print_array(matrix->row_offsets, matrix->nrows);
    printf("cols = \n");
    print_array(matrix->cols, nnz - 1);
}

void free_matrix(matrix_pcoo_t *matrix) {
    free(matrix->rows);
    free(matrix->cols);
    matrix->rows = nullptr;
    matrix->cols = nullptr;
    matrix->nnz = -1;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void free_matrix(matrix_pcsr_t *matrix) {
    free(matrix->row_offsets);
    free(matrix->cols);
    matrix->row_offsets = nullptr;
    matrix->cols = nullptr;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void free_matrix(matrix_rcoo_t *matrix) {
    free(matrix->rows);
    free(matrix->cols);
    free(matrix->weights);
    matrix->rows = nullptr;
    matrix->cols = nullptr;
    matrix->nnz = -1;
    matrix->nrows = -1;
    matrix->ncols = -1;
}

void free_matrix(matrix_rcsr_t *matrix) {
    free(matrix->row_offsets);
    free(matrix->cols);
    free(matrix->weights);
    matrix->row_offsets = nullptr;
    matrix->cols = nullptr;
    matrix->nrows = -1;
    matrix->ncols = -1;
}
