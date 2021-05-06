/****************************************************************************
 *
 * matstorage.h - Functions for converting between different internal storage
 * matrix representations
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
 * --------------------------------------------------------------------------
 *
 * This header file provides functions to convert from the COOrdinate format
 * to the Compressed Sparse Rows. In addition it provides a way of representing
 * adjacency matrices with structures of arrays.
 *
 * The conversion function is based on the one found in the scipy source
 * code, specifically from
 * https://github.com/scipy/scipy/blob/f2ef65dc7f00672496d7de6154744fee55ef95e9/scipy/sparse/sparsetools/coo.h#L33
 *
 * TODO: Remove duplicated edges
 *
 ****************************************************************************/

#pragma once
#ifndef MATSTORAGE_H
#define MATSTORAGE_H

#include <cstdlib>

typedef struct matrix_pcoo_t {
    int nrows;// = ncols since adj matrix is a square matrix
    int nnz;
    int *rows;// row index for each non-zero value
    int *cols;// column index for each non-zero value
} matrix_pcoo_t;

typedef struct matrix_rcoo_t {
    int nrows;// = ncols since adj matrix is a square matrix
    int nnz;
    int *rows;   // row index for each non-zero value
    int *cols;   // column index for each non-zero value
    int *weights;// value of each entry
} matrix_rcoo_t;

typedef struct matrix_pcsr_t {
    int nrows;
    //    int nnz;  // found at row_offsets[nrows]
    int *row_offsets;// offset in columns
    int *cols;       // column index for each non-zero value
} matrix_pcsr_t;

typedef struct matrix_rcsr_t {
    int nrows;
    int *row_offsets;// offset in columns
    int *cols;       // column index for each non-zero value
    int *weights;    // value of each entry
} matrix_rcsr_t;

void check_bc(matrix_pcsr_t g, const float *bc_cpu, const float *bc_gpu);

int check_matrix_pcoo_init(matrix_pcoo_t *matrix);

/**
 * Convert a matrix A, stored in COO format, to a matrix B, stored in the CSR
 * format.
 *
 * At the end (row_offsets, cols) forms a CSR representation,
 * with possible duplicates if they were present in the COO format.
 * This means that duplicate entries in the COO matrix are not merged.
 *
 * Row and column indices *are not* assumed to be ordered.
 *
 * The conversion algorithm has linear complexity.
 * Specifically O(nnz(A) + max(nrows, ncols)).
 *
 * As long as the average number of non-zeroes per row is > 1,
 * this format saves space relative to COO.
 *
 * @param nrows_dense
 * @param m_coo
 * @param m_csr structure representing the matrix in the new format.
 * Row_offsets and cols fields *must not* be preallocated.
 */
void pcoo_to_pcsr(matrix_pcoo_t *m_coo, matrix_pcsr_t *m_csr);

void rcoo_to_rcsr(matrix_pcoo_t *m_coo, matrix_pcsr_t *m_csr);

void print_matrix_coo(matrix_pcoo_t *matrix);

void print_matrix_csr(matrix_pcsr_t *matrix);

void free_matrix_pcoo(matrix_pcoo_t *matrix);

void free_matrix_pcsr(matrix_pcsr_t *matrix);

void free_matrix_rcoo(matrix_rcoo_t *matrix);

void free_matrix_rcsr(matrix_rcsr_t *matrix);

#endif// MATSTORAGE_H
