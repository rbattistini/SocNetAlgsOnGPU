/****************************************************************************
 * @file matstorage.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions to convert from the COOrdinate format to the Compressed
 * Sparse Rows. In addition it provides a way of representing adjacency
 * matrices with structures of arrays.
 *
 * @see https://github.com/scipy/scipy/tree/master/scipy/sparse/sparsetools
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
 * --------------------------------------------------------------------------
 *
 * TODO: Handling duplicated edges by removing them
 *
 ****************************************************************************/

#pragma once
#ifndef MATSTORAGE_H
#define MATSTORAGE_H

#include "common.h"
#include <cmath>

typedef struct matrix_pcoo_t {
    int nrows;
    int ncols;
    int nnz;
    int *rows;// row index for each non-zero value
    int *cols;// column index for each non-zero value
} matrix_pcoo_t;

typedef struct matrix_rcoo_t : matrix_pcoo_t {
    int *weights;// value of each entry
} matrix_rcoo_t;

typedef struct matrix_pcsr_t {
    int nrows;
    int ncols;
    int *row_offsets;// offset in columns
    int *cols;       // column index for each non-zero value
} matrix_pcsr_t;

typedef struct matrix_rcsr_t : matrix_pcsr_t {
    int *weights;// value of each entry
} matrix_rcsr_t;


int check_matrix_pcoo(matrix_pcoo_t *matrix);

int check_matrix_pcsr(matrix_pcsr_t *matrix);

int check_matrix_rcoo(matrix_rcoo_t *matrix);

int check_matrix_rcsr(matrix_rcsr_t *matrix);

void free_matrix_pcoo(matrix_pcoo_t *matrix);

void free_matrix_pcsr(matrix_pcsr_t *matrix);

void free_matrix_rcoo(matrix_rcoo_t *matrix);

void free_matrix_rcsr(matrix_rcsr_t *matrix);

void print_matrix_pcoo(matrix_pcoo_t *matrix);

void print_matrix_pcsr(matrix_pcsr_t *matrix);

/**
 * @brief Expand a compressed row pointer into a row array.
 *
 * @note Output array Bi must be preallocated
 *
 * @param nrows number of rows in A
 * @param row_offsets row pointer
 * @param rows row indices
 */
void expand_row_pointer(int nrows, const int *row_offsets, int *rows);

/**
 * @brief Computes A = B, where A is a pattern matrix in COOrdinate format and
 * B is a pattern matrix in CSR format.
 *
 * @note Duplicate entries in the COO matrix are not merged.
 * @note Row and column indices are not assumed to be ordered.
 * @note Row_offsets and cols fields *must not* be preallocated.
 * @note As long as the average number of non-zeroes per row is > 1, this
 * format saves space relative to COO.
 *
 * @param A sparse pattern matrix in COO format
 * @param B sparse pattern matrix in CSR format
 * @return 0 if successful, 1 otherwise
 */
int coo_to_csr(matrix_pcoo_t *A, matrix_pcsr_t *B);

/**
 * @brief Computes A^T, in which A is an m x n CSR format sparse pattern A.
 *
 * @note Assumes that A is not stored to exploit symmetry.
 *
 * @note Transposition of a CSR A is the equivalent of converting the
 * A to the CSC format.
 *
 * @param A sparse pattern matrix, can be symmetric or unsymmetric
 * @param B  the transpose of A
 * @return 0 if successful, 1 otherwise
 */
int transpose(matrix_pcsr_t *A, matrix_pcsr_t *B);

#endif// MATSTORAGE_H
