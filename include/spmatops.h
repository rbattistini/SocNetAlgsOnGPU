/****************************************************************************
 * @file spmatops.h
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
 *
 ****************************************************************************/

#pragma once
#ifndef SPMATOPS_H
#define SPMATOPS_H

#include "common.h"
#include "matds.h"

/**
 * @brief Implements the Gustavson’s row-wise sparse general matrix-matrix
 * multiplication algorithm.
 *
 * @cite gustavson_two_1978
 *
 * Executes C = A * B.
 *
 * @note The operation is performed only on patterns. It works with sorted or
 * unsorted matrices and symmetric or unsymmetric ones.
 *
 * @param A p x q sparse pattern matrix in CSR format
 * @param B q x r sparse pattern matrix in CSR format
 * @param C p x r sparse pattern matrix in CSR format
 * @return 0 if successful, 1 otherwise
 */
int spgemm(matrix_pcsr_t *A,
           matrix_pcsr_t *B,
           matrix_pcsr_t *C);

/**
 * @brief Implements the SpRef function with SpGEMM as the main subroutine as
 * described by Buluç and Gilbert.
 *
 * @cite buluc_parallel_2012
 *
 * Executes R * A * Q = C.
 *
 * @note The operation is performed only on patterns. It works with sorted or
 * unsorted matrices and symmetric or unsymmetric ones.
 *
 * @param R sparse pattern matrix in CSR format
 * @param A sparse pattern matrix in CSR format
 * @param Q sparse pattern matrix in CSR format
 * @param C sparse pattern matrix in CSR format
 * @return 0 if successful, 1 otherwise
 */
int spref(matrix_pcsr_t *R,
          matrix_pcsr_t *A,
          matrix_pcsr_t *Q,
          matrix_pcsr_t *C);

int get_R_matrix(matrix_pcsr_t *R,
                 const int *vertices,
                 int nvertices,
                 int nrows);

#endif//SPMATOPS_H
