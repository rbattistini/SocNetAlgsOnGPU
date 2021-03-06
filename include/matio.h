/****************************************************************************
 * @file matio.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
* @brief Functions for reading and writing Matrix Market files.
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
 * ---------------------------------------------------------------------------
 *
 *  * Here is an example of the matrix market format:
 *
 * +----------------------------------------------+
 * |%%MatrixMarket matrix coordinate real general | <--- header line
 * |%                                             | <--+
 * |% comments                                    |    |-- 0 or more comments
 * |%                                             | <--+
 * |  M N L                                       | <--- rows, columns, entries
 * |  I1 J1 A(I1, J1)                             | <--+
 * |  I2 J2 A(I2, J2)                             |    |
 * |  I3 J3 A(I3, J3)                             |    |-- L lines
 * |     . . .                                    |    |
 * |  IL JL A(IL, JL)                             | <--+
 * +----------------------------------------------+
 ****************************************************************************/

#pragma once
#ifndef MATIO_H
#define MATIO_H

/*
 * The maximum line length of the Matrix Market format.
 */
#define BUFFER_SIZE 1030

#include "graphs.h"
#include "matds.h"
#include "mmio.h"
#include <climits>
#include <cstring>

int query_gprops(const char *fname, gprops_t *gp);

/**
 * @brief Reads a graph stored in a Matrix Market file and stores it in a
 * COOrdinate format matrix.
 *
 * Lacks duplicated entries elimination.
 * Undirected graphs must be stored as unsymmetric matrices to be recognized.
 * Converts one-based indexing to zero-based indexing.
 *
 * @param[in] fname
 * @param[out] m_coo
 * @return 0 if successful, 1 otherwise
 */
int read_mm_pattern(FILE *f, matrix_pcoo_t *m_coo, gprops_t *gp);

int read_mm_real(FILE *f, matrix_rcoo_t *m_coo, gprops_t *gp);

int read_matrix(const char *fname, matrix_pcoo_t *m_coo, gprops_t *gp);

int write_mm_pattern(FILE *f, matrix_pcoo_t *m_coo, bool directed);

int write_mm_real(FILE *f, matrix_rcoo_t *m_coo, bool directed);

#endif//MATIO_H
