/****************************************************************************
 *
 * matio.h - Functions for reading and writing external matrix storage file
 * formats (COO or edge list)
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
 * ---------------------------------------------------------------------------
 *
 * Note that only two types of files are supported and each with its own
 * restrictions.
 *
 * The only type of graph supported is the undirected and with uniform weights
 * one. Connectedness and the presence of loops are not checked.
 *
 * TODO handle unconnected graphs by extraction of largest SCC
 *
 * TODO allow directed -> undirected transformation on demand
 *
 * TODO handle duplicated entries
 *
 * TODO check for too big numbers (strtol)
 *
 ****************************************************************************/

#ifndef MATIO_H
#define MATIO_H

/*
 * The maximum line length of the Matrix Market format.
 */
#define BUFFER_SIZE 1030

#include "matstorage.h"
#include "graphs.h"
#include <cstdio>
#include <cstdbool>
#include <cstring>
#include <cstdlib>
#include <climits>
extern "C" {
    #include "mmio.h"
}

int query_gprops(const char *fname, gprops_t *gp);

/**
 * @brief Reads a MARKET graph from an input-stream into a CSR sparse format
 *
 * Here is an example of the matrix market format
 * +----------------------------------------------+
 * |%%MatrixMarket matrix coordinate real general | <--- header line
 * |%                                             | <--+
 * |% comments                                    |    |-- 0 or more comment lines
 * |%                                             | <--+
 * |  M N L                                       | <--- rows, columns, entries
 * |  I1 J1 A(I1, J1)                             | <--+
 * |  I2 J2 A(I2, J2)                             |    |
 * |  I3 J3 A(I3, J3)                             |    |-- L lines
 * |     . . .                                    |    |
 * |  IL JL A(IL, JL)                             | <--+
 * +----------------------------------------------+
 *
 * Read a matrix market file. Removes self-loops.
 * Guaranteed ordered cols. Lacks duplicated entries (edges) elimination.
 * Symmetric graphs are meant to be undirected. Non-symmetric matrices are not
 * converted to undirected graphs. Converts one-based to zero-based. Validates
 * indices' ranges. Accepts only adjacency matrices.
 *
 * @param fname
 * @param m_coo
 * @return
 */
int read_matrix_market_pattern(const char *fname, matrix_pcoo_t *m_coo);

int read_matrix_market_real(const char *fname, matrix_rcoo_t *m_coo);

int read_matrix(const char *fname, matrix_pcoo_t *m_coo, gprops_t *gp);

#endif //MATIO_H
