/****************************************************************************
 * @file degree.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions to compute the degree of the vertices of an undirected and
 * unweighted graph stored as a sparse pattern matrix in CSR format.
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

#include "degree.h"

void compute_degrees_undirected(matrix_pcsr_t *g, int *degree) {

    if (!check_matrix_pcsr(g)) {
        ZF_LOGF("The graph is not initialized");
        return;
    }

    int *rows = g->row_offsets;        // row indices of A
    int nrows = g->nrows;              // number of rows in A

    fill(degree, nrows, 0);

    /*
     * Compute number of non-zero entries per row of A.
     */
    for (int i = 0; i < nrows; i++) {
        degree[i] = rows[i + 1] - rows[i];
    }
}

void compute_degrees_undirected(matrix_pcoo_t *g, int *degree) {

    if (!check_matrix_pcoo(g)) {
        ZF_LOGF("The graph is not initialized");
        return;
    }

    int *rows = g->rows;   // row indices of A
    int nnz = g->nnz;      // number of nnz in A
    int nrows = g->nrows;  // number of rows in A

    fill(degree, nrows, 0);

    /*
     * Compute number of non-zero entries per row of A.
     */
    for (int n = 0; n < nnz; n++) {
        degree[rows[n]]++;
    }
}

void compute_degrees_directed(matrix_pcoo_t *g, int *in_degree,
                              int *out_degree) {

    if (g->rows == 0) {
        ZF_LOGF("The graph is not initialized");
        return;
    }

    int *rows = g->rows;  // row indices of A
    int nnz = g->nnz;     // number of nnz in A
    int length = g->nrows;// number of rows and columns in A
    int *cols = g->cols;  // column indices of A

    //    offsets = (int *) malloc((length + 1) * sizeof(*offsets));
    fill(in_degree, length, 0);
    fill(out_degree, length, 0);

    /*
     * Compute number of non-zero entries per column of A.
     */
    for (int n = 0; n < nnz; n++) {
        out_degree[rows[n]]++;
    }

    /*
     * Compute number of non-zero entries per row of A.
     */
    for (int n = 0; n < nnz; n++) {
        in_degree[cols[n]]++;
    }
}
