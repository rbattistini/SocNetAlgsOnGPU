/****************************************************************************
 *
 * matio.h - Functions for reading and writing Matrix Market files (.mm, .mtx)
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
 ****************************************************************************/

#ifndef MATIO_H
#define MATIO_H

extern "C" {
    #include "mmio.h"
};

#include "matstorage.h"

int readMatrixMarketFile(const char *fname, matrix_coo_t *m_coo,
                         const int nedges) {
    FILE *f;
    MM_typecode matcode;
    int nz, m, n;
    int *rows, *columns;

    f = fopen(fname, "r");

    if( f == nullptr) {
        fprintf(stderr, "Could not open %s\n", fname);
        return EXIT_FAILURE;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner\n");
        return EXIT_FAILURE;
    }

    /*
     * Only matrices of reals are supported.
     */
    if (!mm_is_real(matcode) || !mm_is_matrix(matcode)) {
        fprintf(stderr, "This application does not support\n"
                        "Market Matrix type: %s\n",
                mm_typecode_to_str(matcode));
        return EXIT_FAILURE;
    }

    /*
     * Get the shape of the sparse matrix and the number of non zero elements.
     */
    if (mm_read_mtx_crd_size(f, &m, &n, &nz) != 0) {
        fprintf(stderr, "Could not read shape and nnz elements of the sparse"
                        "matrix\n");
        return EXIT_FAILURE;
    }

    /*
     * Allocate memory for the matrix.
     */
    size_t diag_elems = ( (2 * nz) - nedges ) / 2;
    size_t size = (nz - diag_elems) * 2 + diag_elems;
    rows = (int *) malloc(size * sizeof(*rows));
    columns = (int *) malloc(size * sizeof(*columns));

    if( mm_is_symmetric(matcode) ) {
        int cnt = 0, i = 0;
        while(cnt < nz) {
            int tmp_col, tmp_row;

            int err_code = fscanf(f, "%d %d\n", &tmp_col, &tmp_row);
            if(err_code == 0) {
                fprintf(stderr, "Could not read entry %d", i);
                return EXIT_FAILURE;
            }

            /*
             * Reindex from 0 to 1.
             */
            if(tmp_col != tmp_row) {
                columns[i] = tmp_col - 1;
                columns[i + 1] = tmp_row - 1;
                rows[i] = tmp_row - 1;
                rows[i + 1] = tmp_col - 1;
                i += 2;
            } else {
                columns[i] = tmp_col - 1;
                rows[i] = tmp_row - 1;
                i++;
            }

            cnt++;
        };
    } else {

        /*
         * Reindex from 0 to 1.
         */
        for (int i = 0; i < nz; i++) {
            int err_code = fscanf(f, "%d %d\n", &columns[i], &rows[i]);
            rows[i]--;
            columns[i]--;

            if(err_code == 0) {
                fprintf(stderr, "Could not read entry %d", i);
                return EXIT_FAILURE;
            }
        }
    }

    fclose(f);

    m_coo->nnz = nz;
    m_coo->nrows = m; // m = n
    m_coo->rows = rows;
    m_coo->cols = columns;

#ifdef DEBUG
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, m, n, nz);
    for (int i=0; i<nz; i++)
        fprintf(stdout, "%d %d\n", rows[i] + 1, columns[i] + 1);
#endif

    return 0;
}

// TODO
// int writeMatrixMarketFile()

#endif //MATIO_H
