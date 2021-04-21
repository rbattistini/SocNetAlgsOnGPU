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
 * TODO read from SNAP file
 *
 * TODO write to Matrix Market file
 *
 * TODO handle unconnected graphs by extraction of largest SCC
 *
 * REVIEW symmetric matrices parsing in read_matrix_market()
 *
 ****************************************************************************/

#ifndef MATIO_H
#define MATIO_H

extern "C" {
    #include "mmio.h"
};

#include <string>
#include "matstorage.h"

/*
 * Quick and dirty check
 */
std::string get_extension(const std::string& fn) {
    return fn.substr(fn.find_last_of('.') + 1);
}

int read_matrix_market(const char *fname, matrix_coo_t *m_coo) {

    FILE *f;
    MM_typecode matcode;
    int nnz, m, n;
    int *rows, *cols;

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
     * Only matrices of reals or pattern are supported.
     */
    if (!(( mm_is_pattern(matcode) || mm_is_real(matcode) )
        && mm_is_matrix(matcode))) {
        fprintf(stderr, "This application does not support\n"
                        "Market Matrix type: %s\n",
                mm_typecode_to_str(matcode));
        return EXIT_FAILURE;
    }

    /*
     * Get the shape of the sparse matrix and the number of non zero elements.
     */
    if (mm_read_mtx_crd_size(f, &m, &n, &nnz) != 0) {
        fprintf(stderr, "Could not read shape and nnz elements of the sparse"
                        "matrix\n");
        return EXIT_FAILURE;
    }

    /*
     * Allocate memory for the matrix.
     */
//    size_t diag_elems = ((2 * nnz) - nedges ) / 2;
//    size_t size = (nnz - diag_elems) * 2 + diag_elems;
    size_t size = 2 * nnz + 1;
    rows = (int *) malloc(size * sizeof(*rows));
    cols = (int *) malloc(size * sizeof(*cols));

    if( mm_is_symmetric(matcode) ) {
        int cnt = 0, i = 0;
        while(cnt < nnz) {
            int tmp_col, tmp_row;

            int err_code = fscanf(f, "%d %d\n", &tmp_row, &tmp_col);
            if(err_code == 0) {
                fprintf(stderr, "Could not read entry %d", i);
                return EXIT_FAILURE;
            }

            /*
             * Reindex from 0 to 1.
             */
            if(tmp_col != tmp_row) { // if not diagonal
                cols[i] = tmp_col;
                cols[i + 1] = tmp_row;
                rows[i] = tmp_row;
                rows[i + 1] = tmp_col;
                i += 2;
            } else {
                cols[i] = tmp_col;
                rows[i] = tmp_row;
                i++;
            }

            cnt++;
        };

        /*
         * Resize if there is excess space, i. e. if there is at least one
         * value on the matrix diagonal.
         */
        if( 2 * nnz != i ){
            rows = (int *) realloc(rows, i * sizeof(*rows));
            cols = (int *) realloc(cols, i * sizeof(*cols));
        }

        nnz = i;

    } else {

        for (int i = 0; i < nnz; i++) {
            int err_code = fscanf(f, "%d %d\n", &rows[i], &cols[i]);

            if(err_code == 0) {
                fprintf(stderr, "Could not read entry %d", i);
                return EXIT_FAILURE;
            }
        }
    }

    fclose(f);

    m_coo->nnz = nnz;
    m_coo->nrows = m; // m = n
    m_coo->rows = rows;
    m_coo->cols = cols;

    return 0;
}

int read_snap(const char *fname, matrix_coo_t *m_coo) {
    return EXIT_FAILURE;
};

int read_matrix(const char *fname, matrix_coo_t *m_coo)
{
    std::string ext = get_extension(fname);
    int status = 0;

    if (ext == "mtx" || ext == "mm")
        status = read_matrix_market(fname, m_coo);
    else if(ext == "txt")
        status = read_snap(fname, m_coo);
    else {
        fprintf(stderr, "Unsupported file type\n");
        return EXIT_FAILURE;
    }

    return status;
}

#endif //MATIO_H
