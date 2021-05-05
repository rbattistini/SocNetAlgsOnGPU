/****************************************************************************
 *
 * matio.cpp - Functions for reading and writing external matrix storage file
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
 ****************************************************************************/

#include "matio.h"
#include <cassert>

inline int max(int a, int b) {
    return a > b ? a : b;
}

/*
 * Check if a filename has a given extension.
 * Thanks to: https://stackoverflow.com/questions/4849986/how-can-i-check-the-file-extensions-in-c
 */
static int has_extension(const char* name, const char* extension, size_t length)
{
    const char* ldot = strrchr(name, '.');
    if (ldot != nullptr)
    {
        if (length == 0)
            length = strlen(extension);
        return strncmp(ldot + 1, extension, length) == 0;
    }
    return 0;
}

/*
 * Read one line of the file, return TRUE if successful, FALSE if EOF.
 * Thanks to: https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CHOLMOD/Check/cholmod_read.c
 */
static int get_line (FILE *f, char *buf) {
    buf [0] = '\0' ;
    buf [1] = '\0' ;
    buf [BUFFER_SIZE] = '\0' ;
    return (fgets (buf, BUFFER_SIZE, f) != nullptr) ;
}

int query_gprops(const char *fname, gprops_t *gp) {

    FILE *f;
    MM_typecode matcode;
    int nnz, m, n, nitems = 0;
    char buf[BUFFER_SIZE];

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
     * Undirected graphs must be stored as symmetric matrices.
     */
    if(mm_is_symmetric(matcode)) {
        gp->is_directed = false;
    }

    /*
     * Get the shape of the sparse matrix and the number of non zero elements.
     */
    if (mm_read_mtx_crd_size(f, &m, &n, &nnz) != 0) {
        fprintf(stderr, "Could not read shape and nnz elements of the sparse"
                        "matrix\n");
        return EXIT_FAILURE;
    }

    if(m == 0 || n == 0 || nnz == 0) {
        fprintf(stderr, "An empty matrix was given\n");
        return EXIT_FAILURE;
    }

    if(m != n) {
        fprintf(stderr, "An adjacency matrix must be square\n");
        return EXIT_FAILURE;
    }

    int j = 0;

    get_line (f, buf);
    int tmp_col, tmp_row;
    while(buf[j] != '\0') {
        if(buf[j] == ' ' || buf[j] == '\n' || buf[j] == '\t') {
            nitems++;
        }
        j++;
    }

    if(nitems == 2 && mm_is_pattern(matcode)) {
        int tmp_wgh;
        sscanf(buf, "%d %d %d\n", &tmp_row, &tmp_col, &tmp_wgh);
        gp->is_weighted = true;
    } else if(nitems == 3 && mm_is_real(matcode)) {
        sscanf(buf, "%d %d\n", &tmp_row, &tmp_col);
        gp->is_weighted = false;
    } else {
        fprintf(stderr, "Only two or three entries per row are supported\n");
        fprintf(stderr, "Pattern must have two entries per row\n");
        fprintf(stderr, "Real must have three entries per row\n");
        return EXIT_FAILURE;
    }

    gp->is_connected = -1;

    return 0;
}

int read_matrix_market_real(const char *fname, matrix_rcoo_t *m_coo) {

    FILE *f;
    MM_typecode matcode;
    int nnz, m, n, rmax = 0, cmax = 0, nitems = 0;
    bool one_based = true;
    int *rows, *cols, *weights;
    char buf[BUFFER_SIZE];

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
    size_t size = mm_is_symmetric(matcode) ? 2 * nnz + 1 : nnz + 1;

    rows = (int *) malloc(size * sizeof(*rows));
    cols = (int *) malloc(size * sizeof(*cols));
    weights = (int *) malloc(size * sizeof(*weights));
    assert(rows);
    assert(cols);
    assert(weights);

    int i = 0;
    for(int cnt = 1; cnt < nnz; cnt++) {
        int tmp_row, tmp_col, tmp_wgh;

        /*
         * Premature end of file.
         */
        if (!get_line (f, buf))
            return EXIT_FAILURE;

        nitems = sscanf(buf, "%d %d %d\n", &tmp_row, &tmp_col, &tmp_wgh);

        if(cnt == 0 && (tmp_row == 0 || tmp_col == 0))
            one_based = false;

        /*
         * Check number of entries and their value.
         */
        if (nitems != 3 || tmp_row > INT_MAX || tmp_col > INT_MAX) {
            return EXIT_FAILURE;
        }

        /*
         * Discard self-loops.
         */
        if(tmp_col != tmp_row) {
            if( mm_is_symmetric(matcode) ) {
                cols[i] = tmp_col;
                rows[i] = tmp_row;
                weights[i] = tmp_wgh;
                rmax = max(rows[i], rmax);
                cmax = max(cols[i], cmax);
                i++;

                cols[i] = tmp_row;
                rows[i] = tmp_col;
                weights[i] = tmp_wgh;
                rmax = max(rows[i], rmax);
                cmax = max(cols[i], cmax);
                i++;
            } else {
                cols[i] = tmp_col;
                rows[i] = tmp_row;
                weights[i] = tmp_wgh;
                rmax = max(rows[i], rmax);
                cmax = max(cols[i], cmax);
                i++;
            }
        }
    }

    if (one_based ? (rmax > n || cmax > n) : (rmax >= n || cmax >= n)) {
        fprintf(stderr, "Indices out of range\n");
        return EXIT_FAILURE;
    }

    /*
     * If there are self-edges removed and the allocated space is not entirely
     * used reallocate memory.
     */
    if( (mm_is_symmetric(matcode)  && (2 * nnz != i) ) ||
        (!mm_is_symmetric(matcode) && (nnz != i) ) ) {
        nnz = i;
        rows = (int*) realloc(rows, (nnz) * sizeof(int));
        cols = (int*) realloc(cols, (nnz) * sizeof(int));
    } else {
        nnz = i;
    }

    /*
     * Convert to zero-based representation.
     */
    if(one_based) {
        for(i = 0; i < nnz; i++) {
            cols[i]--;
            rows[i]--;
        }
    }

    fclose(f);

    m_coo->nnz = nnz;
    m_coo->nrows = n;
    m_coo->rows = rows;
    m_coo->cols = cols;
    m_coo->weights = weights;

    return 0;
}

int read_matrix_market_pattern(const char *fname, matrix_pcoo_t *m_coo) {

    FILE *f;
    MM_typecode matcode;
    int nnz, m, n, rmax = 0, cmax = 0, nitems = 0;
    bool one_based = true;
    int *rows, *cols;
    char buf[BUFFER_SIZE];

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
     * Get the shape of the sparse matrix and the number of non zero elements.
     */
    if (mm_read_mtx_crd_size(f, &m, &n, &nnz) != 0) {
        fprintf(stderr, "Could not read shape and nnz elements of the sparse"
                        "matrix\n");
        return EXIT_FAILURE;
    }

    int i = 0;

    /*
     * Allocate memory for the matrix.
     */
    size_t size = mm_is_symmetric(matcode) ? 2 * nnz + 1 : nnz + 1;

    rows = (int *) malloc(size * sizeof(*rows));
    cols = (int *) malloc(size * sizeof(*cols));
    assert(rows);
    assert(cols);

    for(int cnt = 1; cnt < nnz; cnt++) {
        int tmp_row, tmp_col;
        /*
         * Premature end of file.
         */
        if (!get_line (f, buf))
            return EXIT_FAILURE;

        nitems = sscanf(buf, "%d %d\n", &tmp_row, &tmp_col);

        if(cnt == 0 && (tmp_row == 0 || tmp_col == 0))
            one_based = false;

        /*
         * Check number of entries and their value.
         */
        if (nitems != 2 || tmp_row > INT_MAX || tmp_col > INT_MAX) {
            return EXIT_FAILURE;
        }

        /*
         * Discard self-loops.
         */
        if(tmp_col != tmp_row) {
            if( mm_is_symmetric(matcode) ) {
                cols[i] = tmp_col;
                rows[i] = tmp_row;
                rmax = max(rows[i], rmax);
                cmax = max(cols[i], cmax);
                i++;

                cols[i] = tmp_row;
                rows[i] = tmp_col;
                rmax = max(rows[i], rmax);
                cmax = max(cols[i], cmax);
                i++;
            } else {
                cols[i] = tmp_col;
                rows[i] = tmp_row;
                rmax = max(rows[i], rmax);
                cmax = max(cols[i], cmax);
                i++;
            }
        }
    }

    if (one_based ? (rmax > n || cmax > n) : (rmax >= n || cmax >= n)) {
        fprintf(stderr, "Indices out of range\n");
        return EXIT_FAILURE;
    }

    /*
     * If there are self-edges removed and the allocated space is not entirely
     * used reallocate memory.
     */
    if( (mm_is_symmetric(matcode)  && (2 * nnz != i) ) ||
        (!mm_is_symmetric(matcode) && (nnz != i) ) ) {
        nnz = i;
        rows = (int*) realloc(rows, (nnz) * sizeof(int));
        cols = (int*) realloc(cols, (nnz) * sizeof(int));
    } else {
        nnz = i;
    }

    /*
     * Convert to zero-based representation.
     */
    if(one_based) {
        for(i = 0; i < nnz; i++) {
            cols[i]--;
            rows[i]--;
        }
    }

    fclose(f);

    m_coo->nnz = nnz;
    m_coo->nrows = n;
    m_coo->rows = rows;
    m_coo->cols = cols;

    return 0;
}

int read_matrix(const char *fname, matrix_pcoo_t *m_coo, gprops_t *gp) {

    if (has_extension(fname, "mtx", strlen(fname)) ||
        has_extension(fname, "mm", strlen(fname))) {

        if(gp->is_weighted) {
            fprintf(stderr, "Weighted graph is not supported\n");
            return EXIT_FAILURE;
        } else {
            read_matrix_market_pattern(fname, m_coo);
        }

    } else {
        fprintf(stderr, "Unsupported file type\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
