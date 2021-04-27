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

/*
 * from https://stackoverflow.com/questions/3437404/min-and-max-in-c
 */
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#ifdef __cplusplus
extern "C" {
    #include "mmio.h"
}
#include <cstdbool>
#include <cstring>
#include <cstdlib>
#include <climits>
#else
#include "mmio.h"
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#endif

#include "matstorage.h"

/*
 * Thanks to: https://stackoverflow.com/questions/4849986/how-can-i-check-the-file-extensions-in-c
 */
static int has_extension(const char* name, const char* extension, size_t length)
{
    const char* ldot = strrchr(name, '.');
    if (ldot != NULL)
    {
        if (length == 0)
            length = strlen(extension);
        return strncmp(ldot + 1, extension, length) == 0;
    }
    return 0;
}

//static int cmpfunc (const void * a, const void * b) {
//    return ( *(int*)a - *(int*)b );
//}
//        qsort(cols, nnz, sizeof(int), cmpfunc);

/*
 * Read one line of the file, return TRUE if successful, FALSE if EOF.
 * Thanks to: https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CHOLMOD/Check/cholmod_read.c
 */
static int get_line (FILE *f, char *buf) {
    buf [0] = '\0' ;
    buf [1] = '\0' ;
    buf [BUFFER_SIZE] = '\0' ;
    return (fgets (buf, BUFFER_SIZE, f) != NULL) ;
}

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
int read_matrix_market(const char *fname, matrix_coo_t *m_coo) {

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

    if(m == 0 || n == 0 || nnz == 0) {
        fprintf(stderr, "An empty matrix was given\n");
        return EXIT_FAILURE;
    }

    if(m != n) {
        fprintf(stderr, "An adjacency matrix must be square\n");
        return EXIT_FAILURE;
    }

    /*
     * Allocate memory for the matrix.
     */
    size_t size = mm_is_symmetric(matcode) ? 2 * nnz + 1 : nnz + 1;

    rows = (int *) malloc(size * sizeof(*rows));
    cols = (int *) malloc(size * sizeof(*cols));
    assert(rows);
    assert(cols);

    int i = 0;

    for(int cnt = 0; cnt < nnz; cnt++) {
        int tmp_col, tmp_row;

        /*
         * Premature end of file.
         */
        if (!get_line (f, buf))
           return EXIT_FAILURE;

        nitems = sscanf (buf, "%d %d\n", &tmp_row, &tmp_col) ;

        /*
         * Check number of entries.
         */
        if (nitems < 2 || nitems > 3
                || tmp_row > INT_MAX || tmp_col > INT_MAX) {
            return EXIT_FAILURE;
        }

        /*
         * Check indexing, only for the first pair/triplet.
         */
        if(cnt == 0) {
            if(tmp_row == 0 || tmp_col == 0)
                one_based = false;
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

int read_matrix(const char *fname, matrix_coo_t *m_coo)
{
    int status;

    if (has_extension(fname, "mtx", strlen(fname)) ||
        has_extension(fname, "mm", strlen(fname)))

        status = read_matrix_market(fname, m_coo);
    else {
        fprintf(stderr, "Unsupported file type\n");
        status = EXIT_FAILURE;
    }

    return status;
}

#endif //MATIO_H
