/****************************************************************************
 * @file matio.cpp
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
 ****************************************************************************/

#include "matio.h"

/**
 * @brief Check if a filename has a given extension.
 *
 * @see https://stackoverflow.com/questions/4849986/how-can-i-check-the-file-extensions-in-c
 */
static int has_extension(const char *name, const char *extension,
                         size_t length) {
    const char *ldot = strrchr(name, '.');
    if (ldot != 0) {
        if (length == 0)
            length = strlen(extension);
        return strncmp(ldot + 1, extension, length) == 0;
    }
    return 0;
}

/*
 * Read one line of the file, return TRUE if successful, FALSE if EOF.
 */
static int get_line(FILE *f, char *buf) {
    buf[0] = '\0';
    buf[1] = '\0';
    buf[BUFFER_SIZE - 1] = '\0';
    return (fgets(buf, BUFFER_SIZE, f) != 0);
}

int query_gprops(const char *fname, gprops_t *gp) {

    FILE *f;
    MM_typecode matcode;
    int nnz, m, n, nitems = 0;
    char buf[BUFFER_SIZE];

    f = fopen(fname, "r");

    if (f == 0) {
        ZF_LOGF("Could not open %s", fname);
        return EXIT_FAILURE;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        ZF_LOGF("Could not process Matrix Market banner");
        return EXIT_FAILURE;
    }

    /*
     * Only matrices of reals or pattern are supported.
     */
    if (!((mm_is_pattern(matcode) || mm_is_real(matcode)) &&
          mm_is_matrix(matcode))) {
        ZF_LOGF("This application does not support Market Matrix type: %s",
                mm_typecode_to_str(matcode));
        return EXIT_FAILURE;
    }

    /*
     * Undirected graphs must be stored as symmetric matrices.
     */
    if (mm_is_symmetric(matcode)) {
        gp->is_directed = false;
    } else {
        gp->is_directed = true;
    }

    /*
     * Get the shape of the sparse matrix and the number of non zero elements.
     */
    if (mm_read_mtx_crd_size(f, &m, &n, &nnz) != 0) {
        ZF_LOGF("Could not read shape and nnz elements of the matrix");
        return EXIT_FAILURE;
    }

    if (m == 0 || n == 0 || nnz == 0) {
        ZF_LOGF("An empty matrix was given");
        return EXIT_FAILURE;
    }

    int j = 0;

    get_line(f, buf);
    int tmp_col, tmp_row;
    while (buf[j] != '\0') {
        if (buf[j] == ' ' || buf[j] == '\n' || buf[j] == '\t') {
            nitems++;
        }
        j++;
    }

    if (nitems == 2 && mm_is_pattern(matcode)) {
        int tmp_wgh;
        sscanf(buf, "%d %d %d\n", &tmp_row, &tmp_col, &tmp_wgh);
        gp->is_weighted = false;
    } else if (nitems == 3 && mm_is_real(matcode)) {
        sscanf(buf, "%d %d\n", &tmp_row, &tmp_col);
        gp->is_weighted = true;
    } else {
        int sup_entr = 0;
        if (mm_is_pattern(matcode))
            sup_entr = 2;
        else if (mm_is_real(matcode))
            sup_entr = 3;

        ZF_LOGF("%d entries given but only %d entries per row are supported",
                nitems, sup_entr);
        return EXIT_FAILURE;
    }

    gp->is_connected = -1;
    close_stream(f);

    return 0;
}

int read_header(FILE *f, MM_typecode *matcode, int *m, int *n, int *nnz) {

    if (mm_read_banner(f, matcode) != 0) {
        ZF_LOGF("Could not process Matrix Market banner");
        return EXIT_FAILURE;
    }

    /*
     * Get the shape of the sparse matrix and the number of non zero elements.
     */
    if (mm_read_mtx_crd_size(f, m, n, nnz) != 0) {
        ZF_LOGF("Could not read shape and nnz elements of the matrix");
        return EXIT_FAILURE;
    }

    if (*m == 0 || *n == 0 || *nnz == 0) {
        ZF_LOGF("An empty matrix was given");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int read_mm(FILE *f, int *nnz, const int *m, const int *n,
            int *rows, int *cols, int *weights, gprops_t *gp) {

    int rmax = 0, cmax = 0, i = 0, nitems;
    bool one_based = true;
    char buf[BUFFER_SIZE];

    if (weights == 0) {
        for (int cnt = 0; cnt < *nnz; cnt++) {
            int tmp_row, tmp_col;

            /*
             * Premature end of file.
             */
            if (!get_line(f, buf))
                return EXIT_FAILURE;

            nitems = sscanf(buf, "%d %d\n", &tmp_row, &tmp_col);

            if (cnt == 0 && (tmp_row == 0 || tmp_col == 0))
                one_based = false;

            /*
             * Check number of entries and their value.
             */
            if (nitems != 2) {
                ZF_LOGF("%d entries given but only 2 entries per row are "
                        "supported on pattern matrices\n",
                        nitems);
                return EXIT_FAILURE;
            }

            if (tmp_row < 0 || tmp_col < 0 ||
                tmp_row > INT_MAX || tmp_col > INT_MAX) {
                ZF_LOGF("Indices out of range");
                return EXIT_FAILURE;
            }

            if (gp->has_self_loops || tmp_col != tmp_row) {
                if (!gp->is_directed && tmp_col != tmp_row) {
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
    } else {
        for (int cnt = 0; cnt < *nnz; cnt++) {
            int tmp_row, tmp_col, tmp_wgh;

            /*
             * Premature end of file.
             */
            if (!get_line(f, buf))
                return EXIT_FAILURE;

            nitems = sscanf(buf, "%d %d %d\n", &tmp_row, &tmp_col, &tmp_wgh);

            if (cnt == 0 && (tmp_row == 0 || tmp_col == 0))
                one_based = false;

            /*
             * Check number of entries and their value.
             */
            if (nitems != 3) {
                ZF_LOGF("%d entries given but only 3 entries per row are "
                        "supported on real matrices\n",
                        nitems);
                return EXIT_FAILURE;
            }

            if (tmp_row < 0 || tmp_col < 0 || tmp_wgh < 0 ||
                tmp_row > INT_MAX || tmp_col > INT_MAX || tmp_wgh > INT_MAX) {
                ZF_LOGF("Indices out of range");
                return EXIT_FAILURE;
            }

            if (gp->has_self_loops && tmp_col != tmp_row) {
                if (!gp->is_directed) {
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
    }

    if (one_based ? (rmax > *m || cmax > *n) : (rmax >= *m || cmax >= *n)) {
        ZF_LOGF("Indices out of range");
        return EXIT_FAILURE;
    }

    /*
     * If there are self-edges removed and the allocated space is not entirely
     * used reallocate memory.
     */
    if ((!gp->is_directed && (2 * *nnz != i)) ||
        (gp->is_directed && (*nnz != i))) {

        *nnz = i;
        rows = (int *) realloc(rows, (*nnz) * sizeof(int));
        cols = (int *) realloc(cols, (*nnz) * sizeof(int));

        if (weights != 0) {
            weights = (int *) realloc(weights, (*nnz) * sizeof(int));
        }
    } else {
        *nnz = i;
    }

    /*
     * Convert to zero-based representation.
     */
    if (one_based) {
        for (i = 0; i < *nnz; i++) {
            cols[i]--;
            rows[i]--;
        }
    }

    return EXIT_SUCCESS;
}

int read_mm_real(FILE *f, matrix_rcoo_t *m_coo, gprops_t *gp) {

    MM_typecode matcode;
    int nnz, m, n;
    int *rows, *cols, *weights;

    if (read_header(f, &matcode, &m, &n, &nnz)) {
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

    if (read_mm(f, &nnz, &m, &n, rows, cols, weights, gp)) {
        return EXIT_FAILURE;
    }

    m_coo->nnz = nnz;
    m_coo->nrows = n;
    m_coo->rows = rows;
    m_coo->cols = cols;
    m_coo->weights = weights;

    return 0;
}

int read_mm_pattern(FILE *f, matrix_pcoo_t *m_coo, gprops_t *gp) {

    MM_typecode matcode;
    int nnz, m, n;
    int *rows, *cols, *weights = 0;

    if (read_header(f, &matcode, &m, &n, &nnz)) {
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

    if (read_mm(f, &nnz, &m, &n, rows, cols, weights, gp)) {
        return EXIT_FAILURE;
    }

    m_coo->nnz = nnz;
    m_coo->nrows = m;
    m_coo->ncols = n;
    m_coo->rows = rows;
    m_coo->cols = cols;

    return 0;
}

int read_matrix(const char *fname, matrix_pcoo_t *m_coo, gprops_t *gp) {

    if (has_extension(fname, "mtx", strlen(fname)) ||
        has_extension(fname, "mm", strlen(fname))) {

        if (gp->is_weighted) {
            ZF_LOGF("Weighted graph is not supported");
            return EXIT_FAILURE;
        } else {
            FILE *f;
            f = fopen(fname, "r");

            if (f == 0) {
                ZF_LOGF("Could not open %s", fname);
                return EXIT_FAILURE;
            }

            if (read_mm_pattern(f, m_coo, gp)) {
                ZF_LOGF("Error reading matrix");
                return EXIT_FAILURE;
            }

            close_stream(f);
        }

    } else {
        ZF_LOGF("Unsupported file type");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int write_mm_pattern(FILE *f, matrix_pcoo_t *m_coo, bool directed) {

    MM_typecode matcode;

    /*
     * Write the banner.
     */
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_pattern(&matcode);

    if (directed)
        mm_set_general(&matcode);
    else
        mm_set_symmetric(&matcode);

    if (mm_write_banner(f, matcode))
        return EXIT_FAILURE;

    /*
     * Write the header and the values.
     */
    if (mm_write_mtx_crd_size(f, m_coo->nrows, m_coo->nrows, m_coo->nnz))
        return EXIT_FAILURE;

    for (int i = 0; i < m_coo->nnz; i++)
        if (fprintf(f, "%d %d\n", m_coo->rows[i] + 1, m_coo->cols[i] + 1) < 0) {
            return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;
}

int write_mm_real(FILE *f, matrix_rcoo_t *m_coo, bool directed) {

    MM_typecode matcode;

    /*
     * Write the banner.
     */
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    if (directed)
        mm_set_general(&matcode);
    else
        mm_set_symmetric(&matcode);

    if (mm_write_banner(f, matcode))
        return EXIT_FAILURE;

    /*
     * Write the header and the values.
     */
    if (mm_write_mtx_crd_size(f, m_coo->nrows, m_coo->nrows, m_coo->nnz))
        return EXIT_FAILURE;

    for (int i = 0; i < m_coo->nnz; i++)
        if (fprintf(f, "%d %d %d\n", m_coo->rows[i] + 1, m_coo->cols[i] + 1,
                    m_coo->weights[i]) < 0) {
            return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;
}
