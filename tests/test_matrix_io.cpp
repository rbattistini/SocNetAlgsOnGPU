/****************************************************************************
 *
 * test_matrix_io.cpp
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

#include "graphs.h"
#include "matio.h"
#include "matstorage.h"
#include "tests.h"

#define MM_HEADER 51

static matrix_pcoo_t coo;
static matrix_pcsr_t csr;
static gprops_t gprops;

static void write_tmp_file(FILE *tmp, char *tmp_mm_header, int *tmp_header,
                    int *tmp_row, int *tmp_col, int *tmp_wgh) {

    fprintf(tmp, "%s\n", tmp_mm_header);
    fprintf(tmp, "%d %d %d\n", tmp_header[0], tmp_header[1], tmp_header[2]);

    if(tmp_wgh == nullptr) {
        for (int i = 0; i < tmp_header[2]; i++) {
            fprintf(tmp, "%d %d\n", tmp_row[i], tmp_col[i]);
        }
    } else {
        for (int i = 0; i < tmp_header[2]; i++) {
            fprintf(tmp, "%d %d %d\n", tmp_row[i], tmp_col[i], tmp_wgh[i]);
        }
    }
}

static inline unsigned random_letter() {
    long l;
    do { l = random(); } while (l>=(RAND_MAX/26)*26);
    return (unsigned)(l % 26);
}

static inline char gen_random_char() {
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[random_letter () % 26];
}

static inline unsigned gen_random_in_range(unsigned lower, unsigned upper) {
    return (rand() % (upper - lower + 1)) + lower;
}

static void write_random_char_to_file(FILE *f, int nit) {

    fseek(f, 0, SEEK_END);
    unsigned file_length = ftell(f);

    for(int i = 0; i < nit; i++) {
        unsigned n = gen_random_in_range(MM_HEADER + 1, file_length);
        fseek(f, n, SEEK_SET);
        putc(gen_random_char(), f);
    }

    rewind(f);
}

TEST_CASE("Test graph loading with properties querying") {

    SUBCASE("undirected unweighted graph") {
        const char* gname = "und_unw_test.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate pattern symmetric";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 7, 5, 8};
        int tmp_col[] = {0, 0, 0, 0, 1, 2, 2, 4, 7};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, nullptr);

        fclose(tmp);
        tmp = fopen(gname, "r");

        query_gprops(gname, &gprops);

        CHECK_EQ(gprops.is_directed, false);
        CHECK_EQ(gprops.is_weighted, false);

        fclose(tmp);
    }

    SUBCASE("undirected weighted graph") {
        const char* gname = "und_wgh_test.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate real symmetric";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 7, 5, 8};
        int tmp_col[] = {0, 0, 0, 0, 1, 2, 2, 4, 7};
        int tmp_weight[] = {1, 2, 3, 5, 2, 4, 3, 2, 1};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, tmp_weight);

        fclose(tmp);
        tmp = fopen(gname, "r");
        query_gprops(gname, &gprops);

        CHECK_EQ(gprops.is_directed, false);
        CHECK_EQ(gprops.is_weighted, true);

        fclose(tmp);
    }

    SUBCASE("directed unweighted graph") {
        const char* gname = "dir_unw_test.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate pattern general";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 7, 5, 8};
        int tmp_col[] = {0, 0, 0, 0, 1, 2, 2, 4, 7};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, nullptr);

        fclose(tmp);
        tmp = fopen(gname, "r");
        query_gprops(gname, &gprops);

        CHECK_EQ(gprops.is_directed, true);
        CHECK_EQ(gprops.is_weighted, false);

        fclose(tmp);
    }

    SUBCASE("directed weighted graph") {
        const char* gname = "dir_wgh_test.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate real general";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 7, 5, 8};
        int tmp_col[] = {0, 0, 0, 0, 1, 2, 2, 4, 7};
        int tmp_weight[] = {1, 2, 3, 5, 2, 4, 3, 2, 1};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, tmp_weight);

        fclose(tmp);
        tmp = fopen(gname, "r");
        query_gprops(gname, &gprops);

        CHECK_EQ(gprops.is_directed, true);
        CHECK_EQ(gprops.is_weighted, true);

        fclose(tmp);
    }
}

TEST_CASE("Test graph loading using mmio with some corner cases") {

    SUBCASE("edges missing for undirected graph without symmetric property") {
        const char* gname = "error1.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate pattern general";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 7, 5, 8};
        int tmp_col[] = {0, 0, 0, 0, 1, 2, 2, 4, 7};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, nullptr);

        fclose(tmp);
        tmp = fopen(gname, "r");

        CHECK_UNARY_FALSE(read_matrix_market_pattern(tmp, &coo));
        CHECK_NE(coo.nnz, 18);

        fclose(tmp);
    }

    SUBCASE("check out of bound indices") {
        const char* gname = "error2.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate pattern symmetric";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 9, 5, 8};
        int tmp_col[] = {0, -1, 0, 0, 1, 2, 2, 4, 7};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, nullptr);

        fclose(tmp);
        tmp = fopen(gname, "r");

        CHECK_UNARY(read_matrix_market_pattern(tmp, &coo));

        fclose(tmp);
    }

    SUBCASE("handling numbers near infinite") {
        const char* gname = "error3.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[] = "%%MatrixMarket matrix coordinate pattern symmetric";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 423423644, 2, 6, 5, 5, 8};
        int tmp_col[] = {0, 0, 0, 0, 1, 2, 2, 4, 729483202};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, nullptr);

        fclose(tmp);
        tmp = fopen(gname, "r");

        CHECK_UNARY(read_matrix_market_pattern(tmp, &coo));

        fclose(tmp);
    }

    SUBCASE("handling non-numeric values") {
        const char* gname = "error4.mtx";
        FILE *tmp = fopen(gname, "w");
        REQUIRE_UNARY(tmp);

        char tmp_mm_header[MM_HEADER] = "%%MatrixMarket matrix coordinate pattern symmetric";
        int tmp_header[] = {9, 9, 9};
        int tmp_row[] = {1, 3, 4, 5, 2, 6, 9, 5, 8};
        int tmp_col[] = {0, -1, 0, 0, 1, 2, 2, 4, 7};

        write_tmp_file(tmp, tmp_mm_header, tmp_header, tmp_row, tmp_col, nullptr);
        write_random_char_to_file(tmp, 3);

        fclose(tmp);
        tmp = fopen(gname, "r");

        CHECK_UNARY(read_matrix_market_pattern(tmp, &coo));

        fclose(tmp);
    }
}

TEST_CASE("Test sparse matrix storage format conversion from COO to CSR") {

    /*
     * Workspace setup for this test.
     */
    int rows[] = {1, 0, 3, 0, 4, 0, 5, 0, 2, 1, 6, 2, 7, 2, 5, 4, 8, 7};
    int cols[] = {0, 1, 0, 3, 0, 4, 0, 5, 1, 2, 2, 6, 2, 7, 4, 5, 7, 8};

    coo.nnz = 18;
    coo.nrows = 9;
    coo.rows = rows;
    coo.cols = cols;

    pcoo_to_pcsr(&coo, &csr);

    int expected_row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
    int expected_cols[] =
            {1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};

    CHECK_EQ(csr.nrows, 9);

    for (int i = 0; i < coo.nrows; i++) {
        CHECK_EQ(csr.row_offsets[i], expected_row_offsets[i]);
    }

    for (int i = 0; i < coo.nnz; i++) {
        CHECK_EQ(csr.cols[i], expected_cols[i]);
    }
}
