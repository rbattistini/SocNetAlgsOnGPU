/****************************************************************************
 * @file test_serial_bc.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
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

#include "tests.h"
#include "bc.h"
#include "graphs.h"

static matrix_pcsr_t csr;
static gprops_t gprops;

TEST_CASE("Test BC computation on an undirected unweighted graph") {

    /*
     * Workspace setup for this test.
     */
    int nrows = 9;
    int row_offsets[] = {0, 4, 6, 9, 10, 12, 14, 15, 17, 18};
    int cols[] = {1, 3, 4, 5, 0, 2, 1, 6, 7, 0, 0, 5, 0, 4, 2, 2, 8, 7};

    float expected_bc_scores[] =
            {17.0, 16.0, 17.0, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0};
    float *bc_scores;

    csr.nrows = nrows;
    csr.row_offsets = row_offsets;
    csr.cols = cols;

    gprops.is_directed = false;

    bc_scores = (float*) malloc(csr.nrows * sizeof(*bc_scores));
    REQUIRE_UNARY(bc_scores);

    get_vertex_betweenness(&csr, bc_scores, gprops.is_directed);

    for(int i = 0; i < nrows; i++) {
        REQUIRE_EQ(bc_scores[i], expected_bc_scores[i]);
    }

    free(bc_scores);
}
