/****************************************************************************
 *
 * test_mod_bc.cpp
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

#include "spvb.h"
#include "graphs.h"
#include "tests.h"
#include "Snap.h"

static matrix_pcoo_t m_coo;
static gprops_t gprops;

TEST_CASE("Test SPVB BC computation on an undirected unweighted graph") {

    /*
     * Workspace setup for this test.
     */
    int rows[] = {1, 0, 3, 0, 4, 0, 5, 0, 2, 1, 6, 2, 7, 2, 5, 4, 8, 7};
    int cols[] = {0, 1, 0, 3, 0, 4, 0, 5, 1, 2, 2, 6, 2, 7, 4, 5, 7, 8};

    m_coo.nnz = 18;
    m_coo.nrows = 9;
    m_coo.rows = rows;
    m_coo.cols = cols;

    float expected_bc_scores[] =
            {17.0, 16.0, 17.0, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0};
    float *bc_scores, *partial_bc_scores;

    PUNGraph g = TUNGraph::New();
    for(int i = 0; i < m_coo.nrows; i++) {
        g->AddNode(i);
    }

    for(int i = 0; i < m_coo.nnz; i++) {
        g->AddEdge(m_coo.rows[i], m_coo.cols[i]);
    }

    int *degree = (int*) malloc(m_coo.nrows * sizeof(*degree));
    REQUIRE_UNARY(degree);

    compute_degrees_undirected(&m_coo, degree);
    gprops.is_directed = false;

    //verified both g and degree

    bc_scores = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
    partial_bc_scores = (float*) malloc(m_coo.nrows * sizeof(*bc_scores));
    REQUIRE_UNARY(bc_scores);
    REQUIRE_UNARY(partial_bc_scores);

    spvb(g, degree, bc_scores, partial_bc_scores, gprops.is_directed);

    for(int i = 0; i < m_coo.nrows; i++) {
        REQUIRE_EQ(bc_scores[i], expected_bc_scores[i]);
    }

    free(bc_scores);
}
