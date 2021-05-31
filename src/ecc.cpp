/****************************************************************************
 * @file ecc.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions to compute the diameter of an undirected graph, the
 * eccentricity of its vertices and to get the graph density.
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

#include "ecc.h"

int get_diameter(matrix_pcsr_t *g) {

    int max_distance;
    int *eccentricity = (int *) malloc(g->nrows * sizeof(*eccentricity));

    if (eccentricity == 0) {
        ZF_LOGF("Could not allocate memory");
        return -1;
    }

    get_vertices_eccentricity(g, eccentricity);
    max_distance = eccentricity[argmax(eccentricity, g->nrows)];
    free(eccentricity);

    return max_distance;
}

int get_vertices_eccentricity(matrix_pcsr_t *g, int *eccentricity) {

    auto *d = (int *) malloc(g->nrows * sizeof(int));
    if(d == 0) {
        ZF_LOGF("Could not allocate memory");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < g->nrows; i++) {
        fill(d, g->nrows, INT_MAX);
        BFS_visit(g, d, i);
        eccentricity[i] = d[argmax(d, g->nrows)];
    }

    free(d);
    return EXIT_SUCCESS;
}

double get_density(double nvertices, double nedges) {
    return nedges / ((nvertices * (nvertices - 1.0f)) / 2.0f);
}
