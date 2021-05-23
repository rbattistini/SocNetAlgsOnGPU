/****************************************************************************
 * @file bc_statistics.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Compute and store statistics such as runtime and teps for the
 * GPU BC computation algorithm.
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

#include "bc_statistics.h"

/*
 * REVIEW
 */
int dump_stats(stats_t *stats, char *fname) {

    if (fname == 0) {
        ZF_LOGE("No filename given");
        return EXIT_FAILURE;
    }

    FILE *f = fopen(fname, "w");

    if (stats->total_time == 0 || stats->load_time == 0 ||
        stats->bc_comp_time == 0 || stats->unload_time == 0 ||
        stats->nedges_traversed == 0 || stats->nrun == 0) {
        ZF_LOGE("Statistics not completely initialized");
        return EXIT_FAILURE;
    }

    if (f != 0) {

        fprintf(f, "\"Total Time\", \"Load Time\", \"Unload Time\","
                   " \"BC Comp Time\", \"Teps\"\n");

        for (int i = 0; i < stats->nrun; i++) {

            double teps = get_bc_teps(stats->nedges_traversed,
                                      stats->total_time[i]);

            fprintf(f, "%.2f, %.2f, %.2f, %.2f, %.2f\n",
                    stats->total_time[i],
                    stats->load_time[i],
                    stats->unload_time[i],
                    stats->bc_comp_time[i],
                    teps);
        }

    } else {
        ZF_LOGE("Failed to dump statistics");
        return EXIT_FAILURE;
    }

    return close_stream(f);
}

/*
 * REVIEW
 */
int dump_bc_scores(int nvertices, const float *bc_scores, char *fname) {

    if (fname == 0) {
        ZF_LOGE("No filename given");
        return EXIT_FAILURE;
    }

    FILE *f = fopen(fname, "w");

    if (bc_scores == 0) {
        ZF_LOGE("Bc scores not initialized");
        return EXIT_FAILURE;
    }

    if (f != 0) {
        fprintf(f, "\"Vertex Id\", \"Bc score\"\n");

        for (int i = 0; i < nvertices; i++) {
            fprintf(f, "%d, %.2f\n", i, bc_scores[i]);
        }

    } else {
        ZF_LOGE("Failed to create output file");
        return EXIT_FAILURE;
    }

    return close_stream(f);
}

void free_stats(stats_t *s) {
    free(s->total_time);
    free(s->load_time);
    free(s->unload_time);
    free(s->bc_comp_time);
    s->nrun = 0;
    s->nedges_traversed = 0;
}
