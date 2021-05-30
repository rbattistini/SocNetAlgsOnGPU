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

double check_score(int nrows,
                   const float *score_cpu,
                   const float *score_gpu) {
    double error = 0;
    double max_error = 0;

    for (int i = 0; i < nrows; i++) {

        double current_error = abs(score_cpu[i] - score_gpu[i]);
        error += current_error * current_error;

        if (current_error > max_error) {
            max_error = current_error;
        }
    }
    error = error / (float) nrows;
    error = sqrt(error);
    return error;
}

void print_stats(stats_t *s) {

    print_separator();

    printf("%15s|%15s|%15s|%15s|%15s|\n", "Total Time", "Load Time",
           "Unload Time", "BC Comp Time", "Teps");

    print_separator();

    double teps = get_bc_teps(s->nedges_traversed, s->total_time);

    printf("%15g|%15g|%15g|%15g|%15g|\n",
           s->total_time,
           s->load_time,
           s->unload_time,
           s->bc_comp_time,
           teps);
}

int append_stats(stats_t *stats, char *fname) {

    if (fname == 0) {
        ZF_LOGE("No filename given");
        return EXIT_FAILURE;
    }

    FILE *f = fopen(fname, "a");

    if (stats->total_time == 0 || stats->load_time == 0 ||
        stats->bc_comp_time == 0 || stats->unload_time == 0 ||
        stats->nedges_traversed == 0) {
        ZF_LOGE("Statistics not completely initialized");
        return EXIT_FAILURE;
    }

    if (f != 0) {

//        fprintf(f, "\"Total Time\", \"Load Time\", \"Unload Time\","
//                   " \"BC Comp Time\", \"Teps\"\n");

        double teps = get_bc_teps(stats->nedges_traversed,
                                  stats->total_time);

        fprintf(f, "%.2f, %.2f, %.2f, %.2f, %.2f\n",
                stats->total_time,
                stats->load_time,
                stats->unload_time,
                stats->bc_comp_time,
                teps);

    } else {
        ZF_LOGE("Failed to dump statistics");
        return EXIT_FAILURE;
    }

    return close_stream(f);
}

int dump_scores(int nvertices,
                const int *degree_scores,
                const float *bc_scores,
                const float *cl_scores,
                char *fname) {

    if (fname == 0) {
        ZF_LOGE("No filename given");
        return EXIT_FAILURE;
    }

    FILE *f = fopen(fname, "w");

    if (degree_scores == 0) {
        ZF_LOGE("Degree centrality scores not initialized");
        return EXIT_FAILURE;
    }

    if (bc_scores == 0) {
        ZF_LOGE("Betweenness centrality scores not initialized");
        return EXIT_FAILURE;
    }

    if (cl_scores == 0) {
        ZF_LOGE("Closeness centrality scores not initialized");
        return EXIT_FAILURE;
    }

    if (f != 0) {
        fprintf(f, "\"Vertex Id\", \"Degree\", \"Betweenness\","
                   " \"Closeness\"\n");

        for (int i = 0; i < nvertices; i++) {
            fprintf(f, "%d, %d, %.2f, %.2f\n", i,
                    degree_scores[i],
                    bc_scores[i],
                    cl_scores[i]);
        }

    } else {
        ZF_LOGE("Failed to create output file");
        return EXIT_FAILURE;
    }

    return close_stream(f);
}
