/****************************************************************************
 * @file sna_bc.cu
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
 ****************************************************************************/

#include "bc.h"
#include "bc_kernels.cuh"
#include "bc_kernels_pitched.cuh"
#include "matio.h"
#include "timing.cuh"
#include <cli.h>

int main(int argc, char *argv[]) {

    params_t params;
    int err_code = parse_args(&params, argc, argv);

    if (err_code == EXIT_FAILURE) {
        return EXIT_FAILURE;
    } else if(err_code == EXIT_WHELP_OR_USAGE) {
        close_stream(stdin);
        close_stream(stdout);
        close_stream(stderr);
        return EXIT_SUCCESS;
    }

    /*
     * Check if required arguments are provided.
     */
    if (params.input_file == 0) {
        ZF_LOGF("Input file required");
        return EXIT_FAILURE;
    }

    if (params.nrun <= 0) {
        ZF_LOGF("Number of runs required");
        return EXIT_FAILURE;
    }

    print_run_config(&params);
    set_device(params.device_id);
    print_gpu_overview(params.device_id);

    if (params.verbose) {
        zf_log_set_output_level(ZF_LOG_INFO);
    } else {
        zf_log_set_output_level(ZF_LOG_ERROR);
    }

    if (get_compute_capability_major() < 6) {
        ZF_LOGW("This program has been tested only with devices "
                "with compute capability at least 6.x");
    }

    matrix_pcsr_t m_csr;
    matrix_pcoo_t m_coo;
    gprops_t gp;
    EventTimer chrono;

    gp.has_self_loops = params.self_loops_allowed;

    /*
     * Load matrix in COO format.
     */
    if (query_gprops(params.input_file, &gp) ||
        read_matrix(params.input_file, &m_coo, &gp)) {

        ZF_LOGF("Could not read matrix %s", params.input_file);
        return EXIT_FAILURE;
    }
    chrono.start();
    coo_to_csr(&m_coo, &m_csr);
    chrono.stop();

    chrono.log("COO to CSR");

    /*
     * Extract the subgraph induced by vertices of the largest cc.
     */
    matrix_pcsr_t g;
    components_t ccs;

    chrono.start();
    get_cc(&m_csr, &ccs);
    chrono.stop();

    gp.is_connected = (ccs.cc_count == 1);

    chrono.log("Get cc");

    if (params.verbose) {
        print_separator();
    }
    print_gprops(&gp);

    if (!gp.is_connected) {
        chrono.start();
        get_largest_cc(&m_csr, &g, &ccs);
        chrono.stop();
        free_matrix_pcsr(&m_csr);
    } else {
        g = m_csr;
    }
    free_ccs(&ccs);

    chrono.log("Extract largest cc");

    /*
     * Allocate memory on the host.
     */
    int *dist;
    unsigned long long *sigma;
    float *bc_gpu, *delta;

    sigma = (unsigned long long *) malloc(g.nrows * sizeof(*sigma));
    dist = (int *) malloc(g.nrows * sizeof(*dist));
    bc_gpu = (float *) malloc(g.nrows * sizeof(*bc_gpu));
    delta = (float *) malloc(g.nrows * sizeof(*delta));

    if (sigma == 0 || dist == 0 || bc_gpu == 0 || delta == 0) {
        ZF_LOGF("Could not allocate memory");
        return EXIT_FAILURE;
    }

    /*
     * Memory allocation of the structure to which gathered statistics are
     * stored.
     */
    stats_t stats;
    if (params.dump_stats != 0) {

        stats.total_time = (float*) malloc(stats.nrun * sizeof(float));
        stats.bc_comp_time = (float*) malloc(stats.nrun * sizeof(float));
        stats.unload_time = (float*) malloc(stats.nrun * sizeof(float));
        stats.load_time = (float*) malloc(stats.nrun * sizeof(float));

        if(stats.total_time == 0 || stats.bc_comp_time == 0 ||
            stats.unload_time == 0 || stats.load_time == 0) {
            ZF_LOGF("Could not allocate memory");
            return EXIT_FAILURE;
        }
    }

    /*
     * BC computation on the GPU.
     */
    for (int run = 0; run < params.nrun; run++) {

        chrono.start();
        compute_bc_gpu(&g, bc_gpu);
        chrono.stop();

        float time_elapsed = chrono.elapsed();
        ZF_LOGI("BC on GPU: %f", time_elapsed);

        chrono.start();
        compute_bc_gpu_wpitched(&g, bc_gpu);
        chrono.stop();
        time_elapsed = chrono.elapsed();
        ZF_LOGI("BC on GPU with pitched memory: %f", time_elapsed);

        if (params.dump_stats != 0) {
            stats.total_time[run] = time_elapsed;
        }
    }

    /*
     * BC computation on the CPU.
     */
    if (params.run_check) {
        auto bc_cpu = (float *) malloc(g.nrows * sizeof(float));
        if (bc_cpu == 0) {
            ZF_LOGF("Could not allocate memory");
            return EXIT_FAILURE;
        }

        chrono.start();
        compute_par_bc_cpu(&g, bc_cpu);
        chrono.stop();

        chrono.log("BC on CPU");

        double error = check_bc(g.nrows, bc_cpu, bc_gpu);
        printf("RMSE error: %g\n", error);
        free(bc_cpu);
    }

    /*
     * Dump scores and statistics if requested.
     */
    if (params.dump_scores != 0) {
        dump_bc_scores(g.nrows, bc_gpu, params.dump_scores);
        free_params(&params);
    }

    if (params.dump_stats != 0) {
        dump_stats(&stats, params.dump_stats);
        free_stats(&stats);
    }

    /*
     * Cleanup.
     */
    free(sigma);
    free(dist);
    free(delta);
    free(bc_gpu);

    if(!gp.is_connected)
        free_matrix_pcsr(&g);

    close_stream(stdin);
    close_stream(stdout);
    close_stream(stderr);

    return EXIT_SUCCESS;
}
