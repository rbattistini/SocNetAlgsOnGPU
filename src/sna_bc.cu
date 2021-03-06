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
#include "bc_ep_kernel.cuh"
#include "bc_statistics.h"
#include "bc_vp_kernel.cuh"
#include "bc_we_kernel.cuh"
#include "bc_we_kernel_nopitch.cuh"
#include "cl.h"
#include "cl_kernels.cuh"
#include "degree.h"
#include "matio.h"
#include <cli.cuh>

int main(int argc, char *argv[]) {

    params_t params;
    int err_code = parse_args(&params, argc, argv);

    if (err_code == EXIT_FAILURE) {
        return EXIT_FAILURE;
    } else if (err_code == EXIT_WHELP_OR_USAGE) {
        return EXIT_SUCCESS;
    }

    /*
     * Check if required arguments are provided.
     */
    if (params.input_file == 0) {
        ZF_LOGF("Input file required");
        return EXIT_FAILURE;
    }

    set_device(params.device_id);

    if (params.verbose) {
        zf_log_set_output_level(ZF_LOG_INFO);
    } else {
        zf_log_set_output_level(ZF_LOG_ERROR);
    }

    if (get_compute_capability_major() < 6) {
        ZF_LOGF("Atomic operations for doubles are available only for compute"
                "capability at least 6.x");
        return EXIT_FAILURE;
    }

    matrix_pcsr_t m_csr;
    matrix_pcoo_t m_coo;
    gprops_t gp;
    double tstart, tend, tstart_coo_to_csr, tend_coo_to_csr, tstart_cc,
            tend_cc, tstart_sub_ex, tend_sub_ex;

    gp.has_self_loops = params.self_loops_allowed;

    /*
     * Load matrix in COO format.
     */
    if (query_gprops(params.input_file, &gp) ||
        read_matrix(params.input_file, &m_coo, &gp)) {

        ZF_LOGF("Could not read matrix %s", params.input_file);
        return EXIT_FAILURE;
    }
    tstart_coo_to_csr = get_time();
    coo_to_csr(&m_coo, &m_csr);
    tend_coo_to_csr = get_time();

    /*
     * Extract the subgraph induced by vertices of the largest cc.
     */
    matrix_pcsr_t g;
    components_t ccs;

    tstart_cc = get_time();
    get_cc(&m_csr, &ccs);
    tend_cc = get_time();

    gp.is_connected = (ccs.cc_count == 1);
    if (!gp.is_connected) {
        tstart_sub_ex = get_time();
        get_largest_cc(&m_csr, &g, &ccs);
        tend_sub_ex = get_time();
        free_matrix_pcsr(&m_csr);
    } else {
        g = m_csr;
    }
    free_ccs(&ccs);

    tstart = get_time();
    auto degree = (int *) malloc(g.nrows * sizeof(int));
    compute_degrees_undirected(&g, degree);
    tend = get_time();

    /*
     * Print overview.
     */
    if(!params.quiet) {
        print_run_config(&params);
        print_gpu_overview(params.device_id);
        print_graph_properties(&gp);
        print_graph_overview(&g, degree);

        ZF_LOGI("COO to CSR executed in: %g s",
                tend_coo_to_csr - tstart_coo_to_csr);
        ZF_LOGI("Connected Component computation executed in: %g s",
                tend_cc - tstart_cc);
        ZF_LOGI("Subgraph extraction from largest cc executed in: %g s",
                tend_sub_ex - tstart_sub_ex);
        ZF_LOGI("Degree computation executed in: %g s",
                tend - tstart);
    } else {
        printf("File: %s, Technique: %d\n",
               params.input_file,
               params.technique);
    }

    /*
     * Allocate memory on the host.
     */
    int *dist;
    unsigned long long *sigma;
    double *bc_gpu, *delta, *cl_gpu;

    sigma = (unsigned long long *) malloc(g.nrows * sizeof(*sigma));
    dist = (int *) malloc(g.nrows * sizeof(*dist));
    bc_gpu = (double *) malloc(g.nrows * sizeof(*bc_gpu));
    cl_gpu = (double *) malloc(g.nrows * sizeof(*cl_gpu));
    delta = (double *) malloc(g.nrows * sizeof(*delta));

    if (sigma == 0 || dist == 0 || bc_gpu == 0 || delta == 0) {
        ZF_LOGF("Could not allocate memory");
        return EXIT_FAILURE;
    }

    /*
     * Memory allocation of the structure to which gathered statistics are
     * stored.
     */
    stats_t stats;
    stats.nedges_traversed = g.nrows * g.row_offsets[g.nrows];

    /*
     * Closeness centrality computation on the GPU.
     */
    compute_cl_gpu_p(&g, cl_gpu, &stats);

    /*
     * BC computation on the GPU.
     */
    ParStrategy technique = params.technique;
    switch (technique) {
        case work_efficient:
            compute_bc_gpu_wep(&g, bc_gpu, &stats);
//            compute_bc_gpu_we(&g, bc_gpu, &stats);
            break;
        case vertex_parallel:
            compute_bc_gpu_vpp(&g, bc_gpu, &stats);
            break;
        case edge_parallel:
            compute_bc_gpu_epp(&g, bc_gpu, &stats);
            break;
        default:
            ZF_LOGE("Invalid technique Id, cannot compute betweenness");
    }

    /*
     * BC and Closeness centrality computation on the CPU.
     */
    if (params.run_check) {
        auto bc_cpu = (double *) malloc(g.nrows * sizeof(double));
        if (bc_cpu == 0) {
            ZF_LOGF("Could not allocate memory");
            return EXIT_FAILURE;
        }

        tstart = get_time();
        compute_ser_bc_cpu(&g, bc_cpu, gp.is_directed);
        tend = get_time();
        stats.cpu_time = tend - tstart;

        double bc_error = check_score(g.nrows, bc_cpu, bc_gpu);

        auto cl_cpu = (double *) malloc(g.nrows * sizeof(double));
        if (cl_cpu == 0) {
            ZF_LOGF("Could not allocate memory");
            return EXIT_FAILURE;
        }

        tstart = get_time();
        compute_cl_cpu(&g, cl_cpu);
        tend = get_time();
        stats.cpu_time = tend - tstart;

        double cl_error = check_score(g.nrows, cl_cpu, cl_gpu);

        if(!params.quiet) {
            printf("Betweenness RMSE error: %g\n", bc_error);
            printf("Closeness RMSE error: %g\n", cl_error);
            free(bc_cpu);
            free(cl_cpu);
        }
    }

    /*
     * Dump scores and statistics if requested.
     */
    if (params.dump_scores != 0) {
        dump_scores(g.nrows, degree, bc_gpu, cl_gpu, params.dump_scores);
        free_params(&params);
    }

    if (params.dump_stats != 0) {
        append_stats(&stats, params.dump_stats, params.technique);
    } else if(!params.quiet){
        print_stats(&stats);
    }

    /*
     * Cleanup.
     */
    free(sigma);
    free(dist);
    free(delta);
    free(bc_gpu);

    if (!gp.is_connected)
        free_matrix_pcsr(&g);

    cudaSafeCall(cudaDeviceReset());
    return EXIT_SUCCESS;
}
