/****************************************************************************
 * @file cli.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Command line argument parsing functions using getopt
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

#pragma once
#ifndef CLI_H
#define CLI_H

#include "zf_log.h"
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#define EXIT_WHELP_OR_USAGE 2
#define NTECHNIQUES 3

/**
 * List all technique used for computing BC on the GPU.
 */
enum ParStrategy {
    vertex_parallel = 1,
    edge_parallel   = 2,
    work_efficient  = 3
};

/**
 * Lists all accepted command line arguments.
 *
 * verify, vf (pbgl run and RMSE error reporting)
 * dump_scores, ds (.csv file, vertex id and score)
 * append_stats, dst (.csv file, ph 1, 2, 3, runtime, teps)
 */
typedef struct params_t {
    int run_check;
    int verbose;
    int quiet;
    int self_loops_allowed;
    int device_id;
    ParStrategy technique;
    char *dump_scores;
    char *dump_stats;
    char *input_file;
} params_t;

int parse_args(params_t *p, int argc, char *argv[]);

/**
 * @brief Dump parameters for the current program configuration to stdout.
 *
 * @param p structure with running configuration parameters
 * @param fname file where the dump happens
 * @return 0 if successful, -1 if the stream was not closed correctly,
 * 1 if another error occurred
 */
void print_run_config(params_t *p);

void free_params(params_t *p);

#endif//CLI_H
