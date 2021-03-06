/****************************************************************************
 * @file cli.cu
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

#include "cli.cuh"

typedef struct commands_t {
    const char *cmd_name;
    const char *cmd_descr;
} commands_t;

static void print_usage(char *app_name) {
    printf("Usage:\n %s\t[-i|--input file] [-t|--technique] [-b|--dump-scores file] \n"
           "\t\t[-s|--dump-stats file] [-v|--verbose] [-c|--check]\n"
           "\t\t[-wsl|--wself-loops] [-d|--device] [-q|--quiet]\n"
           "\t\t[-u|--usage] ][-h|--help]\n",
           app_name);
}

static void print_help() {

    const int nopt = 11;
    static struct commands_t cmds[nopt] = {
            {"(i) input \t= <filename>\t",
                    "input matrix market file"},
            {"(b) dump-scores = <filename>\t",
                    "dump computed bc scores to <filename>"},
            {"(s) dump-stats \t= <filename>\t",
                    "dump stats of the GPU algorithm to <filename>"},
            {"(v) verbose\t\t\t",
                    "print info messages and errors"},
            {"(q) quiet\t\t\t",
                    "print only errors if they occur"},
            {"(c) check\t\t\t",
                    "compare GPU and CPU parallel algorithms"},
            {"(l) wself-loops\t\t\t",
                    "don't remove self loops from the input graph"},
            {"(t) technique\t\t\t",
                    "technique used to distribute work among GPU threads"},
            {"(d) device\t\t\t",
                    "set the device id of the GPU, 0 is the default"},
            {"(u) usage\t\t\t",
                    "print usage of the command"},
            {"(h) help\t\t\t",
                    "print a summary of available commands"}
    };

    printf("Available options:\n\n");

    for (int i = 0; i < nopt; i++)
        printf("%s%s\n", cmds[i].cmd_name, cmds[i].cmd_descr);

    print_separator();

    printf("Available techniques are: \n\n");
    printf("(1) Vertex Parallel\n");
    printf("(2) Edge Parallel\n");
    printf("(3) Work efficient\n");
}

/**
 * @brief Parse char array to long with error handling.
 *
 * @see https://en.cppreference.com/w/c/string/byte/strtol
 *
 * @param p pointer to the null-terminated byte string to be interpreted
 * @param endp pointer to a pointer to character
 * @param base base of the interpreted integer value
 * @return the parsed number or -1 if unsuccessful
 */
static long strtol_wcheck(const char *p, char *endp, int base) {
    int done = 0;
    long i = -1;

    while (!done) {
        /*
         * errno can be set to any non-zero value by a library function call
         * regardless of whether there was an error, so it needs to be cleared
         * in order to check the error set by strtol.
         */
        errno = 0;

        i = strtol(p, &endp, base);
        p = endp;

        if (errno == ERANGE) {
            ZF_LOGE("Range error occurred");
            return -1;
        }

        if (p == endp)
            done = 1;
    }

    return i;
}

/**
 * @brief Concatenate two strings. Uses heap memory.
 *
 * @param src first string to concatenate
 * @param dest second string to concatenate
 * @return the concatenated string
 */
static char *concat(const char *src, const char *dest) {

    size_t res_length = strlen(src) + strlen(dest) + 1;
    char *result = (char *) malloc(res_length * sizeof(char));

    if (result == 0)
        ZF_LOGF("Could not allocate memory");

    strncpy(result, src, res_length);
    snprintf(&result[strlen(result)], res_length, "%s", dest);
    return result;
}

static const char* get_technique_from_id(ParStrategy technique) {
    switch (technique) {
        case work_efficient:
            return "Work Efficient";
        case vertex_parallel:
            return "Vertex Parallel";
        case edge_parallel:
            return "Edge Parallel";
        default:
            ZF_LOGE("Invalid technique");
            return 0;
    }
}

int parse_args(params_t *params, int argc, char *argv[]) {

    int verbose = 0;
    int run_check = 0;
    int show_help = 0;
    int self_loops_allowed = 0;
    int show_usage = 0;
    int quiet = 0;

    char *technique = 0;
    char *dump_scores = 0;
    char *dump_stats = 0;
    char *input_file = 0;
    char *device_id = 0;
    int index;
    int cmd;

    static struct option long_options[] =
            {
                    {"verbose",     no_argument,       0, 'v'},
                    {"quiet",       no_argument,       0, 'q'},
                    {"check",       no_argument,       0, 'c'},
                    {"help",        no_argument,       0, 'h'},
                    {"wself-loops", no_argument,       0, 'l'},
                    {"usage",       no_argument,       0, 'u'},
                    {"device",      no_argument,       0, 'd'},
                    {"dump-scores", required_argument, 0, 'b'},
                    {"dump-stats",  required_argument, 0, 's'},
                    {"technique",   required_argument, 0, 't'},
                    {"input",       required_argument, 0, 'i'},
                    {0, 0,                             0, 0}
            };

    while (true) {

        int option_index = 0;
        cmd = getopt_long(argc, argv, "t:b:s:i:d:uvchql", long_options,
                          &option_index);

        /*
         * Detect the end of the options.
         */
        if (cmd == -1)
            break;

        switch (cmd) {
            case 'h':
                show_help = 1;
                break;
            case 'u':
                show_usage = 1;
                break;
            case 'q':
                quiet = 1;
                break;
            case 'l':
                self_loops_allowed = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 't':
                technique = optarg;
                break;
            case 'c':
                run_check = 1;
                break;
            case 'i':
                input_file = optarg;
                break;
            case 'b':
                dump_scores = optarg;
                break;
            case 's':
                dump_stats = optarg;
                break;
            case 'd':
                device_id = optarg;
                break;
            case '?':
                // getopt_long already printed an error message.
                break;
            default:
                return EXIT_FAILURE;
        }
    }

    /*
     * Whether to print usage and terminate execution.
     */
    if (show_usage) {
        print_usage(argv[0]);
        return EXIT_WHELP_OR_USAGE;
    }

    /*
     * Whether to print an help screen and terminate execution.
     */
    if (show_help) {
        print_help();
        return EXIT_WHELP_OR_USAGE;
    }

    /*
     * Whether to allow evaluation and printing of logger messages with level
     * other than ERROR.
     */
    params->verbose = verbose;

    /*
     * Whether to print overview and info messages regarding the device, the
     * input graph given and logger messages with level other than ERROR.
     */
    params->quiet = quiet;

    /*
     * Whether the reference CPU algorithm should be run and compared to the
     * GPU one in terms of runtime and RMSE error on the bc scores obtained.
     */
    params->run_check = run_check;

    /*
     * Whether self loops are allowed or must be removed.
     */
    params->self_loops_allowed = self_loops_allowed;

    /*
     * Set the id of the device to be used. The default is 0.
     */
    if (device_id != 0) {
        int tmp_id = (int) (strtol_wcheck(device_id, 0, 10));
        if (tmp_id < 0 || tmp_id > get_device_count() - 1) {
            ZF_LOGF("Invalid device id: min is 0, max is %d",
                    get_device_count() - 1);
            return EXIT_FAILURE;
        }
        params->device_id = tmp_id;
    } else {
        params->device_id = 0; // default is 0
    }

    /*
     * Whether to dump bc scores to a file.
     */
    params->dump_scores =
            (dump_scores == 0) ? dump_scores : concat(dump_scores, ".csv");

    /*
     * Whether to dump statistics to stdout or to a file.
     */
    params->dump_stats =
            (dump_stats == 0) ? dump_stats : concat(dump_stats, ".csv");

    /*
     * Define the approach to be used for work-distribution among GPU's threads.
     */
    if (technique != 0) {
        int tmp_technique = (int) (strtol_wcheck(technique, 0, 10));
        if (tmp_technique <= 0 || tmp_technique > NTECHNIQUES) {
            ZF_LOGF("Invalid technique chosen");
            return EXIT_FAILURE;
        }
        params->technique = (ParStrategy) tmp_technique;
    } else {
        ZF_LOGF("Invalid technique given in input");
        return EXIT_FAILURE;
    }

    /*
     * Input file.
     */
    params->input_file = input_file;

    /*
     * Print any remaining command line arguments (not valid options).
     */
    for (index = optind; index < argc; index++)
        ZF_LOGW("Non-recognized argument %s", argv[index]);

    return EXIT_SUCCESS;
}

void print_run_config(params_t *p) {

    print_separator();

    const char *output =
            (p->verbose) ? "verbose" : (p->quiet) ? "quiet" : "normal";

    const char *technique = get_technique_from_id(p->technique);

    printf("Run configuration:\n\n");
    printf("\tInput graph: \t\t%s\n", p->input_file);
    printf("\tStatistic file: \t%s\n", p->dump_stats);
    printf("\tBC scores file: \t%s\n", p->dump_scores);
    printf("\tTechnique: \t\t%s\n", technique);
    printf("\tDevice id: \t\t%d\n", p->device_id);
    printf("\tOutput: \t\t%s\n", output);
    printf("\tWith verification: \t%s\n",
           (p->run_check) ? "enabled" : "disabled");

    print_separator();
}

void free_params(params_t *p) {

    if (p->dump_scores != 0)
        free(p->dump_scores);

    if (p->dump_stats != 0)
        free(p->dump_stats);
}
