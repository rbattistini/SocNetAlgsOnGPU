extern "C" {
#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"
}
#include "matio.h"
#include "matds.h"
#include "graphs.h"

int dump_scores(long nvertices, const double *bc_scores, char *fname) {

    if (fname == 0) {
        fprintf(stderr, "No filename given");
        return EXIT_FAILURE;
    }

    FILE *f = fopen(fname, "w");

    if (bc_scores == 0) {
        fprintf(stderr, "Betweenness centrality scores not initialized");
        return EXIT_FAILURE;
    }

    if (f != 0) {
        fprintf(f, "\"Vertex Id\", \"Betweenness\"\n");

        for (int i = 0; i < nvertices; i++) {
            fprintf(f, "%d, %.2f\n", i, bc_scores[i]);
        }

    } else {
        fprintf(stderr, "Failed to create output file");
        return EXIT_FAILURE;
    }

    fclose(f);
    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {

    char *infilename, *outfilename, *graph_type;
    FILE *fp;
    graph_t *g;
    double *bc_scores;
    int curArgIndex;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                        " (-graph <graph type> -outfile <output filename>)\n\n",
                "eval_vertex_betweenness");

        usage_graph_options();
        exit(-1);
    }

    curArgIndex = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));

    strcpy(outfilename, "output.csv");

    while (curArgIndex < argc) {

        if (strcmp(argv[curArgIndex], "-infile") == 0) {
            strcpy(infilename, argv[++curArgIndex]);
        }

        if (strcmp(argv[curArgIndex], "-outfile") == 0) {
            strcpy(outfilename, argv[++curArgIndex]);
        }

        if (strcmp(argv[curArgIndex], "-graph") == 0) {
            strcpy(graph_type, argv[++curArgIndex]);
        }
        curArgIndex++;
    }

    fp = fopen(infilename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not open input file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);

    fp = fopen(outfilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not write to output file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);

    graph_ext_check(infilename, graph_type);

    fprintf(stdout, "\n");
    fprintf(stdout, "Input Graph File    : %s\n", infilename);
    fprintf(stdout, "Output Graph File   : %s\n\n", outfilename);

    /* Step 2: Generate graph */
    matrix_pcsr_t m_csr;
    matrix_pcoo_t m_coo;
    gprops_t gp;

    if (query_gprops(infilename, &gp) ||
        read_matrix(infilename, &m_coo, &gp)) {

        fprintf(stderr, "Could not read matrix %s", infilename);
        return EXIT_FAILURE;
    }

    coo_to_csr(&m_coo, &m_csr);

    matrix_pcsr_t g_tmp;
    components_t ccs;

    get_cc(&m_csr, &ccs);

    gp.is_connected = (ccs.cc_count == 1);
    if (!gp.is_connected) {
        get_largest_cc(&m_csr, &g_tmp, &ccs);
        free_matrix_pcsr(&m_csr);
    } else {
        g_tmp = m_csr;
    }
    free_ccs(&ccs);

    long nnz = g_tmp.row_offsets[g_tmp.nrows];
    auto rows = (int *) malloc(nnz * sizeof(int));
    expand_row_pointer(g_tmp.nrows, g_tmp.row_offsets, rows);

    g = (graph_t *) malloc(sizeof(graph_t));

    g->n = g_tmp.nrows;
    g->m = nnz;
    g->endV = g_tmp.cols;
    g->numEdges = g_tmp.row_offsets;
    g->edge_id = rows;
    g->undirected = 1;
    g->weight_type = 0;
    g->zero_indexed = 1;

    /* Step 3: Run algorithm */
    bc_scores = (double *) calloc(g->n, sizeof(double));
    vertex_betweenness_centrality(g, bc_scores, g->n);

    /* Step 4: Write results to output file */
    for (int k = 0; k < g->n; k++)
        bc_scores[k] /= 2;

    dump_scores(g->n, bc_scores, outfilename);

    /* Step 5: Clean up */
    free(bc_scores);
    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);

    return 0;
}
