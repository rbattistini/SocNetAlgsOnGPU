/****************************************************************************
 * @file graph_gen.cpp
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Program to generate synthetic undirected and unweighted graphs
 * dataset using the Boost Graph Library. Based on examples available in the
 * BGL documentation.
 *
 * @see https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/erdos_renyi_generator.html
 * @see https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/plod_generator.html
 * @see https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/small_world_generator.html
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
 * ---------------------------------------------------------------------------
 *
 * Compile with:
 * g++ -I/usr/include/boost graph_gen_boost.cpp -o graph_gen_boost
 *
 * Run with:
 * ./graph_gen_boost [filename] [graph_type] [number_vertices]
 *
 ****************************************************************************/

#include "mmio.h"
#include <cstdlib>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/random/linear_congruential.hpp>

/**
 * @brief Parsing with error handling.
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
            fprintf(stderr, "Range error occurred");
            return -1;
        }

        if (p == endp)
            done = 1;
    }

    return i;
}

enum GraphType {
    Random     = 1, // Erdős-Rényi random graph
    SmallWorld = 2, // Watts-Strogatz Small World graph
    ScaleFree  = 3  // Barabási-Albert Scale Free graph
};

typedef boost::adjacency_list<> Graph;

static int write_bgl_to_mtx(FILE *f, const Graph &g) {

    typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;

    /*
     * Write the banner.
     */
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_pattern(&matcode);
    mm_set_symmetric(&matcode);

    if (mm_write_banner(f, matcode)) {
        fprintf(stderr, "Error writing mm banner");
        return EXIT_FAILURE;
    }

    /*
     * Write the header and the values.
     */
    if (mm_write_mtx_crd_size(f,
                              boost::num_vertices(g),
                              boost::num_vertices(g),
                              boost::num_edges(g))) {
        fprintf(stderr, "Error writing header");
        return EXIT_FAILURE;
    }

    IndexMap index = get(boost::vertex_index, g);
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        if(fprintf(f, "%lu %lu\n", index[source(*ei, g)],
                    index[target(*ei, g)]) < 0) {
            fprintf(stderr, "Error writing values");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int main( int argc, char *argv[] ) {

    typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
    typedef boost::small_world_iterator<boost::minstd_rand, Graph> SWGen;
    typedef boost::plod_iterator<boost::minstd_rand, Graph> SFGen;

    boost::minstd_rand gen;
    int v;
    int gtype;
    Graph g = 0;
    char* fname;

    /*
     * Basic configuration parameters.
     */
    const int min_vertices = 20;
    const int max_vertices = 1000;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s [filename] [graph_type] [number_vertices]\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    fname = argv[1];
    gtype = (int) strtol_wcheck(argv[2], nullptr, 10);
    v = (int) strtol_wcheck(argv[3], nullptr, 10);

    if (v < min_vertices || v > max_vertices) {
        fprintf(stderr, "A graph must have at least 20 nodes and no more "
                "than 1000 nodes");
        return EXIT_FAILURE;
    }

    switch (gtype) {
        case Random: {
            // probability of an edge between a pair of nodes
            const float p = 0.05;
            g = Graph(ERGen(gen, v, p), ERGen(), v);
            break;
        }
        case SmallWorld: {
            // number of neighbours to which each node is connected
            const int k = 17 * v;
            // rewiring probability
            const float p = 0.03;
            g = Graph(SWGen(gen, v, k, p), SWGen(), v);
            break;
        }
        case ScaleFree: {
            // controls how steeply is the power law curve
            const float alpha = 2.5;
            // determines the average degree of vertices
            const int beta = 1000;
            g = Graph(SFGen(gen, v, alpha, beta), SFGen(), v);
            break;
        }
        default:
            fprintf(stderr, "Graph types accepted are: \n"
                         "[1]: Random graph \n"
                         "[2]: Small World graph \n"
                         "[3]: Scale Free graph\n");
            return EXIT_FAILURE;
    }

    /*
     * Write BGL graph to a Matrix Market file.
     */
    FILE *f = fopen(fname, "w");

    if (write_bgl_to_mtx(f, g))
        return EXIT_FAILURE;

    fclose(f);
    return EXIT_SUCCESS;
}
