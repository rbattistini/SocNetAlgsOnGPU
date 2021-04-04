/****************************************************************************
 *
 * betweenness.cpp - Serial algorithm for computing betweenness centrality
 *
 * Based on ...
 *
 * Last updated in 2021 by Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * ---------------------------------------------------------------------------
 *
 * Compile with:
 * g++ betweenness.cpp -o betweenness
 *
 * Run with:
 * ./betweenness [input_filename]
 *
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include "Snap.h"

using namespace TSnap;

/**
 * Print graph's basic statistics, namely the name of the graph, the number
 * of edges and of vertices and whether the graph is empty and connected.
 *
 * @tparam PGraph type of graph, directed or undirected
 * @param s pointer to the name of the graph
 * @param G pointer to the object representing the graph
 */
template <class PGraph>
void PrintGStats(const char *s, PGraph G) {
    printf("G's %s statistics: \n"
           "nodes:\t %d\n"
           "edges:\t %d\n"
           "isEmpty:\t %s\n"
           "isConnected:\t %s\n",
           s, G->GetNodes(), G->GetEdges(),
           G->Empty() ? "yes" : "no",
           IsConnected(G) ? "yes" : "no");
}

int main( int argc, char *argv[] ) {

    char* fName;
    PUNGraph g;

    if (argc > 2) {
        std::cerr << "Usage: " << argv[0]<< " [input_filename]" << std::endl;
        return EXIT_FAILURE;
    }

    /*
     * Load SNAP graph stored in .DAT format, as list of edges where nodes are
     * ids.
     */
    fName = argv[1];
    g = LoadEdgeList<PUNGraph>(fName);


    return EXIT_SUCCESS;
}
