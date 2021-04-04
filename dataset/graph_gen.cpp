/****************************************************************************
 *
 * graph_gen.cpp - Graph generation
 *
 * Script to generate synthetic undirected and unweighted graphs datasets
 * using the SNAP Library.
 *
 * Based on examples available in the SNAP documentation.
 *
 * Last updated in 2021 by Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * ---------------------------------------------------------------------------
 *
 * Compile with:
 * g++ graph_gen.cpp -o graph_gen
 *
 * Run with:
 * ./graph_gen [filename] [graph_type] [number_vertices]
 *
 ****************************************************************************/

#include <iostream>
#include <cstdlib>
#include <Snap.h>

enum GraphType {
    Random     = 1, // Erdős-Rényi random graph
    SmallWorld = 2, // Watts-Strogatz Small World graph
    ScaleFree  = 3  // Barabási-Albert Scale Free graph
};

using namespace TSnap;

/*
 * Basic configuration parameters.
 */
const int min_vertices = 20;
const int max_vertices = 1000;

int binCoeff(int n, int k) {
    if (k == 0 || k == n)
        return 1;
    else
        return binCoeff(n - 1, k - 1) + binCoeff(n - 1, k);
}

double fact(double n)
{
    if (n == 0)
        return (1);
    else
        return (n * fact(n - 1));
}

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
    printf("statistics of undirected, unweighted graph in %s: \n"
           "nodes:\t\t %d\n"
           "edges:\t\t %d\n"
           "isEmpty:\t %s\n"
           "isConnected:\t %s\n",
           s, G->GetNodes(), G->GetEdges(),
           G->Empty() ? "yes" : "no",
           IsConnected(G) ? "yes" : "no");
}

int main( int argc, char *argv[] )
{
    int v, gType;
    char* fName;
    TInt::Rnd.PutSeed(0); // random number generator
    PUNGraph g;

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]<< " [filename] [graph_type] "
                                            "[number_vertices]" << std::endl;
        return EXIT_FAILURE;
    }

    fName = argv[1];
    gType = (int) strtol(argv[2], nullptr, 10);
    v = (int) strtol(argv[3], nullptr, 10);

    if(v < min_vertices || v > max_vertices) {
        std::cerr << "A graph must have at least 20 nodes and no more "
                     "than 1000 nodes" << std::endl;
        return EXIT_FAILURE;
    }

    switch (gType) {
        case Random: {
            /*
             * Number of edges, computed so that there is a high probability of
             * generating a connected random graph.
             */
            const int e = (int) (binCoeff(v, 2) * (log(v) / v)) + 2;
            g = GenRndGnm<PUNGraph>(v, e, false, TInt::Rnd);
            break;
        }
        case SmallWorld: {
            // number of neighbours to which each node is connected (Out-Degree)
            const int k = (int) (0.8 * v);
            // rewiring probability
            const float p = 0.03;
            g = GenSmallWorld(v, k, p, TInt::Rnd);
            break;
        }
        case ScaleFree: {
            // number of neighbours to which each node is connected (Out-Degree)
            const int k = (int) (0.8 * v);
            g = GenPrefAttach(v, k);
            break;
        }
        default:
            std::cerr << "Graph types accepted are: \n"
                         "[1]: Random graph \n"
                         "[2]: Small World graph \n"
                         "[3]: Scale Free graph"<< std::endl;
            return EXIT_FAILURE;
    }

    /*
     * Print graph basic statistics.
     */
    PrintGStats(fName, g);

    /*
     * Write SNAP graph to a DAT file, as list of edges where nodes are ids.
     */
    const char *description =
            "Randomly generated graph for algorithms benchmarking";
    SaveEdgeList(g, fName, description);

    return EXIT_SUCCESS;
}
