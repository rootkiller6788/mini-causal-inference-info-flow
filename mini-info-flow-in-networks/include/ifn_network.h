#ifndef IFN_NETWORK_H
#define IFN_NETWORK_H
#include "ifn_core.h"

/* Network structure and information metrics.
 * Characterizing how network topology shapes information flow. */

/* Network creation and topology */
IFN_CausalGraph* ifn_network_random_erdos(int n, double p);
IFN_CausalGraph* ifn_network_barabasi_albert(int n, int m0, int m);
IFN_CausalGraph* ifn_network_watts_strogatz(int n, int k, double beta);
IFN_CausalGraph* ifn_network_from_adjacency(const double* adj, int n);

/* Network metrics */
double ifn_clustering_coefficient(const IFN_CausalGraph* g);
double ifn_average_path_length(const IFN_CausalGraph* g);
int ifn_shortest_path(const IFN_CausalGraph* g, int from, int to,
    int* path, int max_len);
double ifn_modularity(const IFN_CausalGraph* g, const int* communities, int n_comm);
double ifn_assortativity(const IFN_CausalGraph* g);

/* Information-theoretic network measures */
double ifn_network_mutual_information(const IFN_CausalGraph* g,
    const double* signal, int n_nodes);
double ifn_network_transfer_entropy_sum(const IFN_CausalGraph* g,
    const double** signals, int n_nodes, int T, int k, int l);
double ifn_network_efficiency(const IFN_CausalGraph* g);
double ifn_network_redundancy(const IFN_CausalGraph* g);

/* Flow optimization on networks */
void ifn_max_flow(const IFN_CausalGraph* g, int source, int sink, double* flow);
double ifn_min_cut(const IFN_CausalGraph* g, int source, int sink);
void ifn_shortest_path_tree(const IFN_CausalGraph* g, int root, int* parent);

/* Network comparison */
double ifn_network_kl_divergence(const IFN_CausalGraph* g1, const IFN_CausalGraph* g2);
double ifn_spectral_distance(const IFN_CausalGraph* g1, const IFN_CausalGraph* g2);

#endif
