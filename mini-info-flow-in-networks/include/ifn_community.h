#ifndef IFN_COMMUNITY_H
#define IFN_COMMUNITY_H
#include "ifn_core.h"

/* Community detection in information flow networks.
 * Identifies modules where information flows more densely internally
 * than externally. Communities = functional processing units.
 *
 * Algorithms:
 *   - Girvan-Newman: iterative edge betweenness removal
 *   - Louvain: greedy modularity optimization
 *   - Infomap: information-theoretic (minimize random walk description length)
 *
 * Ref: Girvan & Newman (2002) PNAS 99, 7821
 *      Blondel et al. (2008) J. Stat. Mech. P10008
 *      Rosvall & Bergstrom (2008) PNAS 105, 1118
 */

/* Girvan-Newman algorithm: remove edges with highest betweenness.
 * Returns community assignments (array of n_nodes, values 0..n_communities-1).
 * n_communities is set to the number of detected communities. */
void ifn_girvan_newman(const IFN_CausalGraph* g,
    int* community_labels, int* n_communities);

/* Louvain method: greedy modularity maximization.
 * Two-phase iterative algorithm. Phase 1: local node movement.
 * Phase 2: community aggregation. Repeats until no improvement.
 * Returns modularity Q of final partition. */
double ifn_louvain(const IFN_CausalGraph* g,
    int* community_labels, int max_iter);

/* Infomap: information-theoretic community detection.
 * Minimizes the map equation L(M) = q*H(Q) + sum p_i*H(P_i)
 * where L(M) is the description length of a random walker.
 * Returns minimum description length. */
double ifn_infomap(const IFN_CausalGraph* g,
    int* community_labels, int max_iter);

/* Edge betweenness: number of shortest paths that pass through each edge.
 * Used by Girvan-Newman to identify community boundaries. */
void ifn_edge_betweenness(const IFN_CausalGraph* g, double* edge_between);

/* Community mutual information: I(communities; information_flow_pattern).
 * Measures how well community structure predicts information flow. */
double ifn_community_mi(const IFN_CausalGraph* g,
    const int* community_labels, int n_communities,
    const double* flow_pattern, int n_edges);

/* Normalized mutual information between two community partitions.
 * Standard metric for comparing community detection results. */
double ifn_nmi(const int* labels_a, const int* labels_b, int n);

#endif
