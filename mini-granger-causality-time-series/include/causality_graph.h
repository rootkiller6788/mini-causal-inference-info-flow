#ifndef CAUSALITY_GRAPH_H
#define CAUSALITY_GRAPH_H
#include "ts_core.h"
#include <stdbool.h>

/* Causality Graph from Pairwise Granger Tests
 *
 * Definition: Given M time series, construct a directed graph G=(V,E)
 * where edge i→j exists iff variable i Granger-causes variable j
 * at significance level α. The edge weight is -log10(p-value).
 *
 * This implements the "Granger causality graph" approach analogous to
 * Eichler (2013), "Causal inference with multiple time series."
 *
 * Reference: Eichler, M. (2013). Phil. Trans. R. Soc. A, 371(1997).
 */

typedef struct {
    int n_vars;
    int n_lags;
    bool** adjacency;       /* [n_vars][n_vars] adjacency matrix */
    double** edge_weights;  /* [n_vars][n_vars] -log10(p-value) as weight */
    double** f_statistics;  /* [n_vars][n_vars] F-statistics */
    double** p_values;      /* [n_vars][n_vars] p-values */
    int total_edges;
    double significance_level;
} CausalityGraph;

/* Construct causality graph from M-variate time series.
 * Performs M*(M-1) pairwise Granger causality tests.
 * Complexity: O(M² * max_lag * length * dim³). */
CausalityGraph* graph_build(const TimeSeries* ts, int max_lag, double alpha);
void graph_free(CausalityGraph* g);
void graph_print(const CausalityGraph* g);

/* Compute graph-theoretic metrics on the causality graph:
 * in_degree[i] = number of incoming edges to node i
 * out_degree[i] = number of outgoing edges from node i
 * betweenness[i] = normalized betweenness centrality (Brandes 2001)
 * Complexity: O(V * E) for betweenness. */
typedef struct {
    int* in_degree;
    int* out_degree;
    double* betweenness;
    int n_vars;
} GraphMetrics;

GraphMetrics* graph_metrics(const CausalityGraph* g);
void gm_free(GraphMetrics* m);
void gm_print(const GraphMetrics* m);

/* Detect causal hubs: nodes with out_degree > τ * (n_vars-1).
 * Returns array of hub node indices, terminated by -1.
 * τ is the fraction threshold (typically 0.5). */
int* graph_find_hubs(const CausalityGraph* g, double tau_fraction);

/* Check if the causality graph is a DAG (no cycles).
 * Uses Kahn's algorithm for topological sorting. O(V+E). */
bool graph_is_dag(const CausalityGraph* g);

#endif
