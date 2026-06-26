#ifndef IFN_CAUSAL_H
#define IFN_CAUSAL_H
#include "ifn_core.h"

/* Causal graphs and DAG operations. */

IFN_CausalGraph* ifn_causal_create(int n_nodes);
void ifn_causal_free(IFN_CausalGraph* g);
void ifn_causal_add_edge(IFN_CausalGraph* g, int from, int to, double weight);
void ifn_causal_remove_edge(IFN_CausalGraph* g, int from, int to);
bool ifn_causal_has_edge(const IFN_CausalGraph* g, int from, int to);

/* DAG operations */
bool ifn_is_dag(const IFN_CausalGraph* g);
void ifn_topological_sort(const IFN_CausalGraph* g, int* order);
void ifn_transitive_closure(const IFN_CausalGraph* g, int* closure);

/* d-separation: are X and Y d-separated by Z? */
bool ifn_d_separated(const IFN_CausalGraph* g, int x, int y, const int* z, int n_z);

/* Moral graph (connect parents, remove directions) */
void ifn_moralize(const IFN_CausalGraph* g, int* moral);

/* PC algorithm (constraint-based causal discovery) */
IFN_CausalGraph* ifn_pc_algorithm(const IFN_TimeSeries* data, int n_vars,
    int n_bins, double alpha);

/* Information flow on causal graph */
double ifn_causal_flow(const IFN_CausalGraph* g, int source, int target);
void ifn_causal_path_effect(const IFN_CausalGraph* g, int source, int target,
    const int* path, int path_len, double* direct, double* indirect);

/* Do-calculus interventions */
void ifn_do_intervention(IFN_CausalGraph* g, int node, double value);
double ifn_causal_effect(const IFN_CausalGraph* g, int cause, int effect,
    const double* data, int n);

#endif
