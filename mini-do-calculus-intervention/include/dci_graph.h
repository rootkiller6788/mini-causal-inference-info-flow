#ifndef DCI_GRAPH_H
#define DCI_GRAPH_H

#include "dci_core.h"

/* ---- Graph Algorithms for Causal Inference ---- */

/* Path in causal graph */
typedef struct {
    int nodes[DCI_MAX_VARS];
    int length;
    bool is_directed;
    bool is_backdoor;
} DCI_Path;

/* Topological order */
typedef struct {
    int order[DCI_MAX_VARS];
    int n_nodes;
} DCI_TopoOrder;

/* d-separation test */
typedef struct {
    bool is_d_separated;
    int* blocking_set;
    int n_blocking;
} DCI_DSepResult;

/* Path finding */
int dci_find_all_paths(DCI_Graph* g, int from, int to,
    DCI_Path* paths, int max_paths);
int dci_find_directed_paths(DCI_Graph* g, int from, int to,
    DCI_Path* paths, int max_paths);
int dci_find_backdoor_paths(DCI_Graph* g, int from, int to,
    DCI_Path* paths, int max_paths);

/* d-separation (Pearl, 1988) */
DCI_DSepResult dci_d_separation(DCI_Graph* g, int x, int y,
    const int* conditioning, int n_cond);
bool dci_is_d_separated(DCI_Graph* g, int x, int y,
    const int* z, int n_z);

/* Topological sorting */
DCI_TopoOrder dci_topological_sort(DCI_Graph* g);
bool dci_graph_is_dag(DCI_Graph* g);
int dci_graph_cycle_detect(DCI_Graph* g);

/* Graph surgery (for do-calculus) */
void dci_graph_remove_incoming(DCI_Graph* g, int node, DCI_Graph* result);
void dci_graph_remove_outgoing(DCI_Graph* g, int node, DCI_Graph* result);
void dci_graph_mutilate(DCI_Graph* g, const int* x_nodes, int n_x,
    DCI_Graph* result);

/* Conditional independence test (simulation-based) */
bool dci_test_independence(DCI_SCM* scm, int x, int y,
    const int* z, int n_z, int n_samples);

/* Causal ordering (Simon, 1953) */
int dci_causal_ordering(DCI_Graph* g, int* order);

/* Parents/children/ancestors/descendants */
void dci_get_parents(DCI_Graph* g, int node, int* parents, int* n);
void dci_get_children(DCI_Graph* g, int node, int* children, int* n);
void dci_get_ancestors(DCI_Graph* g, int node, int* ancestors, int* n);
void dci_get_descendants(DCI_Graph* g, int node, int* descendants, int* n);
bool dci_is_ancestor(DCI_Graph* g, int ancestor, int node);
bool dci_is_descendant(DCI_Graph* g, int descendant, int node);

void dci_transitive_closure(DCI_Graph* g, DCI_Graph* closure);
int dci_find_minimal_separator(DCI_Graph* g, int x, int y, int* separator);
void dci_moral_graph(DCI_Graph* g, DCI_Graph* moral);
int dci_find_v_structures(DCI_Graph* g, int (*triples)[3], int max_triples);
void dci_graph_skeleton(DCI_Graph* g, DCI_Graph* skeleton);
int dci_connected_components(DCI_Graph* g, int* component_id);
bool dci_collider_bias_check(DCI_Graph* g, int cause, int effect, int collider);

#endif
