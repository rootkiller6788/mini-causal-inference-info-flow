#ifndef DAG_REPRESENTATION_H
#define DAG_REPRESENTATION_H

/*
 * dag_representation.h — DAG Data Structures and Graph Algorithms
 *
 * A Directed Acyclic Graph (DAG) is the fundamental structure in causal
 * inference. Nodes represent variables; directed edges represent direct
 * causal relationships (X -> Y means X directly causes Y).
 *
 * Key operations:
 *   - Topological ordering (causal ordering)
 *   - Cycle detection (DAG validation)
 *   - Path finding (directed paths, back-door paths)
 *   - Ancestor/descendant relations
 *   - d-separation testing
 *   - Moral graph / undirected skeleton
 *
 * References:
 *   - Pearl, "Causality: Models, Reasoning, and Inference", 2009
 *   - Spirtes, Glymour, Scheines, "Causation, Prediction, and Search", 2000
 *   - Koller & Friedman, "Probabilistic Graphical Models", 2009
 */

#include <stdbool.h>
#include <stddef.h>

/* ---- DAG Structure ---- */
/* Adjacency matrix representation: adj[i][j] = 1 if i -> j (row-major).
 * Also stores parent/child lists for efficient traversal. */
typedef struct {
    int     n;          /* number of nodes (variables)               */
    int    *adj;        /* adjacency matrix, n x n, row-major        */
    int    *parents;    /* list of parent indices per node           */
    int    *children;   /* list of child indices per node            */
    int    *n_parents;  /* number of parents for each node           */
    int    *n_children; /* number of children for each node          */
    int   *topo_order;  /* topological ordering (causal order)       */
    char  **names;      /* variable names (optional, max 32 chars)   */
} DAG;

/* ---- Path Structure ---- */
typedef struct {
    int    *nodes;      /* sequence of node indices in path          */
    int     length;     /* number of nodes in path                   */
    bool    is_directed;/* true = directed path, false = any path    */
} Path;

/* ---- d-Separation Result ---- */
typedef struct {
    bool    is_separated; /* whether X and Y are d-separated by Z    */
    Path   *blocking_path;/* a path that shows non-separation (if any)*/
} DSepResult;

/* ---- DAG Lifecycle ---- */
DAG* dag_create(int n);
void dag_free(DAG* dag);

/* Add a directed edge i -> j. Returns 0 on success, -1 if edge creates cycle. */
int dag_add_edge(DAG* dag, int i, int j);

/* Remove a directed edge i -> j. */
int dag_remove_edge(DAG* dag, int i, int j);

/* Check if there is an edge i -> j. */
bool dag_has_edge(const DAG* dag, int i, int j);

/* ---- Structural Properties ---- */

/* Compute topological ordering (Kahn's algorithm).
 * Returns 0 on success, -1 if graph has a cycle (not a DAG). */
int dag_topological_sort(DAG* dag);

/* Check if the graph is acyclic. */
bool dag_is_acyclic(const DAG* dag);

/* Check if there is a directed path from i to j (BFS/DFS). */
bool dag_has_path(const DAG* dag, int from, int to);

/* Find all ancestors of a node (nodes with directed paths TO the node). */
int dag_ancestors(const DAG* dag, int node, int *ancestors, int max_n);

/* Find all descendants of a node (nodes with directed paths FROM the node). */
int dag_descendants(const DAG* dag, int node, int *descendants, int max_n);

/* ---- Path Operations ---- */

/* Find a directed path from  to  if one exists. */
Path* dag_find_directed_path(const DAG* dag, int from, int to);

/* Find all back-door paths between  and .
 * A back-door path is any path that starts with an edge INTO . */
int dag_find_backdoor_paths(const DAG* dag, int from, int to,
                              Path **paths, int max_paths);

/* ---- d-Separation ---- */

/* Test if X and Y are d-separated by conditioning set Z.
 * d-separation is the graphical criterion for conditional independence. */
DSepResult* dag_d_separation(const DAG* dag, int x, int y,
                              const int *z, int n_z);

/* ---- Graph Transformations ---- */

/* Create the moral graph (connect all parents of each node, then undirect). */
DAG* dag_moral_graph(const DAG* dag);

/* Remove outgoing edges from a node (graph surgery for intervention). */
DAG* dag_do_intervention(const DAG* dag, int node);

/* ---- Constants ---- */
#define DAG_MAX_NODES    200
#define DAG_MAX_PATHS     50
#define DAG_MAX_PATH_LEN  50
#define DAG_NAME_LEN      32

#endif /* DAG_REPRESENTATION_H */
