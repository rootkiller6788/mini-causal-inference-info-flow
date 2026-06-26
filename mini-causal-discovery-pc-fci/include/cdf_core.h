#ifndef CDF_CORE_H
#define CDF_CORE_H

/*
 * cdf_core.h — Core data structures for causal discovery (PC & FCI)
 *
 * Causal discovery aims to learn causal graphs from observational data.
 * The PC algorithm (Spirtes-Glymour-Scheines, 1991) assumes causal
 * sufficiency (no latent confounders) and outputs a CPDAG.
 * The FCI algorithm extends PC for latent variables, outputting a PAG.
 *
 * Core concepts:
 *   DAG:    Directed Acyclic Graph (true causal structure)
 *   CPDAG:  Completed Partially Directed Acyclic Graph (PC output)
 *   MAG:    Maximal Ancestral Graph (with latents)
 *   PAG:    Partial Ancestral Graph (FCI output)
 *   d-sep:  d-separation criterion for reading CI from DAGs
 *   SepSet: separating set — variables that make two nodes CI
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ── Constants ───────────────────────────────────────────────────── */

#define CDF_MAX_NODES     64    /* max number of variables */
#define CDF_MAX_DEGREE    32    /* max neighbors per node */
#define CDF_MAX_SAMPLES   10000 /* max dataset size */
#define CDF_PI            3.14159265358979323846
#define CDF_EPS           1e-12
#define CDF_ALPHA_DEFAULT 0.05  /* default significance level */

/* ── Edge Types ──────────────────────────────────────────────────── */

/** Possible edge states in CPDAG / PAG */
typedef enum {
    CDF_EDGE_NONE       = 0,  /* no edge */
    CDF_EDGE_UNDIRECTED = 1,  /* i —— j  (in skeleton/CPDAG) */
    CDF_EDGE_DIRECTED   = 2,  /* i →  j */
    CDF_EDGE_BIDIRECTED = 3,  /* i ↔ j  (in PAG, indicates latent confounder) */
    CDF_EDGE_PARTIAL_I  = 4,  /* i ∘→ j (partially directed, tail uncertain) */
    CDF_EDGE_PARTIAL_J  = 5,  /* i ←∘ j */
    CDF_EDGE_NONDIR     = 6,  /* i ∘—∘ j (nondirected, both ends uncertain) */
    CDF_EDGE_UNKNOWN    = 7   /* unclassified */
} CdfEdgeType;

/* ── Edge ────────────────────────────────────────────────────────── */

/** Single edge in the causal graph */
typedef struct {
    int         u;         /* source node index */
    int         v;         /* target node index */
    CdfEdgeType type;      /* edge type */
} CdfEdge;

/* ── SepSet (Separating Set) ─────────────────────────────────────── */

/** Records the set of variables that d-separate / CI-ize a pair */
typedef struct {
    int     u;             /* node u */
    int     v;             /* node v */
    int     *set;          /* [size] array of separating variables */
    int      size;         /* size of separating set */
    double   p_value;      /* p-value of the maximal CI test */
} CdfSepSet;

/* ── Dataset ─────────────────────────────────────────────────────── */

/**
 * CdfDataset — N samples × p variables, column-major storage.
 *
 * data[i * N + j] = value of variable i at sample j.
 */
typedef struct {
    double *data;          /* [p * N] column-major: var-i at sample-j = data[i*N + j] */
    char   **var_names;   /* [p] variable names */
    int     N;            /* number of samples */
    int     p;            /* number of variables */
} CdfDataset;

/* ── Adjacency List Representation ───────────────────────────────── */

/** Neighbor list for one node */
typedef struct {
    int  neighbors[CDF_MAX_DEGREE];
    int  n_neighbors;
} CdfNeighborList;

/* ── Causal Graph ────────────────────────────────────────────────── */

/**
 * CdfGraph — causal graph representation.
 *
 * Uses adjacency matrix: edges[p * u + v] = type of edge u→v.
 * For undirected edges, both directions must be symmetric.
 * For directed edge u→v: edges[u*p+v] = DIRECTED, edges[v*p+u] = NONE.
 */
typedef struct {
    CdfEdgeType *edges;    /* [p * p] adjacency matrix, row-major */
    CdfSepSet   *sepsets;  /* [p * p] separation sets for non-adjacent pairs */
    int         *sepset_count; /* [p * p] has sepset entry if sepset_count[idx] >= 0 */
    int          p;        /* number of nodes */
    bool         is_cpdag; /* true if this is a CPDAG (PC output) */
    bool         is_pag;   /* true if this is a PAG (FCI output) */
} CdfGraph;

/* ── PC Algorithm Configuration ──────────────────────────────────── */

typedef struct {
    double alpha;           /* significance level for CI tests */
    int    max_cond_size;   /* maximum conditioning set size (-1 = unlimited) */
    bool   use_fisher_z;    /* use Fisher's z-test (true) or mutual information (false) */
    bool   conservative;    /* use conservative orientation (Ramsey et al.) */
    bool   verbose;         /* print progress */
} CdfPCConfig;

/* ── FCI Algorithm Configuration ─────────────────────────────────── */

typedef struct {
    double alpha;           /* significance level */
    int    max_cond_size;   /* max conditioning set size */
    int    max_path_length; /* max path length for Possible-D-SEP search */
    bool   use_fisher_z;
    bool   conservative;
    bool   complete_rule_set; /* apply all 10 FCI orientation rules */
    bool   verbose;
} CdfFCIConfig;

/* ── Lifecycle ───────────────────────────────────────────────────── */

/** Create a dataset from raw data array. */
CdfDataset* cdf_dataset_create(const double *data, int N, int p);

/** Free a dataset. */
void cdf_dataset_free(CdfDataset *ds);

/** Create an empty causal graph with p nodes. */
CdfGraph* cdf_graph_create(int p);

/** Free a causal graph. */
void cdf_graph_free(CdfGraph *g);

/** Add an edge to the graph. */
void cdf_graph_add_edge(CdfGraph *g, int u, int v, CdfEdgeType type);

/** Remove an edge from the graph. */
void cdf_graph_remove_edge(CdfGraph *g, int u, int v);

/** Check if an edge exists between u and v (any type). */
bool cdf_graph_has_edge(const CdfGraph *g, int u, int v);

/** Get the edge type from u to v. */
CdfEdgeType cdf_graph_edge_type(const CdfGraph *g, int u, int v);

/** Check if there is a directed path from u to v (reachability). */
bool cdf_graph_reachable(const CdfGraph *g, int u, int v);

/** Create default PC configuration. */
CdfPCConfig cdf_pc_config_default(void);

/** Create default FCI configuration. */
CdfFCIConfig cdf_fci_config_default(void);

/** Get variable names from dataset (NULL if not set). */
const char* cdf_dataset_var_name(const CdfDataset *ds, int i);

/** Count export functions. */
int cdf_version(void);

#endif /* CDF_CORE_H */