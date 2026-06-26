#ifndef CDF_PC_H
#define CDF_PC_H

/*
 * cdf_pc.h — PC Algorithm (Peter-Clark)
 *
 * The PC algorithm (Spirtes & Glymour, 1991) learns the causal
 * structure from observational data under the assumptions of:
 *   1. Causal Markov Condition (CMC)
 *   2. Faithfulness
 *   3. Causal Sufficiency (no unobserved confounders)
 *
 * It outputs a Completed Partially Directed Acyclic Graph (CPDAG)
 * representing the Markov equivalence class of the true DAG.
 *
 * Algorithm phases:
 *   Phase 0: Initialize complete undirected graph
 *   Phase 1: Skeleton discovery — remove edges via CI tests
 *   Phase 2: V-structure orientation — detect colliders
 *   Phase 3: Meek's rules — propagate orientations (R1-R4)
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "cdf_core.h"
#include "cdf_graph.h"

/* ── PC Result ───────────────────────────────────────────────────── */

typedef struct {
    CdfGraph  *graph;          /* output CPDAG */
    int        n_ci_tests;     /* total CI tests performed */
    int        n_edges_removed;/* edges removed in skeleton phase */
    int        n_vstructs;     /* v-structures detected */
    int        n_oriented;     /* edges oriented by Meek rules */
    double     elapsed_sec;    /* runtime (approximate) */
} CdfPCResult;

/* ── PC Algorithm ────────────────────────────────────────────────── */

/**
 * Run the full PC algorithm on a dataset.
 *
 * @param ds      dataset (N samples × p variables)
 * @param config  PC configuration (NULL for defaults)
 * @return        PC result with CPDAG (caller must cdf_pc_result_free)
 */
CdfPCResult* cdf_pc_run(const CdfDataset *ds, const CdfPCConfig *config);

/** Free PC result. */
void cdf_pc_result_free(CdfPCResult *res);

/**
 * PC Phase 0: Initialize the complete undirected graph over p nodes.
 * Returns a fully connected skeleton.
 */
CdfGraph* cdf_pc_init_graph(int p);

/**
 * PC Phase 1: Skeleton adjacency search.
 * Removes edges (i,j) if ∃ S ⊆ Adj(i)\{j} such that i ⊥ j | S.
 * Stores SepSet(i,j) for each removed edge.
 *
 * @param g       complete graph (modified in-place)
 * @param ds      dataset
 * @param config  PC config
 * @return        number of CI tests performed
 */
int cdf_pc_skeleton_phase(CdfGraph *g, const CdfDataset *ds,
                           const CdfPCConfig *config);

/**
 * PC Phase 2: V-structure detection.
 * For each unshielded triple (i—k—j) where i,j not adjacent:
 *   If k ∉ SepSet(i,j), orient i → k ← j.
 *
 * @param g  skeleton (modified with v-structures)
 * @return   number of v-structures found
 */
int cdf_pc_vstructure_phase(CdfGraph *g);

/**
 * PC Phase 3: Apply Meek's orientation rules (R1-R4).
 * Exhaustively applies rules until no more orientations possible.
 *
 * @param g  graph with skeleton + v-structures (modified)
 * @return   number of edges oriented
 */
int cdf_pc_orientation_phase(CdfGraph *g);

/**
 * Conservative PC: only orient v-structures if the CI tests
 * are unambiguous (all conditioning sets of the same size
 * give the same independence judgment).
 *
 * @param g       skeleton
 * @param ds      dataset
 * @param config  PC config
 * @return        number of conservative v-structures
 */
int cdf_pc_conservative_vstructures(CdfGraph *g, const CdfDataset *ds,
                                     const CdfPCConfig *config);

/* ── PC Utilities ────────────────────────────────────────────────── */

/** Print PC algorithm result summary. */
void cdf_pc_print_result(const CdfPCResult *res);

/** Compute the structural Hamming distance (SHD) between two CPDAGs.
 *  Counts: extra edges + missing edges + orientation differences. */
int cdf_pc_shd(const CdfGraph *g1, const CdfGraph *g2);

/** Check if a CPDAG is valid (no directed cycles, edges consistent). */
bool cdf_pc_cpdag_is_valid(const CdfGraph *g);

/** Extract the undirected skeleton from a CPDAG. */
void cdf_pc_extract_skeleton(const CdfGraph *cpdag, CdfGraph *skeleton);

/** Count directed edges in the graph. */
int cdf_pc_count_directed(const CdfGraph *g);

/** Count undirected edges in the graph. */
int cdf_pc_count_undirected(const CdfGraph *g);

#endif /* CDF_PC_H */