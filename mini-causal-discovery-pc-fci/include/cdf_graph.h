#ifndef CDF_GRAPH_H
#define CDF_GRAPH_H

/*
 * cdf_graph.h — Graph Operations for Causal Discovery
 *
 * Implements fundamental graph-theoretic operations needed by
 * the PC and FCI algorithms:
 *   - d-separation / m-separation
 *   - Skeleton construction via adjacency search
 *   - V-structure (collider) detection
 *   - Path finding and reachability
 *   - Ancestral relationships
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "cdf_core.h"

/* ── Path and Reachability ───────────────────────────────────────── */

/** Check if there is an undirected path from u to v in the graph.
 *  Uses BFS, ignoring edge directions. */
bool cdf_graph_has_undirected_path(const CdfGraph *g, int u, int v);

/** Check if there is a directed path from u to v.
 *  Uses BFS along directed edges only. */
bool cdf_graph_has_directed_path(const CdfGraph *g, int u, int v);

/** Check if v is an ancestor of u (including itself).
 *  v is an ancestor of u if there is a directed path v → ... → u. */
bool cdf_graph_is_ancestor(const CdfGraph *g, int v, int u);

/** Check if v is a possible ancestor of u in a PAG.
 *  Includes partially directed paths (∘→, →). */
bool cdf_graph_is_possible_ancestor(const CdfGraph *g, int v, int u);

/** Get the adjacency set (neighbors) of node v.
 *  Returns all nodes connected to v by any edge type. */
void cdf_graph_neighbors(const CdfGraph *g, int v, int *neighbors, int *n_nei);

/** Get parents of node v (nodes u with u → v). */
void cdf_graph_parents(const CdfGraph *g, int v, int *parents, int *n_par);

/** Get children of node v (nodes u with v → u). */
void cdf_graph_children(const CdfGraph *g, int v, int *children, int *n_chi);

/** Get the set Adj(v) \{u} — neighbors excluding a specific node. */
void cdf_graph_adj_except(const CdfGraph *g, int v, int exclude,
                           int *result, int *n_res);

/* ── d-Separation ────────────────────────────────────────────────── */

/**
 * Check if X and Y are d-separated by set Z in a DAG.
 *
 * d-separation criterion (Pearl, 1988):
 *   X ⊥_d Y | Z iff Z blocks every path between X and Y.
 *   A path is blocked by Z if:
 *     - It contains a chain i→m→j or fork i←m→j with m ∈ Z, OR
 *     - It contains a collider i→m←j where m ∉ Z and no descendant
 *       of m is in Z.
 *
 * Implementation: moralized graph BFS.
 *
 * @param g  causal DAG
 * @param x  node X
 * @param y  node Y
 * @param Z  conditioning set
 * @param nZ size of Z
 */
bool cdf_graph_d_separated(const CdfGraph *g, int x, int y,
                            const int *Z, int nZ);

/**
 * Check m-separation in a MAG (for FCI).
 * Similar to d-separation but handles bidirected edges (↔).
 * A path with ↔ u ↔ is treated as inducing dependence.
 */
bool cdf_graph_m_separated(const CdfGraph *g, int x, int y,
                            const int *Z, int nZ);

/* ── Skeleton ────────────────────────────────────────────────────── */

/**
 * Learn the skeleton (undirected graph) from data using the
 * adjacency search phase of PC.
 *
 * Starting with a complete undirected graph, iteratively remove
 * edges (u,v) if ∃ Z ⊆ Adj(u)\{v} (|Z| = ℓ) such that u ⊥ v | Z.
 *
 * This function populates g->edges with the skeleton and sets
 * g->sepsets for each removed edge.
 *
 * @param g     graph (will be modified with skeleton)
 * @param ds    dataset
 * @param config PC configuration
 * @return      number of CI tests performed
 */
int cdf_graph_skeleton_pc(CdfGraph *g, const CdfDataset *ds,
                           const CdfPCConfig *config);

/* ── V-Structures (Colliders) ────────────────────────────────────── */

/**
 * Detect v-structures (unshielded colliders) in the skeleton.
 *
 * A triple (u, v, w) forms a v-structure u → v ← w if:
 *   - u—v—w is in the skeleton (adjacent pairs)
 *   - u and w are NOT adjacent
 *   - v ∉ SepSet(u,w) (i.e., v is not in the separating set)
 *
 * Orients edges as u → v ← w.
 *
 * @param g     skeleton (modified to oriented CPDAG)
 * @return      number of v-structures found
 */
int cdf_graph_find_vstructures(CdfGraph *g);

/* ── Triangle Completion ─────────────────────────────────────────── */

/** Check if three nodes form a triangle (all pairwise adjacent). */
bool cdf_graph_is_triangle(const CdfGraph *g, int a, int b, int c);

/** Check if edge u—v is "shielded" (there exists w adjacent to both). */
bool cdf_graph_is_shielded(const CdfGraph *g, int u, int v);

/* ── Possible-D-SEP (FCI) ────────────────────────────────────────── */

/**
 * Find Possible-D-SEP set for a node in FCI.
 *
 * Possible-D-SEP(v) = {w | there is a possibly directed path from v to w
 *   in the current PAG, OR v and w are in the same discriminating path}.
 *
 * This is used in the FCI adjacency phase to find additional
 * conditioning sets that may not be in the initial adjacency set.
 *
 * @param g          current PAG
 * @param v          target node
 * @param max_len    max path length
 * @param pdsep      [p] output: 1 if node is in PDS(v)
 */
void cdf_graph_possible_d_sep(const CdfGraph *g, int v, int max_len,
                               int *pdsep);

/* ── Print ───────────────────────────────────────────────────────── */

/** Print the graph adjacency matrix. */
void cdf_graph_print(const CdfGraph *g);

/** Print a single edge type as a string. */
const char* cdf_graph_edge_char(CdfEdgeType type);

/** Print graph statistics: nodes, edges, edge type counts. */
void cdf_graph_stats(const CdfGraph *g);

#endif /* CDF_GRAPH_H */