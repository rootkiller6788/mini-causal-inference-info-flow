#ifndef CDF_FCI_H
#define CDF_FCI_H

/*
 * cdf_fci.h — FCI Algorithm (Fast Causal Inference)
 *
 * The FCI algorithm (Spirtes et al., 1999; Zhang, 2008) extends PC
 * to handle latent confounders and selection bias. It outputs a
 * Partial Ancestral Graph (PAG) representing the equivalence class
 * of Maximal Ancestral Graphs (MAGs).
 *
 * Key differences from PC:
 *   1. Does NOT assume causal sufficiency (allows latent variables)
 *   2. Uses Possible-D-SEP for additional conditioning sets
 *   3. Produces PAGs with ∘→, ∘—∘, ↔ edge types
 *   4. Applies 10 orientation rules (R1-R10) instead of 4
 *
 * FCI algorithm phases (Zhang, 2008):
 *   Phase 0: Init complete undirected graph
 *   Phase 1: PC-like adjacency search → skeleton
 *   Phase 2: V-structure orientation
 *   Phase 3: Possible-D-SEP search + re-execute CI tests
 *   Phase 4: Re-apply adjacency search using PDS sets
 *   Phase 5: Apply FCI orientation rules (R0-R10)
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "cdf_core.h"
#include "cdf_graph.h"

/* ── FCI Result ──────────────────────────────────────────────────── */

typedef struct {
    CdfGraph  *graph;            /* output PAG */
    int        n_ci_tests;       /* total CI tests */
    int        n_edges_removed;  /* edges removed */
    int        n_vstructs;       /* v-structures found */
    int        n_oriented;       /* edges oriented by rules */
    int        n_pds_tests;      /* CI tests in PDS phase */
    double     elapsed_sec;      /* runtime */
} CdfFCIResult;

/* ── FCI Algorithm ───────────────────────────────────────────────── */

/**
 * Run the full FCI algorithm.
 *
 * @param ds      dataset
 * @param config  FCI configuration (NULL for defaults)
 * @return        FCI result with PAG (caller must cdf_fci_result_free)
 */
CdfFCIResult* cdf_fci_run(const CdfDataset *ds, const CdfFCIConfig *config);

/** Free FCI result. */
void cdf_fci_result_free(CdfFCIResult *res);

/**
 * FCI Phase 1: Initial adjacency search (like PC skeleton).
 * Uses adjacency sets as initial conditioning candidates.
 *
 * @param g       complete graph
 * @param ds      dataset
 * @param config  FCI config
 * @return        CI tests performed
 */
int cdf_fci_initial_skeleton(CdfGraph *g, const CdfDataset *ds,
                              const CdfFCIConfig *config);

/**
 * FCI Phase 2: V-structure detection + initial orientation.
 * Same as PC but produces ∘→ markings for non-adjacent pairs
 * that need further investigation.
 */
int cdf_fci_vstructures(CdfGraph *g);

/**
 * FCI Phase 3: Possible-D-SEP computation and CI testing.
 * For each node v, find PDS(v) and test CI with conditioning
 * sets from PDS(v). This can remove additional edges.
 *
 * @param g       current PAG
 * @param ds      dataset
 * @param config  FCI config
 * @return        CI tests performed
 */
int cdf_fci_pds_phase(CdfGraph *g, const CdfDataset *ds,
                       const CdfFCIConfig *config);

/**
 * FCI Phase 4: Apply full FCI orientation rules (R0-R10).
 * Includes the expanded rule set for latent variable structures.
 */
int cdf_fci_orientation_phase(CdfGraph *g);

/* ── FCI-Specific Operations ─────────────────────────────────────── */

/**
 * RFCI (Really Fast Causal Inference): a simplified FCI variant
 * that avoids the expensive PDS phase by using only adjacency sets.
 * Less accurate but much faster for high-dimensional data.
 */
CdfFCIResult* cdf_fci_rfci_run(const CdfDataset *ds,
                                const CdfFCIConfig *config);

/** Convert a PAG to its implied MAG skeleton for visualization. */
void cdf_fci_pag_to_mag_skeleton(const CdfGraph *pag, CdfGraph *mag);

/** Check if a PAG is valid (no directed cycles excluding ↔ edges). */
bool cdf_fci_pag_is_valid(const CdfGraph *g);

/** Print FCI result. */
void cdf_fci_print_result(const CdfFCIResult *res);

/** Compare two PAGs: count edge type differences. */
int cdf_fci_compare_pags(const CdfGraph *g1, const CdfGraph *g2);

#endif /* CDF_FCI_H */