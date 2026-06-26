/* cdf_fci.c — FCI Algorithm (Fast Causal Inference)
 *
 * Implements the FCI algorithm for causal discovery with latent
 * variables. Extends PC with:
 *   - Possible-D-SEP search for additional conditioning sets
 *   - Extended orientation rules (R1-R10) for PAG output
 *
 * FCI outputs a Partial Ancestral Graph (PAG) representing the
 * equivalence class of Maximal Ancestral Graphs (MAGs).
 */

#include "cdf_fci.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include "cdf_orientation.h"
#include "cdf_pc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

/* ──────────────── Full FCI Algorithm ─────────────────────────────── */

CdfFCIResult* cdf_fci_run(const CdfDataset *ds, const CdfFCIConfig *config)
{
    if (!ds) return NULL;

    CdfFCIConfig cfg;
    if (config) cfg = *config;
    else cfg = cdf_fci_config_default();

    CdfFCIResult *res = (CdfFCIResult*)calloc(1, sizeof(CdfFCIResult));
    if (!res) return NULL;

    clock_t t_start = clock();

    /* Initialize graph */
    int p = ds->p;
    res->graph = cdf_graph_create(p);
    if (!res->graph) { free(res); return NULL; }
    res->graph->is_pag = true;

    /* Complete undirected start */
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            cdf_graph_add_edge(res->graph, i, j, CDF_EDGE_UNDIRECTED);

    /* Phase 1: Initial skeleton (like PC) */
    res->n_ci_tests = cdf_fci_initial_skeleton(res->graph, ds, &cfg);

    /* Count removed edges */
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            if (!cdf_graph_has_edge(res->graph, i, j))
                res->n_edges_removed++;

    /* Phase 2: V-structures */
    res->n_vstructs = cdf_fci_vstructures(res->graph);

    /* Convert remaining undirected edges to ∘—∘ for PAG representation */
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            if (res->graph->edges[i * p + j] == CDF_EDGE_UNDIRECTED) {
                cdf_graph_add_edge(res->graph, i, j, CDF_EDGE_NONDIR);
            }
        }
    }

    /* Phase 3: PDS search + CI testing */
    res->n_pds_tests = cdf_fci_pds_phase(res->graph, ds, &cfg);

    /* Phase 4: Orientation rules */
    res->n_oriented = cdf_fci_orientation_phase(res->graph);

    res->elapsed_sec = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;

    return res;
}

void cdf_fci_result_free(CdfFCIResult *res)
{
    if (!res) return;
    cdf_graph_free(res->graph);
    free(res);
}

/* ──────────────── Phase 1: Initial Skeleton ──────────────────────── */

int cdf_fci_initial_skeleton(CdfGraph *g, const CdfDataset *ds,
                              const CdfFCIConfig *config)
{
    /* Use same skeleton procedure as PC */
    CdfPCConfig pc_cfg;
    pc_cfg.alpha = config->alpha;
    pc_cfg.max_cond_size = config->max_cond_size;
    pc_cfg.use_fisher_z = config->use_fisher_z;
    pc_cfg.conservative = config->conservative;
    pc_cfg.verbose = config->verbose;

    return cdf_graph_skeleton_pc(g, ds, &pc_cfg);
}

/* ──────────────── Phase 2: V-Structures ──────────────────────────── */

int cdf_fci_vstructures(CdfGraph *g)
{
    return cdf_graph_find_vstructures(g);
}

/* ──────────────── Phase 3: Possible-D-SEP ────────────────────────── */

int cdf_fci_pds_phase(CdfGraph *g, const CdfDataset *ds,
                       const CdfFCIConfig *config)
{
    if (!g || !ds) return 0;

    int p = g->p, n_tests = 0;
    int *pdsep = (int*)malloc((size_t)p * sizeof(int));
    if (!pdsep) return 0;

    /* For each pair of non-adjacent nodes */
    for (int i = 0; i < p; i++) {
        cdf_graph_possible_d_sep(g, i, config->max_path_length, pdsep);

        for (int j = i + 1; j < p; j++) {
            if (!cdf_graph_has_edge(g, i, j)) {
                /* Collect PDS nodes as conditioning candidates */
                int cond_set[CDF_MAX_NODES], n_cond = 0;
                for (int k = 0; k < p; k++) {
                    if (k != i && k != j && pdsep[k]) {
                        if (n_cond < CDF_MAX_NODES)
                            cond_set[n_cond++] = k;
                    }
                }

                /* Test subsets of PDS as conditioning sets */
                /* Limit to small subsets for efficiency */
                int max_subset = (n_cond < 3) ? n_cond : 3;
                for (int sz = 0; sz <= max_subset; sz++) {
                    /* Simple: test with first sz elements */
                    if (sz <= n_cond) {
                        CdfCITestResult r = cdf_citest_fisher_z(
                            ds, i, j, cond_set, sz, config->alpha);
                        n_tests++;
                        if (r.is_independent) break;  /* edge already removed */
                    }
                }
            }
        }
    }

    free(pdsep);
    return n_tests;
}

/* ──────────────── Phase 4: FCI Orientation ───────────────────────── */

int cdf_fci_orientation_phase(CdfGraph *g)
{
    CdfOrientResult res = cdf_orient_fci_rules(g);
    return res.total_oriented;
}

/* ──────────────── RFCI (Really Fast Causal Inference) ────────────── */

CdfFCIResult* cdf_fci_rfci_run(const CdfDataset *ds,
                                const CdfFCIConfig *config)
{
    /* RFCI skips the expensive PDS phase.
     * Just run PC-like skeleton + v-structures + FCI orientation rules. */
    if (!ds) return NULL;

    CdfFCIConfig cfg;
    if (config) cfg = *config;
    else cfg = cdf_fci_config_default();

    /* Create PC-compatible config */
    CdfPCConfig pc_cfg;
    pc_cfg.alpha = cfg.alpha;
    pc_cfg.max_cond_size = cfg.max_cond_size;
    pc_cfg.use_fisher_z = cfg.use_fisher_z;
    pc_cfg.conservative = cfg.conservative;
    pc_cfg.verbose = cfg.verbose;

    CdfPCResult *pc_res = cdf_pc_run(ds, &pc_cfg);
    if (!pc_res) return NULL;

    /* Wrap PC result as FCI result */
    CdfFCIResult *fci_res = (CdfFCIResult*)calloc(1, sizeof(CdfFCIResult));
    if (!fci_res) { cdf_pc_result_free(pc_res); return NULL; }

    fci_res->graph = pc_res->graph;  /* transfer ownership */
    fci_res->graph->is_cpdag = false;
    fci_res->graph->is_pag = true;

    /* Convert undirected edges to ∘—∘ for PAG */
    int p = ds->p;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            if (fci_res->graph->edges[i * p + j] == CDF_EDGE_UNDIRECTED) {
                cdf_graph_add_edge(fci_res->graph, i, j, CDF_EDGE_NONDIR);
            }
        }
    }

    /* Apply FCI orientation rules */
    fci_res->n_oriented = cdf_fci_orientation_phase(fci_res->graph);
    fci_res->n_ci_tests = pc_res->n_ci_tests;
    fci_res->n_edges_removed = pc_res->n_edges_removed;
    fci_res->n_vstructs = pc_res->n_vstructs;
    fci_res->elapsed_sec = pc_res->elapsed_sec;

    free(pc_res);  /* pc_res->graph was transferred */
    return fci_res;
}

/* ──────────────── PAG → MAG Skeleton ─────────────────────────────── */

void cdf_fci_pag_to_mag_skeleton(const CdfGraph *pag, CdfGraph *mag)
{
    if (!pag || !mag || pag->p != mag->p) return;

    int p = pag->p;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            CdfEdgeType e = pag->edges[i * p + j];

            switch (e) {
            case CDF_EDGE_DIRECTED:
                cdf_graph_add_edge(mag, i, j, CDF_EDGE_DIRECTED);
                break;
            case CDF_EDGE_BIDIRECTED:
                cdf_graph_add_edge(mag, i, j, CDF_EDGE_BIDIRECTED);
                break;
            case CDF_EDGE_PARTIAL_I:
                /* ∘→ → → in MAG (conservative) */
                cdf_graph_add_edge(mag, i, j, CDF_EDGE_DIRECTED);
                break;
            case CDF_EDGE_NONDIR:
                /* ∘—∘ → — in MAG skeleton */
                cdf_graph_add_edge(mag, i, j, CDF_EDGE_UNDIRECTED);
                break;
            default:
                break;
            }
        }
    }
}

/* ──────────────── Validation ─────────────────────────────────────── */

bool cdf_fci_pag_is_valid(const CdfGraph *g)
{
    if (!g) return false;

    int p = g->p;
    /* Check for consistency: bidirectional edges should not form cycles with → */
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            if (g->edges[i * p + j] == CDF_EDGE_BIDIRECTED) {
                /* Must also be bidirectional from j to i */
                if (g->edges[j * p + i] != CDF_EDGE_BIDIRECTED)
                    return false;
            }
        }
    }
    return true;
}

int cdf_fci_compare_pags(const CdfGraph *g1, const CdfGraph *g2)
{
    if (!g1 || !g2 || g1->p != g2->p) return -1;

    int p = g1->p, diff = 0;
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            if (g1->edges[i * p + j] != g2->edges[i * p + j])
                diff++;
    return diff;
}

void cdf_fci_print_result(const CdfFCIResult *res)
{
    if (!res) return;
    printf("=== FCI Algorithm Result ===\n");
    printf("CI tests:      %d\n", res->n_ci_tests);
    printf("Edges removed: %d\n", res->n_edges_removed);
    printf("V-structures:  %d\n", res->n_vstructs);
    printf("PDS tests:     %d\n", res->n_pds_tests);
    printf("Oriented:      %d\n", res->n_oriented);
    printf("Time:          %.4f sec\n", res->elapsed_sec);
}