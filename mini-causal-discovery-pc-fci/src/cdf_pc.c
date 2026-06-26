/* cdf_pc.c — PC Algorithm (Peter-Clark)
 *
 * Implements the complete PC algorithm for causal discovery:
 *   Phase 1: Skeleton learning via CI tests
 *   Phase 2: V-structure (collider) detection
 *   Phase 3: Meek's orientation rules (R1-R4)
 *
 * PC assumes: Causal Markov, Faithfulness, Sufficiency
 * Output: Completed Partially Directed Acyclic Graph (CPDAG)
 */

#include "cdf_pc.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include "cdf_orientation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

/* ──────────────── Full PC Algorithm ──────────────────────────────── */

CdfPCResult* cdf_pc_run(const CdfDataset *ds, const CdfPCConfig *config)
{
    if (!ds) return NULL;

    CdfPCConfig cfg;
    if (config) cfg = *config;
    else cfg = cdf_pc_config_default();

    CdfPCResult *res = (CdfPCResult*)calloc(1, sizeof(CdfPCResult));
    if (!res) return NULL;

    clock_t t_start = clock();

    /* Phase 0: Initialize complete graph */
    res->graph = cdf_pc_init_graph(ds->p);
    if (!res->graph) { free(res); return NULL; }
    res->graph->is_cpdag = true;

    /* Phase 1: Skeleton */
    res->n_ci_tests = cdf_pc_skeleton_phase(res->graph, ds, &cfg);
    res->n_edges_removed = 0;
    /* Count edges removed */
    int p = ds->p;
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            if (!cdf_graph_has_edge(res->graph, i, j))
                res->n_edges_removed++;

    /* Phase 2: V-structures */
    res->n_vstructs = cdf_pc_vstructure_phase(res->graph);

    /* Phase 3: Meek rules */
    res->n_oriented = cdf_pc_orientation_phase(res->graph);

    res->elapsed_sec = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;

    return res;
}

void cdf_pc_result_free(CdfPCResult *res)
{
    if (!res) return;
    cdf_graph_free(res->graph);
    free(res);
}

/* ──────────────── Phase 0: Initialize ────────────────────────────── */

CdfGraph* cdf_pc_init_graph(int p)
{
    CdfGraph *g = cdf_graph_create(p);
    if (!g) return NULL;

    /* Complete undirected graph */
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            cdf_graph_add_edge(g, i, j, CDF_EDGE_UNDIRECTED);

    g->is_cpdag = true;
    return g;
}

/* ──────────────── Phase 1: Skeleton ──────────────────────────────── */

int cdf_pc_skeleton_phase(CdfGraph *g, const CdfDataset *ds,
                           const CdfPCConfig *config)
{
    return cdf_graph_skeleton_pc(g, ds, config);
}

/* ──────────────── Phase 2: V-Structures ──────────────────────────── */

int cdf_pc_vstructure_phase(CdfGraph *g)
{
    return cdf_graph_find_vstructures(g);
}

/* ──────────────── Phase 3: Meek Rules ────────────────────────────── */

int cdf_pc_orientation_phase(CdfGraph *g)
{
    CdfOrientResult res = cdf_orient_pc_rules(g);
    return res.total_oriented;
}

/* ──────────────── Conservative PC ────────────────────────────────── */

int cdf_pc_conservative_vstructures(CdfGraph *g, const CdfDataset *ds,
                                     const CdfPCConfig *config)
{
    if (!g || !ds) return 0;

    int p = g->p, n_cons = 0;

    /* Conservative v-structure: only orient if all conditioning sets
     * of the same size give the same independence/sep-set judgement */

    for (int v = 0; v < p; v++) {
        int adj_v[CDF_MAX_DEGREE], n_adj_v;
        cdf_graph_neighbors(g, v, adj_v, &n_adj_v);

        for (int ui = 0; ui < n_adj_v; ui++) {
            int u = adj_v[ui];
            for (int wi = ui + 1; wi < n_adj_v; wi++) {
                int w = adj_v[wi];
                if (cdf_graph_has_edge(g, u, w)) continue;  /* shielded */

                /* Test all subsets of Adj(u)∩Adj(w) of appropriate size */
                int adj_u_w[CDF_MAX_DEGREE], n_uw;
                /* Intersection of adj(u) and adj(w) */
                int adj_u[CDF_MAX_DEGREE], n_adj_u;
                int adj_w[CDF_MAX_DEGREE], n_adj_w;
                cdf_graph_neighbors(g, u, adj_u, &n_adj_u);
                cdf_graph_neighbors(g, w, adj_w, &n_adj_w);

                n_uw = 0;
                for (int a = 0; a < n_adj_u; a++) {
                    for (int b = 0; b < n_adj_w; b++) {
                        if (adj_u[a] == adj_w[b] && adj_u[a] != v) {
                            adj_u_w[n_uw++] = adj_u[a];
                        }
                    }
                }

                /* Conservative check: test CI for all subsets of adj_u_w
                 * of size up to min(n_uw, config->max_cond_size).
                 * Only orient if all tests agree on independence.
                 * This implements Ramsey et al.'s conservative PC. */
                int max_test = config ? ((config->max_cond_size >= 0 &&
                    config->max_cond_size < n_uw) ? config->max_cond_size : n_uw) : 0;
                if (max_test < 0) max_test = 0;

                int all_indep = -1; /* unset; 0=dep, 1=indep */
                int ambiguous = 0;

                /* For small n_uw, brute-force test subsets of size up to max_test */
                /* Test empty set */
                {
                    CdfCITestResult r = cdf_citest_fisher_z(ds, u, w, NULL, 0,
                                                            config ? config->alpha : 0.05);
                    int result = r.is_independent ? 1 : 0;
                    if (all_indep == -1) all_indep = result;
                    else if (all_indep != result) { ambiguous = 1; }
                }

                /* Test singletons */
                for (int si = 0; si < n_uw && !ambiguous && max_test >= 1; si++) {
                    int ss[1] = {adj_u_w[si]};
                    CdfCITestResult r = cdf_citest_fisher_z(ds, u, w, ss, 1,
                                                            config ? config->alpha : 0.05);
                    int result = r.is_independent ? 1 : 0;
                    if (all_indep == -1) all_indep = result;
                    else if (all_indep != result) { ambiguous = 1; break; }
                }

                /* Test pairs */
                for (int si = 0; si < n_uw && !ambiguous && max_test >= 2; si++) {
                    for (int sj = si + 1; sj < n_uw && !ambiguous; sj++) {
                        int ss[2] = {adj_u_w[si], adj_u_w[sj]};
                        CdfCITestResult r = cdf_citest_fisher_z(ds, u, w, ss, 2,
                                                                config ? config->alpha : 0.05);
                        int result = r.is_independent ? 1 : 0;
                        if (all_indep == -1) all_indep = result;
                        else if (all_indep != result) { ambiguous = 1; break; }
                    }
                }

                /* Conservative decision:
                 * all_indep == 1: all tests say independent → v IS in SepSet → NOT collider
                 * all_indep == 0: all tests say dependent → v NOT in SepSet → IS collider
                 * ambiguous: mixed results → do not orient (conservative principle)
                 */
                if (!ambiguous && all_indep == 0) {
                    /* Reliable collider: orient u→v←w */
                    cdf_graph_add_edge(g, u, v, CDF_EDGE_DIRECTED);
                    cdf_graph_add_edge(g, w, v, CDF_EDGE_DIRECTED);
                    g->edges[v * p + u] = CDF_EDGE_NONE;
                    g->edges[v * p + w] = CDF_EDGE_NONE;
                    n_cons++;
                }
            }
        }
    }

    return n_cons;
}

/* ──────────────── Results & Validation ───────────────────────────── */

void cdf_pc_print_result(const CdfPCResult *res)
{
    if (!res) return;
    printf("=== PC Algorithm Result ===\n");
    printf("CI tests:      %d\n", res->n_ci_tests);
    printf("Edges removed: %d\n", res->n_edges_removed);
    printf("V-structures:  %d\n", res->n_vstructs);
    printf("Oriented:      %d\n", res->n_oriented);
    printf("Time:          %.4f sec\n", res->elapsed_sec);
}

int cdf_pc_shd(const CdfGraph *g1, const CdfGraph *g2)
{
    if (!g1 || !g2 || g1->p != g2->p) return -1;

    int p = g1->p, shd = 0;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            CdfEdgeType e1 = g1->edges[i * p + j];
            CdfEdgeType e2 = g2->edges[i * p + j];

            if (e1 == CDF_EDGE_NONE && e2 != CDF_EDGE_NONE) shd++;       /* missing */
            else if (e1 != CDF_EDGE_NONE && e2 == CDF_EDGE_NONE) shd++;   /* extra */
            else if (e1 != e2 && e1 != CDF_EDGE_NONE && e2 != CDF_EDGE_NONE) {
                /* Both edges present but different direction */
                /* Undirected vs directed counts as orientation error */
                shd++;
            }
        }
    }
    return shd;
}

bool cdf_pc_cpdag_is_valid(const CdfGraph *g)
{
    if (!g) return false;

    int p = g->p;

    /* Check no directed cycles */
    for (int i = 0; i < p; i++) {
        bool *vis = (bool*)calloc((size_t)p, sizeof(bool));
        bool *rec = (bool*)calloc((size_t)p, sizeof(bool));
        if (!vis || !rec) { free(vis); free(rec); return false; }

        int *stack = (int*)malloc((size_t)p * sizeof(int));
        int *state = (int*)calloc((size_t)p, sizeof(int));
        if (!stack || !state) { free(vis); free(rec); free(stack); free(state); return false; }

        int top = 0;
        stack[top++] = i;
        bool has_cycle = false;

        while (top > 0 && !has_cycle) {
            int cur = stack[top - 1];
            if (state[cur] == 0) {
                state[cur] = 1;
                for (int w = 0; w < p; w++) {
                    if (g->edges[cur * p + w] == CDF_EDGE_DIRECTED) {
                        if (state[w] == 1) { has_cycle = true; break; }
                        if (state[w] == 0) stack[top++] = w;
                    }
                }
            } else if (state[cur] == 1) {
                state[cur] = 2;
                top--;
            } else {
                top--;
            }
        }

        free(vis); free(rec); free(stack); free(state);
        if (has_cycle) return false;
    }

    /* Check edge consistency: no contradictory edge pairs */
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            if (i == j) continue;
            CdfEdgeType e_ij = g->edges[i * p + j];
            CdfEdgeType e_ji = g->edges[j * p + i];

            /* Directed edges must be one-way */
            if (e_ij == CDF_EDGE_DIRECTED && e_ji == CDF_EDGE_DIRECTED)
                return false;
            /* Bidirected must be symmetric */
            if (e_ij == CDF_EDGE_BIDIRECTED && e_ji != CDF_EDGE_BIDIRECTED)
                return false;
            /* Undirected must be symmetric */
            if (e_ij == CDF_EDGE_UNDIRECTED && e_ji != CDF_EDGE_UNDIRECTED)
                return false;
            /* Nondirected must be symmetric */
            if (e_ij == CDF_EDGE_NONDIR && e_ji != CDF_EDGE_NONDIR)
                return false;
        }
    }

    return true;
}

/* ── Utility: extract the skeleton from a CPDAG ──────────────────── */

void cdf_pc_extract_skeleton(const CdfGraph *cpdag, CdfGraph *skeleton)
{
    if (!cpdag || !skeleton || cpdag->p != skeleton->p) return;

    int p = cpdag->p;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            if (cdf_graph_has_edge(cpdag, i, j))
                cdf_graph_add_edge(skeleton, i, j, CDF_EDGE_UNDIRECTED);
        }
    }
}

/* ── Utility: count directed edges in a CPDAG ───────────────────── */

int cdf_pc_count_directed(const CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;
    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++)
            if (g->edges[i * p + j] == CDF_EDGE_DIRECTED)
                count++;
    return count;
}

/* ── Utility: count undirected edges ─────────────────────────────── */

int cdf_pc_count_undirected(const CdfGraph *g)
{
    if (!g) return 0;
    int p = g->p, count = 0;
    for (int i = 0; i < p; i++)
        for (int j = i + 1; j < p; j++)
            if (g->edges[i * p + j] == CDF_EDGE_UNDIRECTED)
                count++;
    return count;
}