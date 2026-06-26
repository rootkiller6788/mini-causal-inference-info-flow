/* test_cdf.c — Tests for PC/FCI causal discovery module
 *
 * Tests: graph operations, CI tests, skeleton, v-structures,
 * orientation rules, PC algorithm, FCI algorithm.
 *
 * All tests use standard assert() from <assert.h>.
 */

#include "cdf_core.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include "cdf_orientation.h"
#include "cdf_pc.h"
#include "cdf_fci.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define TOL 1e-6

/* Generate simple Gaussian data from a known DAG: X1 → X2 → X3 */
static double* make_chain3(int N) {
    double *data = (double*)malloc((size_t)3 * N * sizeof(double));
    if (!data) return NULL;
    for (int i = 0; i < N; i++) {
        double e1 = (double)(i * 7 % 100) / 100.0 - 0.5;
        double e2 = (double)(i * 13 % 100) / 100.0 - 0.5;
        double e3 = (double)(i * 17 % 100) / 100.0 - 0.5;
        data[0 * N + i] = e1;              /* X1 = e1 */
        data[1 * N + i] = 0.7 * e1 + e2;   /* X2 = 0.7*X1 + e2 */
        data[2 * N + i] = 0.6 * (0.7 * e1 + e2) + e3; /* X3 = 0.6*X2 + e3 */
    }
    return data;
}

int main(void)
{
    printf("=== mini-causal-discovery-pc-fci Tests ===\n\n");

    /* ── Test 1: Constants ────────────────────────────────────────── */
    printf("[1] Constants...\n");
    assert(CDF_MAX_NODES >= 32);
    assert(CDF_ALPHA_DEFAULT > 0.0);
    assert(CDF_ALPHA_DEFAULT < 1.0);
    assert(cdf_version() > 0);
    printf("    PASS\n");

    /* ── Test 2: Graph Lifecycle ──────────────────────────────────── */
    printf("[2] Graph lifecycle...\n");
    CdfGraph *g = cdf_graph_create(5);
    assert(g != NULL);
    assert(g->p == 5);

    cdf_graph_add_edge(g, 0, 1, CDF_EDGE_UNDIRECTED);
    assert(cdf_graph_has_edge(g, 0, 1));
    assert(cdf_graph_has_edge(g, 1, 0));  /* undirected symmetric */

    cdf_graph_add_edge(g, 1, 2, CDF_EDGE_DIRECTED);
    assert(cdf_graph_has_edge(g, 1, 2));
    assert(cdf_graph_edge_type(g, 1, 2) == CDF_EDGE_DIRECTED);
    assert(!cdf_graph_has_edge(g, 2, 1)); /* directed one-way */

    cdf_graph_remove_edge(g, 0, 1);
    assert(!cdf_graph_has_edge(g, 0, 1));

    int neighbors[5], n;
    cdf_graph_neighbors(g, 1, neighbors, &n);
    assert(n >= 1);  /* should have at least node 2 */

    cdf_graph_free(g);
    printf("    PASS\n");

    /* ── Test 3: Reachability ─────────────────────────────────────── */
    printf("[3] Reachability...\n");
    CdfGraph *g2 = cdf_graph_create(4);
    cdf_graph_add_edge(g2, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(g2, 1, 2, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(g2, 2, 3, CDF_EDGE_DIRECTED);

    assert(cdf_graph_reachable(g2, 0, 3));
    assert(!cdf_graph_reachable(g2, 3, 0));

    assert(cdf_graph_is_ancestor(g2, 0, 3));
    assert(!cdf_graph_is_ancestor(g2, 3, 0));

    cdf_graph_free(g2);
    printf("    PASS\n");

    /* ── Test 4: d-Separation ─────────────────────────────────────── */
    printf("[4] d-separation...\n");
    CdfGraph *dag = cdf_graph_create(4);
    /* Chain: 0→1→2→3 */
    cdf_graph_add_edge(dag, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(dag, 1, 2, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(dag, 2, 3, CDF_EDGE_DIRECTED);

    /* 0 and 2 should be d-separated by {1} */
    int Z1[] = {1};
    assert(cdf_graph_d_separated(dag, 0, 2, Z1, 1));
    /* 0 and 3 should be d-separated by {1} or {2} */
    assert(cdf_graph_d_separated(dag, 0, 3, Z1, 1));
    int Z2[] = {2};
    assert(cdf_graph_d_separated(dag, 0, 3, Z2, 1));

    /* Fork: 0←1→2 */
    CdfGraph *fork = cdf_graph_create(3);
    cdf_graph_add_edge(fork, 1, 0, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(fork, 1, 2, CDF_EDGE_DIRECTED);
    assert(cdf_graph_d_separated(fork, 0, 2, Z1, 1));

    /* Collider: 0→1←2 */
    CdfGraph *coll = cdf_graph_create(3);
    cdf_graph_add_edge(coll, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(coll, 2, 1, CDF_EDGE_DIRECTED);
    /* 0 and 2 are NOT d-separated by {} (empty set) — they are marginally dependent */
    assert(!cdf_graph_d_separated(coll, 0, 2, NULL, 0));
    /* 0 and 2 ARE d-separated by {1} — conditioning on collider opens path */
    assert(!cdf_graph_d_separated(coll, 0, 2, Z1, 1));

    cdf_graph_free(dag);
    cdf_graph_free(fork);
    cdf_graph_free(coll);
    printf("    PASS\n");

    /* ── Test 5: CI Tests ─────────────────────────────────────────── */
    printf("[5] CI tests...\n");
    /* Fisher z-transform */
    double r = 0.5;
    double z = cdf_citest_fisher_z_transform(r);
    double r_back = cdf_citest_fisher_z_inverse(z);
    assert(fabs(r - r_back) < TOL);

    /* Normal p-value */
    double pv = cdf_citest_normal_pvalue(1.96);
    assert(fabs(pv - 0.05) < 0.01);

    printf("    PASS\n");

    /* ── Test 6: Partial Correlation ──────────────────────────────── */
    printf("[6] Partial correlation...\n");
    /* Generate chain data X→Y→Z and check partial corr */
    int N = 200;
    double *data = make_chain3(N);
    assert(data != NULL);
    CdfDataset *ds = cdf_dataset_create(data, N, 3);
    assert(ds != NULL);

    /* ρ_{XZ|Y} should be ~0 for chain X→Y→Z */
    int Z_y[] = {1};
    double prho = cdf_citest_partial_corr(ds, 0, 2, Z_y, 1);
    assert(fabs(prho) < 0.3);  /* should be small */

    free(data);
    cdf_dataset_free(ds);
    printf("    PASS\n");

    /* ── Test 7: Skeleton ─────────────────────────────────────────── */
    printf("[7] Skeleton...\n");
    data = make_chain3(N);
    ds = cdf_dataset_create(data, N, 3);
    CdfGraph *skel = cdf_graph_create(3);

    CdfPCConfig cfg = cdf_pc_config_default();
    cfg.alpha = 0.05;
    int n_ci = cdf_graph_skeleton_pc(skel, ds, &cfg);
    assert(n_ci > 0);
    /* For chain X→Y→Z, edges X-Y and Y-Z should remain, X-Z should be removed */
    assert(cdf_graph_has_edge(skel, 0, 1));
    assert(cdf_graph_has_edge(skel, 1, 2));
    /* X-Z may or may not be removed depending on data */
    printf("    Skeleton CI tests: %d\n", n_ci);

    cdf_graph_free(skel);
    cdf_dataset_free(ds);
    free(data);
    printf("    PASS\n");

    /* ── Test 8: V-Structures ─────────────────────────────────────── */
    printf("[8] V-structures...\n");
    CdfGraph *vg = cdf_graph_create(3);
    /* Manually create unshielded collider: 0—1—2 with 0,2 not adjacent */
    cdf_graph_add_edge(vg, 0, 1, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(vg, 1, 2, CDF_EDGE_UNDIRECTED);
    /* Set SepSet(0,2) = {1} — so 1 IS in the separating set,
     * meaning we should NOT orient as v-structure */
    int idx = 0 * 3 + 2;
    vg->sepsets[idx].u = 0;
    vg->sepsets[idx].v = 2;
    vg->sepsets[idx].size = 1;
    vg->sepsets[idx].set = (int*)malloc(sizeof(int));
    vg->sepsets[idx].set[0] = 1;
    vg->sepset_count[idx] = 1;

    int n_v = cdf_graph_find_vstructures(vg);
    assert(n_v == 0);  /* 1 is in SepSet, so not a v-structure */
    cdf_graph_free(vg);

    /* Now create a fresh graph: 1 NOT in SepSet → should be v-structure */
    CdfGraph *vg2 = cdf_graph_create(3);
    cdf_graph_add_edge(vg2, 0, 1, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(vg2, 1, 2, CDF_EDGE_UNDIRECTED);
    int idx2 = 0 * 3 + 2;
    vg2->sepset_count[idx2] = 0;  /* empty sepset → v-structure */
    n_v = cdf_graph_find_vstructures(vg2);
    assert(n_v == 1);
    cdf_graph_free(vg2);
    printf("    PASS\n");

    /* ── Test 9: Orientation Rules ────────────────────────────────── */
    printf("[9] Orientation rules...\n");
    CdfGraph *og = cdf_graph_create(4);
    /* Set up: 0→1—2 with 0,2 not adjacent → R1 should orient 1→2 */
    cdf_graph_add_edge(og, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(og, 1, 2, CDF_EDGE_UNDIRECTED);

    int n_r1 = cdf_orient_rule_r1(og);
    assert(n_r1 >= 1);
    assert(cdf_graph_edge_type(og, 1, 2) == CDF_EDGE_DIRECTED);

    cdf_graph_free(og);
    printf("    PASS\n");

    /* ── Test 10: PC Algorithm ────────────────────────────────────── */
    printf("[10] PC algorithm...\n");
    data = make_chain3(N);
    ds = cdf_dataset_create(data, N, 3);

    CdfPCResult *pc = cdf_pc_run(ds, NULL);
    assert(pc != NULL);
    assert(pc->graph != NULL);
    assert(pc->graph->is_cpdag);

    cdf_pc_print_result(pc);

    /* CPDAG validity */
    assert(cdf_pc_cpdag_is_valid(pc->graph));

    cdf_pc_result_free(pc);
    cdf_dataset_free(ds);
    free(data);
    printf("    PASS\n");

    /* ── Test 11: FCI Algorithm ───────────────────────────────────── */
    printf("[11] FCI algorithm...\n");
    data = make_chain3(N);
    ds = cdf_dataset_create(data, N, 3);

    CdfFCIResult *fci = cdf_fci_run(ds, NULL);
    assert(fci != NULL);
    assert(fci->graph != NULL);
    assert(fci->graph->is_pag);

    cdf_fci_print_result(fci);

    /* PAG validity */
    assert(cdf_fci_pag_is_valid(fci->graph));

    cdf_fci_result_free(fci);
    cdf_dataset_free(ds);
    free(data);
    printf("    PASS\n");

    /* ── Test 12: RFCI ────────────────────────────────────────────── */
    printf("[12] RFCI algorithm...\n");
    data = make_chain3(N);
    ds = cdf_dataset_create(data, N, 3);

    CdfFCIResult *rfci = cdf_fci_rfci_run(ds, NULL);
    assert(rfci != NULL);
    assert(rfci->graph != NULL);

    cdf_fci_result_free(rfci);
    cdf_dataset_free(ds);
    free(data);
    printf("    PASS\n");

    /* ── Test 13: Graph Statistics / Edge Types ────────────────────── */
    printf("[13] Graph stats & edge types...\n");
    CdfGraph *g13 = cdf_graph_create(4);
    cdf_graph_add_edge(g13, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(g13, 0, 2, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(g13, 1, 3, CDF_EDGE_BIDIRECTED);
    cdf_graph_add_edge(g13, 2, 3, CDF_EDGE_PARTIAL_I);

    assert(cdf_graph_edge_char(CDF_EDGE_DIRECTED) != NULL);
    assert(cdf_graph_edge_char(CDF_EDGE_NONE) != NULL);

    /* Shielded check */
    assert(cdf_graph_is_shielded(g13, 0, 3));  /* 0—2—3 shields 0-3 */

    cdf_graph_free(g13);
    printf("    PASS\n");

    /* ── Test 14: Null Safety ─────────────────────────────────────── */
    printf("[14] Null pointer safety...\n");
    assert(cdf_dataset_create(NULL, 10, 3) == NULL);
    assert(cdf_graph_create(0) == NULL);
    assert(cdf_graph_create(CDF_MAX_NODES + 1) == NULL);
    assert(!cdf_graph_has_edge(NULL, 0, 1));
    assert(!cdf_graph_d_separated(NULL, 0, 1, NULL, 0));
    assert(cdf_pc_run(NULL, NULL) == NULL);
    assert(cdf_fci_run(NULL, NULL) == NULL);
    cdf_graph_free(NULL);
    cdf_dataset_free(NULL);
    cdf_pc_result_free(NULL);
    cdf_fci_result_free(NULL);
    cdf_citest_result_free(NULL);
    printf("    PASS\n");

    /* ── Test 15: Edge Cases ──────────────────────────────────────── */
    printf("[15] Edge cases...\n");
    /* Self-loop check: not allowed */
    CdfGraph *g15 = cdf_graph_create(3);
    cdf_graph_add_edge(g15, 0, 0, CDF_EDGE_DIRECTED);
    assert(!cdf_graph_has_edge(g15, 0, 0));  /* ignored */

    /* Bidirectional consistency in PAG */
    cdf_graph_add_edge(g15, 0, 1, CDF_EDGE_BIDIRECTED);
    assert(cdf_graph_edge_type(g15, 1, 0) == CDF_EDGE_BIDIRECTED);

    /* Orientation without creating cycles */
    assert(!cdf_orient_would_create_cycle(g15, 0, 1));  /* no existing path 1→0 */

    /* Count unoriented edges */
    int n_un = cdf_orient_count_unoriented(g15);
    assert(n_un >= 0);

    cdf_graph_free(g15);

    /* PAG→MAG skeleton conversion */
    CdfGraph *pag = cdf_graph_create(3);
    cdf_graph_add_edge(pag, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(pag, 0, 2, CDF_EDGE_NONDIR);
    CdfGraph *mag = cdf_graph_create(3);
    cdf_fci_pag_to_mag_skeleton(pag, mag);
    assert(cdf_graph_edge_type(mag, 0, 1) == CDF_EDGE_DIRECTED);

    cdf_graph_free(pag);
    cdf_graph_free(mag);
    printf("    PASS\n");

    printf("\n=== ALL TESTS PASSED (15 assert groups) ===\n");
    return 0;
}