/* demo_causal_graph.c - Causal Graph Visualization Demo
 *
 * Demonstrates the PC and FCI algorithms on synthetic datasets.
 * Shows: skeleton learning, v-structure detection, orientation rules,
 * and d-separation verification. Generates an ASCII visualization
 * of the learned causal graph.
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

#define NSAMP 300

/* Generate data from a known 4-node DAG:
 *   0 -> 1 -> 3
 *   0 -> 2 -> 3
 * (a simple "diamond" structure: X0 causes X1 and X2, which cause X3)
 */
static double* gen_diamond(int N)
{
    double *d = (double*)calloc((size_t)4 * N, sizeof(double));
    if (!d) return NULL;
    for (int i = 0; i < N; i++) {
        double e0 = ((double)((i * 7  + 11) % 997) / 997.0 - 0.5);
        double e1 = ((double)((i * 13 + 23) % 997) / 997.0 - 0.5);
        double e2 = ((double)((i * 17 + 37) % 997) / 997.0 - 0.5);
        double e3 = ((double)((i * 19 + 43) % 997) / 997.0 - 0.5);
        d[0 * N + i] = e0;
        d[1 * N + i] = 0.7 * e0 + e1;
        d[2 * N + i] = 0.6 * e0 + e2;
        d[3 * N + i] = 0.5 * (0.7 * e0 + e1) + 0.4 * (0.6 * e0 + e2) + e3;
    }
    return d;
}

/* Generate data from collider: 0 -> 1 <- 2, 0 -> 3, 2 -> 3 */
static double* gen_collider_fork(int N)
{
    double *d = (double*)calloc((size_t)4 * N, sizeof(double));
    if (!d) return NULL;
    for (int i = 0; i < N; i++) {
        double e0 = ((double)((i * 7  + 11) % 997) / 997.0 - 0.5);
        double e1 = ((double)((i * 13 + 23) % 997) / 997.0 - 0.5);
        double e2 = ((double)((i * 17 + 37) % 997) / 997.0 - 0.5);
        double e3 = ((double)((i * 19 + 43) % 997) / 997.0 - 0.5);
        d[0 * N + i] = e0;
        d[2 * N + i] = e2;
        d[1 * N + i] = 0.7 * e0 + 0.6 * e2 + e1;  /* collider */
        d[3 * N + i] = 0.5 * e0 + 0.4 * e2 + e3;
    }
    return d;
}

/* Print a simple node diagram */
static void print_diagram(const CdfGraph *g, const char *title,
                           const char **names)
{
    int p = g->p;
    printf("\n  %s:\n", title);
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            CdfEdgeType e = g->edges[i * g->p + j];
            const char *sym = cdf_graph_edge_char(e);
            if (e != CDF_EDGE_NONE) {
                printf("    %s %s %s\n",
                       names ? names[i] : "?", sym,
                       names ? names[j] : "?");
            }
        }
    }
}

/* Analyse d-separation properties on the learned graph */
static void check_dseps(const CdfGraph *g, const char **names)
{
    (void)names;
    printf("\n  d-separation analysis:\n");

    /* Check a few key pairs */
    struct { int x, y, z[2], nz; const char *desc; } checks[] = {
        {0, 3, {1},    1, "0 _|_ 3 | {1}? (chain blocking)"},
        {0, 2, {1},    1, "0 _|_ 2 | {1}? (collider opening?)"},
        {1, 3, {0},    1, "1 _|_ 3 | {0}? (fork blocking)"},
        {0, 0, {0},    0, NULL}  /* sentinel */
    };

    for (int c = 0; checks[c].desc; c++) {
        bool dsep = cdf_graph_d_separated(g, checks[c].x, checks[c].y,
                                           checks[c].z, checks[c].nz);
        printf("    %s -> %s\n", checks[c].desc, dsep ? "YES" : "NO");
    }
}

int main(void)
{
    const char *names4[] = {"X0", "X1", "X2", "X3"};
    int N = NSAMP;

    printf("==========================================\n");
    printf("  Causal Discovery Demo - PC & FCI       \n");
    printf("==========================================\n\n");

    /* Demo 1: Diamond DAG */
    printf("[Demo 1] Diamond DAG: X0->X1, X0->X2, X1->X3, X2->X3\n");
    double *data1 = gen_diamond(N);
    CdfDataset *ds1 = cdf_dataset_create(data1, N, 4);
    free(data1);

    if (ds1) {
        CdfPCConfig pc_cfg = cdf_pc_config_default();
        pc_cfg.alpha = 0.05;
        CdfPCResult *pc = cdf_pc_run(ds1, &pc_cfg);
        if (pc) {
            print_diagram(pc->graph, "PC-learned CPDAG", names4);
            printf("  CI tests: %d, V-structures: %d, Oriented: %d\n",
                   pc->n_ci_tests, pc->n_vstructs, pc->n_oriented);
            cdf_pc_result_free(pc);
        }
        cdf_dataset_free(ds1);
    }

    /* Demo 2: Collider detection */
    printf("\n[Demo 2] Collider + Fork: X0->X1<-X2, X0->X3, X2->X3\n");
    double *data2 = gen_collider_fork(N);
    CdfDataset *ds2 = cdf_dataset_create(data2, N, 4);
    free(data2);

    if (ds2) {
        CdfPCConfig pc_cfg = cdf_pc_config_default();
        pc_cfg.alpha = 0.05;
        CdfPCResult *pc = cdf_pc_run(ds2, &pc_cfg);
        if (pc) {
            print_diagram(pc->graph, "PC-learned CPDAG", names4);
            printf("  CI tests: %d, V-structures: %d, Oriented: %d\n",
                   pc->n_ci_tests, pc->n_vstructs, pc->n_oriented);
            check_dseps(pc->graph, names4);
            cdf_pc_result_free(pc);
        }

        /* Also run FCI */
        CdfFCIConfig fci_cfg = cdf_fci_config_default();
        fci_cfg.alpha = 0.05;
        CdfFCIResult *fci = cdf_fci_run(ds2, &fci_cfg);
        if (fci) {
            print_diagram(fci->graph, "FCI-learned PAG", names4);
            cdf_fci_result_free(fci);
        }
        cdf_dataset_free(ds2);
    }

    /* Demo 3: Orientation rule trace */
    printf("\n[Demo 3] Orientation Rules Trace:\n");
    CdfGraph *g = cdf_graph_create(4);
    cdf_graph_add_edge(g, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(g, 0, 2, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(g, 1, 2, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(g, 2, 3, CDF_EDGE_UNDIRECTED);

    printf("  Before: 0->1, 0--2, 1--2, 2--3\n");
    CdfOrientResult r = cdf_orient_pc_rules(g);
    printf("  After PC rules: R1=%d, R2=%d, R3=%d, R4=%d, total=%d\n",
           r.rules_applied[1], r.rules_applied[2],
           r.rules_applied[3], r.rules_applied[4],
           r.total_oriented);
    cdf_graph_print(g);
    cdf_graph_free(g);

    printf("\n==========================================\n");
    printf("  Demo complete.\n");
    printf("==========================================\n");
    return 0;
}
