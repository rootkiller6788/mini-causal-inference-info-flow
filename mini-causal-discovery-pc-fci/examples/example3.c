/* example3.c — PC vs FCI comparison and orientation rules
 *
 * Demonstrates: PC on a simple DAG, FCI on the same data,
 * and compares the outputs. Also shows the effect of
 * Meek's orientation rules (R1-R4).
 */

#include "cdf_core.h"
#include "cdf_citest.h"
#include "cdf_graph.h"
#include "cdf_orientation.h"
#include "cdf_pc.h"
#include "cdf_fci.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    printf("=== Example 3: PC vs FCI Comparison ===\n\n");

    int N = 400, p = 4;

    /* Generate data from DAG: X1 → X2, X1 → X3, X2 → X4, X3 → X4 */
    double *data = (double*)malloc((size_t)p * N * sizeof(double));
    if (!data) { printf("Alloc failed\n"); return 1; }

    for (int i = 0; i < N; i++) {
        double e1 = ((double)(i * 7 % 1000) / 1000.0 - 0.5);
        double e2 = ((double)(i * 13 % 1000) / 1000.0 - 0.5);
        double e3 = ((double)(i * 17 % 1000) / 1000.0 - 0.5);
        double e4 = ((double)(i * 19 % 1000) / 1000.0 - 0.5);

        data[0 * N + i] = e1;                        /* X1 */
        data[1 * N + i] = 0.8 * data[0 * N + i] + e2; /* X2 = 0.8*X1 + e2 */
        data[2 * N + i] = 0.6 * data[0 * N + i] + e3; /* X3 = 0.6*X1 + e3 */
        data[3 * N + i] = 0.5 * data[1 * N + i]       /* X4 = 0.5*X2 + 0.4*X3 + e4 */
                         + 0.4 * data[2 * N + i] + e4;
    }

    CdfDataset *ds = cdf_dataset_create(data, N, p);
    free(data);
    if (!ds) { printf("Dataset failed\n"); return 1; }

    printf("True DAG: X1 → X2, X1 → X3, X2 → X4, X3 → X4\n\n");

    /* ── PC Algorithm ── */
    printf("=== PC Algorithm ===\n");
    CdfPCConfig pc_cfg = cdf_pc_config_default();
    pc_cfg.alpha = 0.05;

    CdfPCResult *pc = cdf_pc_run(ds, &pc_cfg);
    if (pc) {
        cdf_pc_print_result(pc);
        printf("PC output CPDAG:\n");
        cdf_graph_print(pc->graph);
    }

    /* ── FCI Algorithm ── */
    printf("\n=== FCI Algorithm ===\n");
    CdfFCIConfig fci_cfg = cdf_fci_config_default();
    fci_cfg.alpha = 0.05;

    CdfFCIResult *fci = cdf_fci_run(ds, &fci_cfg);
    if (fci) {
        cdf_fci_print_result(fci);
        printf("FCI output PAG:\n");
        cdf_graph_print(fci->graph);
    }

    /* ── Comparison ── */
    if (pc && fci && pc->graph && fci->graph) {
        int diff = cdf_fci_compare_pags(pc->graph, fci->graph);
        printf("\nPC-FCI edge differences: %d\n", diff);
        printf("(With no latents, FCI should nearly match PC)\n");
    }

    /* ── Orientation Rule Demo ── */
    printf("\n=== Orientation Rules Demo ===\n");
    CdfGraph *demo = cdf_graph_create(4);
    /* Set up a structure that exercises R1, R2, R3 */
    cdf_graph_add_edge(demo, 0, 1, CDF_EDGE_DIRECTED);
    cdf_graph_add_edge(demo, 0, 2, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(demo, 1, 2, CDF_EDGE_UNDIRECTED);
    cdf_graph_add_edge(demo, 2, 3, CDF_EDGE_UNDIRECTED);

    printf("Before orientation:\n");
    cdf_graph_print(demo);

    CdfOrientResult orient = cdf_orient_pc_rules(demo);
    cdf_orient_print_result(&orient);

    printf("After orientation:\n");
    cdf_graph_print(demo);

    cdf_graph_free(demo);

    /* Cleanup */
    if (pc) cdf_pc_result_free(pc);
    if (fci) cdf_fci_result_free(fci);
    cdf_dataset_free(ds);
    printf("\nDone.\n");
    return 0;
}