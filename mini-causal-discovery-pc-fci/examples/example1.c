/* example1.c — PC Algorithm on synthetic chain data
 *
 * Demonstrates the PC algorithm discovering a simple causal chain
 * X1 → X2 → X3 from Gaussian data.
 */

#include "cdf_core.h"
#include "cdf_citest.h"
#include "cdf_pc.h"
#include "cdf_graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    printf("=== Example 1: PC Algorithm on Chain ===\n\n");

    int N = 500, p = 3;

    /* Generate data: X1 → X2 → X3 */
    double *data = (double*)malloc((size_t)p * N * sizeof(double));
    if (!data) { printf("Alloc failed\n"); return 1; }

    for (int i = 0; i < N; i++) {
        double e1 = ((double)(i * 7 % 1000) / 1000.0 - 0.5);
        double e2 = ((double)(i * 13 % 1000) / 1000.0 - 0.5);
        double e3 = ((double)(i * 17 % 1000) / 1000.0 - 0.5);
        data[0 * N + i] = e1;
        data[1 * N + i] = 0.8 * e1 + e2;   /* X2 depends on X1 */
        data[2 * N + i] = 0.6 * data[1 * N + i] + e3;
    }

    CdfDataset *ds = cdf_dataset_create(data, N, p);
    free(data);
    if (!ds) { printf("Dataset failed\n"); return 1; }

    /* Configure PC algorithm */
    CdfPCConfig cfg = cdf_pc_config_default();
    cfg.alpha = 0.05;
    cfg.verbose = false;

    printf("Running PC on %d samples, %d variables...\n", N, p);

    CdfPCResult *res = cdf_pc_run(ds, &cfg);
    if (!res) { printf("PC failed\n"); cdf_dataset_free(ds); return 1; }

    cdf_pc_print_result(res);
    printf("\nLearned CPDAG:\n");
    cdf_graph_print(res->graph);
    printf("\nGraph stats:\n");
    cdf_graph_stats(res->graph);

    cdf_pc_result_free(res);
    cdf_dataset_free(ds);
    printf("\nDone.\n");
    return 0;
}