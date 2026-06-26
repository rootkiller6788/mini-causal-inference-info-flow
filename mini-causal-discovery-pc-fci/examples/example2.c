/* example2.c — FCI Algorithm with latent confounder
 *
 * Demonstrates FCI discovering a causal structure with a
 * latent confounder (unobserved common cause).
 *
 * True model: X1 ← L → X2 (L is latent), X1 → X3 ← X2
 * FCI should output: X1 ↔ X2 (bidirected, indicating latent),
 *                     X1 → X3 ← X2
 */

#include "cdf_core.h"
#include "cdf_citest.h"
#include "cdf_fci.h"
#include "cdf_graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    printf("=== Example 2: FCI with Latent Confounder ===\n\n");

    int N = 500, p = 3;

    /* True model: L (unobserved) → X1, L → X2, X1 → X3, X2 → X3
     * Observed: X1, X2, X3 (p=3)
     * X1 and X2 are correlated through L but NOT directly connected.
     */
    double *data = (double*)malloc((size_t)p * N * sizeof(double));
    if (!data) { printf("Alloc failed\n"); return 1; }

    for (int i = 0; i < N; i++) {
        double L_val = ((double)(i * 11 % 1000) / 1000.0 - 0.5);
        double e1 = ((double)(i * 7 % 1000) / 1000.0 - 0.5);
        double e2 = ((double)(i * 13 % 1000) / 1000.0 - 0.5);
        double e3 = ((double)(i * 17 % 1000) / 1000.0 - 0.5);

        data[0 * N + i] = 0.8 * L_val + e1;     /* X1 = 0.8*L + e1 */
        data[1 * N + i] = 0.7 * L_val + e2;     /* X2 = 0.7*L + e2 */
        data[2 * N + i] = 0.5 * data[0 * N + i]  /* X3 = 0.5*X1 + 0.4*X2 + e3 */
                         + 0.4 * data[1 * N + i] + e3;
    }

    CdfDataset *ds = cdf_dataset_create(data, N, p);
    free(data);
    if (!ds) { printf("Dataset failed\n"); return 1; }

    /* Configure FCI */
    CdfFCIConfig cfg = cdf_fci_config_default();
    cfg.alpha = 0.05;
    cfg.max_path_length = 3;

    printf("Running FCI on %d samples, %d variables...\n", N, p);
    printf("(True: X1 ←[L]→ X2, X1 → X3 ← X2)\n\n");

    CdfFCIResult *res = cdf_fci_run(ds, &cfg);
    if (!res) { printf("FCI failed\n"); cdf_dataset_free(ds); return 1; }

    cdf_fci_print_result(res);
    printf("\nLearned PAG:\n");
    cdf_graph_print(res->graph);

    /* Show edge type legend */
    printf("\nEdge legend:\n");
    printf("  →  directed\n");
    printf("  ↔  bidirected (latent confounder)\n");
    printf("  ∘→ partially directed\n");
    printf("  ∘—∘ nondirected\n");

    cdf_fci_result_free(res);
    cdf_dataset_free(ds);
    printf("\nDone.\n");
    return 0;
}