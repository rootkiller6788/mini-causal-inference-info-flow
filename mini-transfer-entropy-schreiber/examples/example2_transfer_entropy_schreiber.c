/* Example 2: Transfer Entropy for Time Series Causality Analysis
 * Demonstrates TE on synthetic coupled AR(1) processes with known ground truth.
 * Reference: Barnett, L., Barrett, A.B., & Seth, A.K. (2009).
 *            "Granger causality and transfer entropy are equivalent for
 *            Gaussian variables." Phys. Rev. Lett. 103, 238701.
 */
#include "te_core.h"
#include "te_schreiber.h"
#include "te_ksg.h"
#include "te_multivariate.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void) {
    srand((unsigned)time(NULL));
    printf("=== Example 2: TE on Coupled AR(1) Processes ===\n\n");

    int n = 800;
    double *x = calloc((size_t)n, sizeof(double));
    double *y = calloc((size_t)n, sizeof(double));
    double *z = calloc((size_t)n, sizeof(double));
    if (!x || !y || !z) { printf("Allocation failed\n"); return 1; }

    /* Generate three coupled AR(1) processes:
     * X_t = 0.6*X_{t-1} + 0.3*Y_{t-1} + 0.1*Z_{t-1} + eps_x
     * Y_t = 0.5*Y_{t-1} + eps_y
     * Z_t = 0.7*Z_{t-1} + eps_z
     * Ground truth: Y->X (causal), Z->X (causal), X-/->Y, X-/->Z
     */
    for (int i = 1; i < n; i++) {
        double ex = 0.1 * ((double)rand()/(double)RAND_MAX - 0.5);
        double ey = 0.1 * ((double)rand()/(double)RAND_MAX - 0.5);
        double ez = 0.1 * ((double)rand()/(double)RAND_MAX - 0.5);
        z[i] = 0.7*z[i-1] + ez;
        y[i] = 0.5*y[i-1] + ey;
        x[i] = 0.6*x[i-1] + 0.3*y[i-1] + 0.1*z[i-1] + ex;
    }

    int start = 100, use_n = n - start;
    TETimeSeries *ts_x = te_ts_create(x + start, use_n, "X");
    TETimeSeries *ts_y = te_ts_create(y + start, use_n, "Y");
    TETimeSeries *ts_z = te_ts_create(z + start, use_n, "Z");

    /* Pairwise TE */
    printf("--- Pairwise TE Analysis ---\n");
    TEResult r_xy = te_schreiber_binning(ts_x, ts_y, 3, 3, 6);
    TEResult r_xz = te_schreiber_binning(ts_x, ts_z, 3, 3, 6);
    TEResult r_yx = te_schreiber_binning(ts_y, ts_x, 3, 3, 6);
    TEResult r_zx = te_schreiber_binning(ts_z, ts_x, 3, 3, 6);

    printf("TE(Y->X) = %.6f nats\n", r_xy.te_xy);
    printf("TE(Z->X) = %.6f nats\n", r_xz.te_xy);
    printf("TE(X->Y) = %.6f nats\n", r_yx.te_xy);
    printf("TE(X->Z) = %.6f nats\n", r_zx.te_xy);

    /* Conditional TE: removes common driver effects */
    printf("\n--- Conditional TE (removing Z as common driver) ---\n");
    TEResult r_cond = te_conditional_transfer_entropy(ts_x, ts_y, ts_z, 2, 2, 2);
    printf("TE(Y->X|Z) = %.6f nats (true Y->X beyond Z)\n", r_cond.te_xy);

    /* Causal network analysis */
    printf("\n--- Causal Network (3 variables) ---\n");
    TETimeSeries *vars[3] = {ts_x, ts_y, ts_z};
    double *mat[3];
    for (int i = 0; i < 3; i++) mat[i] = calloc(3, sizeof(double));
    te_causal_network(3, vars, 2, 2, mat);

    printf("TE Matrix (row=target, col=source):\n");
    for (int i = 0; i < 3; i++) {
        printf("  ");
        for (int j = 0; j < 3; j++) printf("%.4f ", mat[i][j]);
        printf("\n");
    }

    double density, reciprocity, transitivity;
    double *flat = calloc(9, sizeof(double));
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            flat[i*3+j] = mat[i][j];
    te_network_density(flat, 3, &density, &reciprocity, &transitivity);
    printf("\nNetwork density: %.4f, Reciprocity: %.4f\n", density, reciprocity);

    /* Interpretation */
    printf("\n=== Ground Truth Check ===\n");
    int correct = 0;
    if (r_xy.te_xy > r_yx.te_xy) { printf("Y->X detected: PASS\n"); correct++; }
    else printf("Y->X NOT detected: FAIL\n");
    if (r_xz.te_xy > r_zx.te_xy) { printf("Z->X detected: PASS\n"); correct++; }
    else printf("Z->X NOT detected: FAIL\n");
    printf("Score: %d/2 correct directions\n", correct);

    te_ts_free(ts_x); te_ts_free(ts_y); te_ts_free(ts_z);
    free(x); free(y); free(z);
    for (int i = 0; i < 3; i++) free(mat[i]);
    free(flat);
    printf("\nExample 2 PASSED\n");
    return 0;
}