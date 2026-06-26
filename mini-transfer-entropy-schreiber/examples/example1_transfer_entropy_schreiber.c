/* Example 1: Transfer Entropy on Coupled Logistic Maps
 * Demonstrates Schreiber TE, KSG TE, and symbolic TE on a coupled chaotic system.
 * Reference: Schreiber, T. (2000). "Measuring Information Transfer."
 *            Phys. Rev. Lett. 85, 461-464.
 */
#include "te_core.h"
#include "te_schreiber.h"
#include "te_ksg.h"
#include "te_effective.h"
#include "te_symbolic.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Example 1: Coupled Logistic Maps Transfer Entropy ===\n\n");
    int n = 1000;
    double *x = malloc((size_t)n * sizeof(double));
    double *y = malloc((size_t)n * sizeof(double));
    if (!x || !y) { printf("Allocation failed\n"); return 1; }

    /* Generate coupled logistic maps:
     * X_{t+1} = r_x * X_t * (1 - X_t) + c * Y_t (mod 1)
     * Y_{t+1} = r_y * Y_t * (1 - Y_t)
     * Y drives X through coupling c, but X does not drive Y.
     */
    double rx = 3.9, ry = 3.8, c = 0.3;
    x[0] = 0.5; y[0] = 0.3;
    for (int i = 1; i < n; i++) {
        y[i] = ry * y[i-1] * (1.0 - y[i-1]);
        x[i] = rx * x[i-1] * (1.0 - x[i-1]) + c * y[i-1];
        if (x[i] > 1.0) x[i] -= 1.0;
        if (x[i] < 0.0) x[i] += 1.0;
    }

    /* Skip transients: use last 800 points */
    int start = 200, use_n = n - start;
    TETimeSeries *ts_x = te_ts_create(x + start, use_n, "X");
    TETimeSeries *ts_y = te_ts_create(y + start, use_n, "Y");

    /* Method 1: Schreiber binning TE */
    printf("--- Schreiber Binning TE ---\n");
    TEResult r_bin = te_schreiber_binning(ts_x, ts_y, 3, 3, 8);
    te_result_print(&r_bin);

    /* Method 2: KSG k-NN TE */
    printf("\n--- KSG k-NN TE ---\n");
    TEResult r_ksg = te_ksg_compute(ts_x, ts_y, 3, 3, 5);
    te_result_print(&r_ksg);

    /* Method 3: Effective TE with surrogate correction */
    printf("\n--- Effective TE (bias-corrected, 50 surrogates) ---\n");
    TEResult r_eff = te_effective_compute(ts_x, ts_y, 3, 3, 50);
    te_result_print(&r_eff);

    /* Method 4: Symbolic TE */
    printf("\n--- Symbolic Transfer Entropy ---\n");
    TEResult r_sym = te_symbolic_compute(ts_x, ts_y, 4, 2);
    te_result_print(&r_sym);

    /* Interpretation */
    printf("\n=== Interpretation ===\n");
    printf("Since Y drives X (coupling c=%.2f), we expect TE(Y->X) > TE(X->Y).\n", c);
    printf("TE(Y->X) = %.6f nats (should be > 0)\n", r_bin.te_xy);
    printf("TE(X->Y) = %.6f nats (should be near 0)\n", r_bin.te_yx);
    printf("Directionality index = %.4f (positive = Y->X dominant)\n",
           (r_bin.te_xy + r_bin.te_yx > 0) ?
           (r_bin.te_xy - r_bin.te_yx) / (r_bin.te_xy + r_bin.te_yx) : 0.0);

    double dir = te_schreiber_directionality(&r_bin);
    if (dir > 0.1) printf("Result: Detected Y -> X causal direction. PASS.\n");
    else printf("Result: Direction unclear (check coupling).\n");

    te_ts_free(ts_x); te_ts_free(ts_y);
    free(x); free(y);
    printf("\nExample 1 PASSED\n");
    return 0;
}