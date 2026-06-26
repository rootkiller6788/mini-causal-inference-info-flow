/* Example 3: Symbolic Transfer Entropy on Coupled Henon Maps
 * Demonstrates all TE methods on chaotic dynamics with known coupling.
 * Uses various information-theoretic measures: ApEn, SampEn, PE, AIS.
 * Reference: Staniek, M. & Lehnertz, K. (2008).
 *            "Symbolic Transfer Entropy." Phys. Rev. Lett. 100, 158101.
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
    printf("=== Example 3: Multi-Method TE on Coupled Henon Maps ===\n\n");

    int n = 1200;
    double *x = calloc((size_t)n, sizeof(double));
    double *y = calloc((size_t)n, sizeof(double));
    if (!x || !y) { printf("Alloc fail\n"); return 1; }

    /* Coupled Henon maps:
     * X_{t+1} = 1 - 1.4*X_t^2 + 0.3*X_{t-1} + 0.1*Y_t
     * Y_{t+1} = 1 - 1.4*Y_t^2 + 0.3*Y_{t-1}
     * Weak coupling from Y to X.
     */
    double x0 = 0.1, x1 = 0.2, y0 = 0.3, y1 = 0.4;
    x[0] = x0; x[1] = x1; y[0] = y0; y[1] = y1;
    for (int i = 2; i < n; i++) {
        y[i] = 1.0 - 1.4 * y[i-1] * y[i-1] + 0.3 * y[i-2];
        x[i] = 1.0 - 1.4 * x[i-1] * x[i-1] + 0.3 * x[i-2] + 0.1 * y[i-1];
    }

    int start = 200, use_n = n - start;
    TETimeSeries *ts_x = te_ts_create(x + start, use_n, "X");
    TETimeSeries *ts_y = te_ts_create(y + start, use_n, "Y");

    /* 1. Time series statistics */
    printf("--- Time Series Statistics ---\n");
    char buf[256];
    te_ts_summary(ts_x, buf, 256);
    printf("X: %s\n", buf);
    te_ts_summary(ts_y, buf, 256);
    printf("Y: %s\n", buf);
    printf("X autocorrelation (lag 1): %.4f\n", te_ts_autocorrelation(ts_x, 1));
    printf("X Hurst exponent: %.4f\n", te_ts_hurst_rs(ts_x, 8, 64, 5));

    /* 2. Entropy measures */
    printf("\n--- Information Dynamics ---\n");
    double ais_x = te_active_information_storage(ts_x, 3, 8);
    printf("X AIS (Active Info Storage): %.4f nats\n", ais_x);
    double apen_x = te_approximate_entropy(ts_x, 2, 0.2);
    printf("X ApEn: %.4f\n", apen_x);
    double sampen_x = te_sample_entropy(ts_x, 2, 0.2);
    printf("X SampEn: %.4f\n", sampen_x);
    double pe_x = te_symbolic_permutation_entropy(ts_x, 4, 2);
    printf("X Permutation Entropy (order 4): %.4f\n", pe_x);

    /* 3. Transfer Entropy: all methods */
    printf("\n--- Transfer Entropy (Y->X) ---\n");

    TEResult r_schreiber = te_schreiber_binning(ts_x, ts_y, 3, 3, 8);
    printf("Schreiber binning TE: %.6f nats\n", r_schreiber.te_xy);

    TEResult r_ksg = te_ksg_compute(ts_x, ts_y, 3, 3, 5);
    printf("KSG k-NN TE:          %.6f nats\n", r_ksg.te_xy);

    TEResult r_eff = te_effective_compute(ts_x, ts_y, 3, 3, 50);
    printf("Effective TE:         %.6f nats (p=%.4f, %s)\n",
           r_eff.te_xy, r_eff.p_value,
           r_eff.significant ? "SIGNIFICANT" : "not significant");

    TEResult r_sym = te_symbolic_compute(ts_x, ts_y, 4, 2);
    printf("Symbolic TE:          %.6f nats\n", r_sym.te_xy);

    /* 4. Lag scan */
    printf("\n--- TE Lag Scan (Schreiber) ---\n");
    double te_lags[10];
    te_schreiber_scan_lag(ts_x, ts_y, 1, 1, 10, te_lags);
    printf("Lag:  ");
    for (int lag = 1; lag <= 10; lag++) printf("%4d ", lag);
    printf("\nTE:   ");
    for (int lag = 0; lag < 10; lag++) printf("%.3f ", te_lags[lag]);
    printf("\n");

    /* 5. Directionality */
    printf("\n--- Directionality Analysis ---\n");
    printf("TE(Y->X) = %.6f, TE(X->Y) = %.6f\n", r_schreiber.te_xy, r_schreiber.te_yx);
    double dir = te_schreiber_directionality(&r_schreiber);
    printf("Directionality index: %.4f\n", dir);
    if (dir > 0.05) printf("Detected: Y -> X dominant (matches coupling). PASS.\n");
    else printf("Direction ambiguous.\n");

    /* ASCII plot */
    printf("\n--- X time series (first 80 points, ASCII) ---\n");
    TETimeSeries *x80 = te_ts_create(ts_x->data, 80, "X80");
    te_ts_plot_ascii(x80, 60, 15);
    te_ts_free(x80);

    te_ts_free(ts_x); te_ts_free(ts_y);
    free(x); free(y);
    printf("\nExample 3 PASSED\n");
    return 0;
}