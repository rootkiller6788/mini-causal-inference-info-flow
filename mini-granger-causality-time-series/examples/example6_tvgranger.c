/* Example 6: Time-Varying Granger Causality.
 * Demonstrates detection of changing causal relationships over time. */
#include "ts_core.h"
#include "time_varying_granger.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(void) {
    printf("=== Time-Varying Granger Causality ===\n");
    srand(123);
    int T = 500;
    /* Generate non-stationary coupled series:
     * First half: X -> Y (coupled), Second half: decoupled */
    double* x = malloc(T * sizeof(double));
    double* y = malloc(T * sizeof(double));
    x[0] = 0.0; y[0] = 0.0;
    for (int t = 1; t < T; t++) {
        x[t] = 0.6*x[t-1] + 0.2*((double)rand()/RAND_MAX-0.5);
        if (t < T/2)
            y[t] = 0.3*x[t-1] + 0.5*y[t-1] + 0.2*((double)rand()/RAND_MAX-0.5);
        else
            y[t] = 0.5*y[t-1] + 0.2*((double)rand()/RAND_MAX-0.5);
    }
    printf("Series: coupled (t<250), decoupled (t>=250)\n\n");
    TVGrangerResult* tvg = tvg_create();
    tvg_compute(tvg, x, y, T, 100, 25, 4, 0.05, WINDOW_SLIDING);
    tvg_print(tvg);
    /* Changepoint detection */
    int* cps = tvg_detect_changepoints(tvg, 0.7);
    printf("\nDetected changepoints: ");
    for (int i = 0; cps && cps[i] >= 0; i++)
        printf("t=%.0f ", tvg->time_points[cps[i]]);
    printf("\n");
    free(cps);
    /* Significant interval */
    int start, end;
    tvg_significant_interval(tvg, &start, &end);
    printf("Significant interval: [%.0f, %.0f]\n",
        start>=0 ? tvg->time_points[start] : -1.0,
        end>=0 ? tvg->time_points[end] : -1.0);
    tvg_free(tvg); free(x); free(y);
    return 0;
}