/* Example 4: FEVD and Causal Ordering Demo.
 * Demonstrates forecast error variance decomposition to understand
 * the dynamic causal structure of a multivariate time series. */
#include "ts_core.h"
#include "fevd.h"
#include <stdio.h>
#include <stdlib.h>
int main(void) {
    printf("=== FEVD: Causal Contribution Analysis ===\n");
    /* Simulate a bivariate system: Y1(t) = 0.7*Y1(t-1) + 0.2*Y2(t-1) + eps1 */
    /*                              Y2(t) =         0.5*Y2(t-1) + eps2   */
    /* So Y1 causes Y2 only indirectly through lagged Y2. */
    TimeSeries* ts = ts_simulate_var1_coupled(0.7, 0.2, 0.0, 0.5, 300);
    int best_p; double best_aic;
    var_order_select(ts, 8, 2, &best_p, &best_aic);
    printf("Optimal lag: %d (AIC=%.2f)\n", best_p, best_aic);
    VARModel* var = var_create(2, best_p);
    var_fit(var, ts);
    printf("VAR(%d) fitted: sigma2=%.6f\n", var->p, var->sigma2);
    FEVDResult* fevd = fevd_compute(var, 10);
    if (fevd) {
        printf("\nFEVD at horizon 1:\n");
        fevd_print(fevd, 0);
        printf("FEVD at horizon 5:\n");
        fevd_print(fevd, 4);
        printf("FEVD at horizon 10:\n");
        fevd_print(fevd, 9);
        printf("\nInterpretation:\n");
        printf("  Row i, Column j = proportion of var i variance from shock j\n");
        printf("  Diagonal dominance = own shocks matter most\n");
        printf("  Off-diagonal growth = increasing causal influence over time\n");
        fevd_free(fevd);
    }
    /* Generalized FEVD (order-invariant) */
    FEVDResult* gfevd = fevd_generalized(var, 5);
    if (gfevd) {
        printf("\nGeneralized FEVD (Pesaran-Shin 1998), horizon 5:\n");
        fevd_print(gfevd, 4);
        fevd_free(gfevd);
    }
    var_free(var); ts_free(ts);
    return 0;
}