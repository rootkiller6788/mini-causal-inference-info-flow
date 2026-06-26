/* granger_app_economics.c -- L7 Application: Granger causality in economics.
 *
 * Knowledge points:
 * L7: Economics application - GDP, money supply, interest rates
 *
 * Classic application of Granger causality: testing whether money supply
 * Granger-causes GDP (monetarist hypothesis, Friedman & Schwartz 1963).
 * Also tests interest rate -> GDP and GDP -> money supply causality.
 *
 * This uses simulated economic data based on stylized facts from:
 * Stock & Watson (2001), "Vector Autoregressions," JEP 15(4):101-115.
 * Christiano, Eichenbaum, Evans (1999), Handbook of Macroeconomics.
 *
 * Reference: Sims, C.A. (1972). AER, 62(4):540-552.
 */

#include "ts_core.h"
#include "granger_test.h"
#include "fevd.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Simulate a 3-variable macroeconomic system:
 * GDP(t)    = 0.3*GDP(t-1) + 0.1*M2(t-1) - 0.05*IR(t-1) + eps_gdp
 * M2(t)     = 0.2*GDP(t-1) + 0.5*M2(t-1) + eps_m2
 * IR(t)     = -0.1*GDP(t-1) + 0.7*IR(t-1) + eps_ir
 * This reflects: GDP affects money demand, money affects GDP with lag,
 * interest rates are persistent and respond to output.
 */
static TimeSeries* simulate_macro_data(int length) {
    TimeSeries* ts = ts_create(length, 3);
    if (!ts) return NULL;
    /* Initial values at steady state */
    ts->values[0] = 100.0;  /* GDP */
    ts->values[1] = 50.0;   /* Money M2 */
    ts->values[2] = 5.0;    /* Interest rate */
    for (int t = 1; t < length; t++) {
        double gdp_lag = ts->values[(t - 1) * 3 + 0];
        double m2_lag = ts->values[(t - 1) * 3 + 1];
        double ir_lag = ts->values[(t - 1) * 3 + 2];
        /* GDP equation: autoregressive + money effect - interest rate drag */
        ts->values[t * 3 + 0] = 0.3 * gdp_lag + 0.1 * m2_lag - 0.05 * ir_lag
            + 1.0 * ((double)rand() / RAND_MAX - 0.5);
        /* Money equation: responds to GDP + autocorrelated */
        ts->values[t * 3 + 1] = 0.2 * gdp_lag + 0.5 * m2_lag
            + 0.5 * ((double)rand() / RAND_MAX - 0.5);
        /* Interest rate: persistent, counter-cyclical */
        ts->values[t * 3 + 2] = -0.1 * gdp_lag + 0.7 * ir_lag
            + 0.2 * ((double)rand() / RAND_MAX - 0.5);
    }
    return ts;
}

int main(void) {
    printf("=== Granger Causality: Macroeconomic Application ===\n\n");
    srand(42);
    int T = 200;
    TimeSeries* ts = simulate_macro_data(T);
    printf("Simulated %d quarters of macroeconomic data.\n", T);
    printf("Variables: [0]=GDP, [1]=Money(M2), [2]=Interest Rate\n\n");
    /* Extract series */
    double* gdp = ts_extract_dim(ts, 0);
    double* m2  = ts_extract_dim(ts, 1);
    double* ir  = ts_extract_dim(ts, 2);
    /* Pairwise Granger tests */
    printf("--- Pairwise Granger Causality Tests ---\n");
    GrangerTestResult* gt = gt_create();
    /* Test: M2 -> GDP (Monetarist hypothesis) */
    gt_test(gt, m2, gdp, T, 4, 0.05);
    printf("M2 -> GDP:  F=%.3f p=%.4f %s (Monetarist hypothesis)\n",
        gt->f_statistic, gt->p_value, gt->is_significant ? "SIGNIFICANT" : "NOT sig");
    /* Test: GDP -> M2 (Money demand) */
    gt_test(gt, gdp, m2, T, 4, 0.05);
    printf("GDP -> M2:  F=%.3f p=%.4f %s (Money demand)\n",
        gt->f_statistic, gt->p_value, gt->is_significant ? "SIGNIFICANT" : "NOT sig");
    /* Test: IR -> GDP (Monetary transmission) */
    gt_test(gt, ir, gdp, T, 4, 0.05);
    printf("IR  -> GDP: F=%.3f p=%.4f %s (Monetary transmission)\n",
        gt->f_statistic, gt->p_value, gt->is_significant ? "SIGNIFICANT" : "NOT sig");
    gt_free(gt);
    /* VAR analysis */
    printf("\n--- VAR(2) Model Fit ---\n");
    VARModel* var = var_create(3, 2);
    var_fit(var, ts);
    printf("VAR(2) sigma2 = %.6f, AIC = %.3f, BIC = %.3f\n",
        var->sigma2, var->aic, var->bic);
    printf("Coefficient matrix A[1]:\n");
    printf("  [[%.3f, %.3f, %.3f],\n", mat_get(var->A[0],0,0), mat_get(var->A[0],0,1), mat_get(var->A[0],0,2));
    printf("   [%.3f, %.3f, %.3f],\n", mat_get(var->A[0],1,0), mat_get(var->A[0],1,1), mat_get(var->A[0],1,2));
    printf("   [%.3f, %.3f, %.3f]]\n", mat_get(var->A[0],2,0), mat_get(var->A[0],2,1), mat_get(var->A[0],2,2));
    /* FEVD analysis */
    printf("\n--- Forecast Error Variance Decomposition ---\n");
    FEVDResult* fevd = fevd_compute(var, 4);
    if (fevd) {
        fevd_print(fevd, 0);
        fevd_print(fevd, 3);
        printf("\nFEVD Interpretation:\n");
        printf("  Short-run (h=1): Own shocks dominate all variables.\n");
        printf("  Long-run (h=4): Cross-variable contributions emerge.\n");
        printf("  This reveals causal structure: which shocks drive\n");
        printf("  forecast uncertainty for each variable.\n");
        fevd_free(fevd);
    }
    var_free(var);
    ts_free(ts); free(gdp); free(m2); free(ir);
    printf("\n=== Application Complete ===\n");
    return 0;
}