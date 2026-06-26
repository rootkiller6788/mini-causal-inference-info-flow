/* transfer_entropy_schreiber_app2.c -- L7 Application: Climate Teleconnection Analysis
 * Implements TE-based detection of causal links between climate time series
 * (e.g., ENSO -> global temperature, AMO -> regional precipitation).
 *
 * Reference: Runge, J. et al. (2019) "Detecting and quantifying causal
 *            associations in large nonlinear time series datasets."
 *            Science Advances 5(11), eaau4996.
 *
 * Also: Hlinka, J. et al. (2013) "Reliability of inference of directed climate
 *       networks using the IBS method." Clim. Dyn.
 */
#include "te_core.h"
#include "te_schreiber.h"
#include "te_ksg.h"
#include "te_effective.h"
#include "te_multivariate.h"
#include "te_symbolic.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Generate synthetic ENSO-like index (Nino 3.4) with ~2-7 year periodicity */
static void generate_enso(double *data, int n, double *phase) {
    for (int i = 0; i < n; i++) {
        /* Primary 4-year oscillation + seasonal modulation */
        double seasonal = sin(2.0 * 3.141592653589793 * (double)i / 12.0);
        double interannual = 1.5 * sin(2.0 * 3.141592653589793 * (double)i / 48.0
                                        + (*phase));
        /* Stochastic component (weather noise) */
        double noise = 0.3 * ((double)rand() / (double)RAND_MAX - 0.5);
        data[i] = seasonal + interannual + noise;
        *phase += 0.00001 * ((double)rand() / (double)RAND_MAX - 0.5);
    }
}

/* Generate regional temperature anomaly influenced by ENSO with lag */
static void generate_temp_response(const double *enso, double *temp,
                                    int n, double sensitivity, int lag_months) {
    double baseline = 0.0;
    for (int i = 0; i < n; i++) {
        double enso_effect = 0.0;
        if (i >= lag_months) enso_effect = sensitivity * enso[i - lag_months];
        /* AR(1) process for internal variability */
        double ar_noise = 0.5 * baseline + 0.2 * ((double)rand() / (double)RAND_MAX - 0.5);
        baseline = ar_noise;
        temp[i] = baseline + enso_effect;
    }
}

/* Climate network link detection using TE with significance testing.
 * Implements the PCMCI framework (Runge et al. 2019) simplified variant.
 */
int te_climate_link_detection(const double *var1, const double *var2,
                               int n_samples, double tau_max,
                               double threshold, int *best_lag,
                               double *te_strength, double *confidence) {
    TETimeSeries *ts1 = te_ts_create(var1, n_samples, "Var1");
    TETimeSeries *ts2 = te_ts_create(var2, n_samples, "Var2");
    if (!ts1 || !ts2) return -1;

    *best_lag = -1;
    *te_strength = 0.0;
    *confidence = 0.0;
    int link_found = 0;

    for (int lag = 1; lag <= (int)tau_max && lag < n_samples / 10; lag++) {
        TEResult r = te_schreiber_binning(ts1, ts2, 3, 3, 6);
        if (r.te_xy > *te_strength) {
            *te_strength = r.te_xy;
            *best_lag = lag;
        }
    }

    if (*te_strength > threshold) {
        /* Bootstrap confidence */
        double *boot_vals = malloc(100 * sizeof(double));
        if (boot_vals) {
            int n = n_samples;
            for (int b = 0; b < 100; b++) {
                double *xr = malloc((size_t)n * sizeof(double));
                double *yr = malloc((size_t)n * sizeof(double));
                if (!xr || !yr) { free(xr); free(yr); continue; }
                for (int i = 0; i < n; i++) {
                    int idx = rand() % n;
                    xr[i] = var1[idx]; yr[i] = var2[idx];
                }
                TETimeSeries xb = {.data=xr,.length=n,.name="xr"};
                TETimeSeries yb = {.data=yr,.length=n,.name="yr"};
                TEResult rb = te_schreiber_binning(&xb, &yb, 3, 3, 6);
                boot_vals[b] = rb.te_xy;
                free(xr); free(yr);
            }
            /* Count nulls above observed */
            int count = 0;
            for (int b = 0; b < 100; b++)
                if (boot_vals[b] > *te_strength) count++;
            *confidence = 1.0 - (double)count / 100.0;
            free(boot_vals);
        }
        if (*confidence > 0.95) link_found = 1;
    }

    te_ts_free(ts1);
    te_ts_free(ts2);
    return link_found;
}

/* Detect teleconnections in a multi-variable climate system */
void te_climate_teleconnections(TETimeSeries **variables, int n_vars,
                                 double tau_max, double threshold,
                                 double **te_matrix, int *n_links) {
    *n_links = 0;
    for (int i = 0; i < n_vars; i++) {
        for (int j = 0; j < n_vars; j++) {
            if (i == j) { te_matrix[i][j] = 0.0; continue; }
            int best_lag;
            double strength, conf;
            int link = te_climate_link_detection(
                variables[i]->data, variables[j]->data,
                variables[i]->length, tau_max, threshold,
                &best_lag, &strength, &conf);
            te_matrix[i][j] = link ? strength : 0.0;
            if (link) (*n_links)++;
        }
    }
}

#ifdef TE_APP2_MAIN
int main(void) {
    printf("=== Climate Teleconnection Analysis (L7 Application) ===\n\n");

    int n = 600;  /* 50 years of monthly data */
    double *enso = calloc((size_t)n, sizeof(double));
    double *temp = calloc((size_t)n, sizeof(double));
    double *precip = calloc((size_t)n, sizeof(double));
    if (!enso || !temp || !precip) { printf("Memory error\n"); return 1; }

    /* Generate synthetic climate data */
    double phase = 0.0;
    generate_enso(enso, n, &phase);
    generate_temp_response(enso, temp, n, 0.8, 3);
    generate_temp_response(enso, precip, n, -0.4, 6);

    /* Basic statistics */
    TETimeSeries *ts_enso = te_ts_create(enso, n, "ENSO");
    TETimeSeries *ts_temp = te_ts_create(temp, n, "Temperature");
    TETimeSeries *ts_precip = te_ts_create(precip, n, "Precipitation");

    printf("ENSO std: %.3f, mean: %.3f\n",
           te_ts_std(ts_enso), te_ts_mean(ts_enso));
    printf("Temp std: %.3f, mean: %.3f\n",
           te_ts_std(ts_temp), te_ts_mean(ts_temp));

    /* TE analysis */
    int best_lag;
    double te_strength, confidence;
    int link = te_climate_link_detection(
        temp, enso, n, 12.0, 0.01, &best_lag, &te_strength, &confidence);

    printf("\nENSO -> Temperature TE:\n");
    printf("  Strength:   %.6f nats\n", te_strength);
    printf("  Best lag:   %d months\n", best_lag);
    printf("  Confidence: %.2f%%\n", confidence * 100.0);
    printf("  Link found: %s\n", link ? "YES" : "NO");

    /* Cross-correlation for comparison */
    double xc = te_ts_cross_correlation(ts_temp, ts_enso, best_lag > 0 ? best_lag : 3);
    printf("  Cross-corr: %.4f\n", xc);

    printf("\n=== Climate TE analysis complete ===\n");
    te_ts_free(ts_enso); te_ts_free(ts_temp); te_ts_free(ts_precip);
    free(enso); free(temp); free(precip);
    return 0;
}
#endif