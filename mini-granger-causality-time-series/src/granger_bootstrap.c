/* granger_bootstrap.c -- Bootstrap inference for Granger causality.
 *
 * Parametric and non-parametric bootstrap methods for constructing
 * confidence intervals and correcting finite-sample bias in
 * Granger causality statistics.
 *
 * Knowledge points:
 * L4: Bootstrap consistency theorem (Efron & Tibshirani 1993)
 * L5: Parametric bootstrap (residual resampling)
 * L5: Non-parametric bootstrap (block bootstrap for time series)
 * L5: Bias-corrected confidence intervals
 * L8: Advanced inference for Granger causality
 *
 * Reference: Efron & Tibshirani (1993). An Introduction to the Bootstrap.
 *            Hacker & Hatemi-J (2006), Applied Economics, 38(13):1489-1500.
 */

#include "granger_test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Parametric bootstrap Granger test.
 * Fits restricted VAR model (Y on Y-lags only) and generates B surrogate
 * datasets by resampling residuals, refitting the model each time.
 * Returns bootstrap p-value: fraction of surrogate F-stats >= observed. */
typedef struct {
    double observed_F;
    double bootstrap_p_value;
    double* surrogate_F;
    int n_bootstrap;
    double ci_lower;
    double ci_upper;
    double bias_correction;
} BootstrapGrangerResult;

BootstrapGrangerResult* boot_create(int n_bootstrap) {
    BootstrapGrangerResult* b = calloc(1, sizeof(BootstrapGrangerResult));
    if (b) {
        b->n_bootstrap = n_bootstrap;
        b->surrogate_F = malloc((size_t)n_bootstrap * sizeof(double));
    }
    return b;
}

void boot_free(BootstrapGrangerResult* b) {
    if (!b) return;
    free(b->surrogate_F);
    free(b);
}

/* Residual-based parametric bootstrap.
 * The restricted model (no X) is fitted, residuals are collected,
 * then resampled to generate B bootstrap time series under H0.
 * The F-statistic is recomputed on each bootstrap dataset.
 * Complexity: O(B * (length * max_lag + max_lag^3)). */
int boot_parametric(BootstrapGrangerResult* b, const double* x, const double* y,
    int length, int max_lag, double alpha) {
    (void)alpha;
    if (!b || !x || !y || length < max_lag + 4) return -1;
    int B = b->n_bootstrap;
    /* Fit restricted model: Y_t = c + sum_{i=1}^p a_i Y_{t-i} + eps_t */
    int p = max_lag;
    int n = length - p, k_r = 1 + p, k_ur = 1 + 2 * p;
    Matrix* X_r = mat_create(n, k_r);
    double* Y_vec = malloc((size_t)n * sizeof(double));
    double* beta_restricted = malloc((size_t)k_ur * sizeof(double));
    double* resid_r = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        X_r->data[i * k_r + 0] = 1.0;
        for (int lag = 0; lag < p; lag++)
            X_r->data[i * k_r + 1 + lag] = y[i + p - 1 - lag];
        Y_vec[i] = y[i + p];
    }
    double rss_r = ols_estimate(X_r, Y_vec, beta_restricted, resid_r, n, k_r);

    /* Center residuals */
    double mean_eps = 0.0;
    for (int i = 0; i < n; i++) mean_eps += resid_r[i];
    mean_eps /= n;
    for (int i = 0; i < n; i++) resid_r[i] -= mean_eps;

    /* Observed F from unrestricted model */
    Matrix* X_ur = mat_create(n, k_ur);
    double* beta_unrestricted = malloc((size_t)k_ur * sizeof(double));
    for (int i = 0; i < n; i++) {
        X_ur->data[i * k_ur + 0] = 1.0;
        for (int lag = 0; lag < p; lag++) {
            X_ur->data[i * k_ur + 1 + lag] = y[i + p - 1 - lag];
            X_ur->data[i * k_ur + 1 + p + lag] = x[i + p - 1 - lag];
        }
    }
    double rss_ur = ols_estimate(X_ur, Y_vec, beta_unrestricted, resid_r, n, k_ur);
    if (rss_ur < 1e-15) rss_ur = 1e-15;
    double F_obs = ((rss_r - rss_ur) / p) / (rss_ur / (n - k_ur));
    if (F_obs < 0.0) F_obs = 0.0;
    b->observed_F = F_obs;

    /* Bootstrap loop */
    double* y_boot = malloc((size_t)(length) * sizeof(double));
    double* resid_surr = malloc((size_t)n * sizeof(double));
    for (int b_iter = 0; b_iter < B; b_iter++) {
        /* Resample residuals with replacement */
        for (int i = 0; i < n; i++)
            resid_surr[i] = resid_r[rand() % n];
        /* Generate bootstrap time series under H0 (no X effect) */
        for (int i = 0; i < p; i++) y_boot[i] = y[i];
        for (int i = 0; i < n; i++) {
            double y_hat = beta_restricted[0];
            for (int lag = 0; lag < p; lag++)
                y_hat += beta_restricted[1 + lag] * y_boot[i + p - 1 - lag];
            y_boot[i + p] = y_hat + resid_surr[i];
        }
        /* Recompute Granger F on bootstrap sample */
        int n_eff = length - p;
        double* Y_b = malloc((size_t)n_eff * sizeof(double));
        Matrix* X_ur_b = mat_create(n_eff, k_ur);
        Matrix* X_r_b = mat_create(n_eff, k_r);
        double* beta_b = malloc((size_t)k_ur * sizeof(double));
        double* resid_b = malloc((size_t)n_eff * sizeof(double));
        for (int i = 0; i < n_eff; i++) {
            X_r_b->data[i * k_r + 0] = 1.0;
            X_ur_b->data[i * k_ur + 0] = 1.0;
            for (int lag = 0; lag < p; lag++) {
                X_r_b->data[i * k_r + 1 + lag] = y_boot[i + p - 1 - lag];
                X_ur_b->data[i * k_ur + 1 + lag] = y_boot[i + p - 1 - lag];
                X_ur_b->data[i * k_ur + 1 + p + lag] = x[i + p - 1 - lag];
            }
            Y_b[i] = y_boot[i + p];
        }
        double rss_r_b = ols_estimate(X_r_b, Y_b, beta_b, resid_b, n_eff, k_r);
        double rss_ur_b = ols_estimate(X_ur_b, Y_b, beta_b, resid_b, n_eff, k_ur);
        if (rss_ur_b < 1e-15) rss_ur_b = 1e-15;
        double F_b = ((rss_r_b - rss_ur_b) / p) / (rss_ur_b / (n_eff - k_ur));
        b->surrogate_F[b_iter] = (F_b > 0) ? F_b : 0.0;
        mat_free(X_ur_b); mat_free(X_r_b);
        free(Y_b); free(beta_b); free(resid_b);
    }
    /* Compute bootstrap p-value */
    int n_exceed = 0;
    for (int b_iter = 0; b_iter < B; b_iter++)
        if (b->surrogate_F[b_iter] >= F_obs) n_exceed++;
    b->bootstrap_p_value = (double)(n_exceed + 1) / (B + 1);

    /* Percentile confidence interval */
    double* sorted = malloc((size_t)B * sizeof(double));
    memcpy(sorted, b->surrogate_F, (size_t)B * sizeof(double));
    for (int i = 0; i < B - 1; i++)
        for (int j = i + 1; j < B; j++)
            if (sorted[i] > sorted[j]) {
                double tmp = sorted[i]; sorted[i] = sorted[j]; sorted[j] = tmp;
            }
    int lo_idx = (int)(B * 0.025), hi_idx = (int)(B * 0.975);
    b->ci_lower = sorted[lo_idx];
    b->ci_upper = sorted[hi_idx];
    free(sorted);
    /* Bias correction: BC = F_est - (mean(F_boot) - F_est) */
    double mean_boot = 0.0;
    for (int b_iter = 0; b_iter < B; b_iter++) mean_boot += b->surrogate_F[b_iter];
    mean_boot /= B;
    b->bias_correction = F_obs - (mean_boot - F_obs);

    mat_free(X_r); mat_free(X_ur);
    free(Y_vec); free(beta_restricted); free(beta_unrestricted); free(resid_r);
    free(y_boot); free(resid_surr);
    return 0;
}

/* Block bootstrap for time series: preserves temporal dependence.
 * Uses the moving block bootstrap (Kunsch 1989).
 * Block length l is chosen as l ≈ n^(1/3). */
int boot_block(BootstrapGrangerResult* b, const double* x, const double* y,
    int length, int max_lag, double alpha, int block_len) {
    (void)alpha;
    if (!b || !x || !y || length < max_lag + 4) return -1;
    int p = max_lag, B = b->n_bootstrap;
    int n = length - p, k_r = 1 + p, k_ur = 1 + 2 * p;
    if (block_len <= 0) block_len = (int)pow((double)n, 1.0 / 3.0);
    if (block_len < 2) block_len = 2;
    int n_blocks = (n + block_len - 1) / block_len;

    /* Observed F */
    Matrix* X_r = mat_create(n, k_r);
    Matrix* X_ur = mat_create(n, k_ur);
    double* Y_vec = malloc((size_t)n * sizeof(double));
    double* beta_r = malloc((size_t)k_ur * sizeof(double));
    double* resid_r = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        X_r->data[i * k_r + 0] = 1.0;
        X_ur->data[i * k_ur + 0] = 1.0;
        for (int lag = 0; lag < p; lag++) {
            X_r->data[i * k_r + 1 + lag] = y[i + p - 1 - lag];
            X_ur->data[i * k_ur + 1 + lag] = y[i + p - 1 - lag];
            X_ur->data[i * k_ur + 1 + p + lag] = x[i + p - 1 - lag];
        }
        Y_vec[i] = y[i + p];
    }
    double rss_r_obs = ols_estimate(X_r, Y_vec, beta_r, resid_r, n, k_r);
    double rss_ur_obs = ols_estimate(X_ur, Y_vec, beta_r, resid_r, n, k_ur);
    if (rss_ur_obs < 1e-15) rss_ur_obs = 1e-15;
    b->observed_F = ((rss_r_obs - rss_ur_obs) / p) / (rss_ur_obs / (n - k_ur));
    if (b->observed_F < 0.0) b->observed_F = 0.0;

    /* Block bootstrap */
    double* x_boot = malloc((size_t)length * sizeof(double));
    double* y_boot = malloc((size_t)length * sizeof(double));
    for (int b_iter = 0; b_iter < B; b_iter++) {
        for (int blk = 0; blk < n_blocks; blk++) {
            int src = rand() % (n - block_len + 1);
            int dst = blk * block_len;
            for (int k = 0; k < block_len && dst + k < n; k++) {
                y_boot[dst + k + p] = y[src + k + p];
                x_boot[dst + k + p] = x[src + k + p];
            }
        }
        for (int k = 0; k < p; k++) {
            y_boot[k] = y[k]; x_boot[k] = x[k];
        }
        int n_eff = length - p;
        double F_b;
        {
            Matrix* X_r_b = mat_create(n_eff, k_r);
            Matrix* X_ur_b = mat_create(n_eff, k_ur);
            double* Y_b = malloc((size_t)n_eff * sizeof(double));
            double* beta = malloc((size_t)k_ur * sizeof(double));
            double* resid = malloc((size_t)n_eff * sizeof(double));
            for (int i = 0; i < n_eff; i++) {
                X_r_b->data[i * k_r + 0] = 1.0; X_ur_b->data[i * k_ur + 0] = 1.0;
                for (int lag = 0; lag < p; lag++) {
                    X_r_b->data[i * k_r + 1 + lag] = y_boot[i + p - 1 - lag];
                    X_ur_b->data[i * k_ur + 1 + lag] = y_boot[i + p - 1 - lag];
                    X_ur_b->data[i * k_ur + 1 + p + lag] = x_boot[i + p - 1 - lag];
                }
                Y_b[i] = y_boot[i + p];
            }
            double rss_r_b = ols_estimate(X_r_b, Y_b, beta, resid, n_eff, k_r);
            double rss_ur_b = ols_estimate(X_ur_b, Y_b, beta, resid, n_eff, k_ur);
            if (rss_ur_b < 1e-15) rss_ur_b = 1e-15;
            F_b = ((rss_r_b - rss_ur_b) / p) / (rss_ur_b / (n_eff - k_ur));
            if (F_b < 0) F_b = 0;
            mat_free(X_r_b); mat_free(X_ur_b);
            free(Y_b); free(beta); free(resid);
        }
        b->surrogate_F[b_iter] = F_b;
    }
    int n_exceed = 0;
    for (int b_iter = 0; b_iter < B; b_iter++)
        if (b->surrogate_F[b_iter] >= b->observed_F) n_exceed++;
    b->bootstrap_p_value = (double)(n_exceed + 1) / (B + 1);

    free(x_boot); free(y_boot);
    mat_free(X_r); mat_free(X_ur);
    free(Y_vec); free(beta_r); free(resid_r);
    return 0;
}

void boot_print(const BootstrapGrangerResult* b) {
    if (!b) return;
    printf("Bootstrap Granger (B=%d):\n", b->n_bootstrap);
    printf("  F_obs=%.4f, boot-p=%.4f, CI=[%.4f, %.4f]\n",
        b->observed_F, b->bootstrap_p_value, b->ci_lower, b->ci_upper);
    printf("  bias_corrected_F=%.4f\n", b->bias_correction);
}