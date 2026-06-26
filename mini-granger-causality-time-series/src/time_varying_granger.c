/* time_varying_granger.c -- Time-Varying Granger Causality via sliding windows.
 *
 * Computes Granger causality over rolling windows to detect temporal
 * changes in causal relationships. Essential for non-stationary time
 * series where causality may emerge, strengthen, weaken, or vanish over time.
 *
 * Knowledge points:
 * L2: Time-varying Granger causality concept
 * L3: Window functions (sliding, expanding, rolling)
 * L5: Recursive VAR estimation across windows
 * L5: CUSUM changepoint detection
 * L6: Detecting temporal causal transitions
 * L8: Non-stationary causal inference
 *
 * Reference: Hesse et al. (2003), J. Neurosci. Methods, 124(1):89-95.
 *            Seth (2010), J. Neurosci. Methods, 186(2):262-273.
 */

#include "time_varying_granger.h"
#include "granger_test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

TVGrangerResult* tvg_create(void) {
    return calloc(1, sizeof(TVGrangerResult));
}

void tvg_free(TVGrangerResult* r) {
    if (!r) return;
    free(r->time_points);
    free(r->f_statistics);
    free(r->p_values);
    free(r->lag_orders);
    free(r->aic_values);
    free(r->bic_values);
    free(r);
}

static int count_windows(int length, int window_len, int step, WindowType type) {
    if (type == WINDOW_SLIDING || type == WINDOW_ROLLING)
        return (length - window_len) / step + 1;
    if (type == WINDOW_EXPANDING)
        return length - window_len + 1;
    return 0;
}

static void window_bounds(int w_idx, int length, int window_len, int step,
    WindowType type, int* start, int* end) {
    if (type == WINDOW_EXPANDING) {
        *start = 0;
        *end = window_len + w_idx;
    } else {
        *start = w_idx * step;
        *end = *start + window_len;
    }
    if (*end > length) *end = length;
}

int tvg_compute(TVGrangerResult* r, const double* x, const double* y,
    int length, int window_len, int step, int max_lag,
    double alpha, WindowType type) {
    if (!r || !x || !y || length < window_len + 2 || step < 1 || max_lag < 1)
        return -1;
    int nw = count_windows(length, window_len, step, type);
    if (nw <= 0) return -1;
    r->n_windows = nw;
    r->window_length = window_len;
    r->step_size = step;
    r->time_points = malloc((size_t)nw * sizeof(double));
    r->f_statistics = malloc((size_t)nw * sizeof(double));
    r->p_values = malloc((size_t)nw * sizeof(double));
    r->lag_orders = malloc((size_t)nw * sizeof(int));
    r->aic_values = malloc((size_t)nw * sizeof(double));
    r->bic_values = malloc((size_t)nw * sizeof(double));
    if (!r->time_points || !r->f_statistics || !r->p_values ||
        !r->lag_orders || !r->aic_values || !r->bic_values) return -1;
    for (int w = 0; w < nw; w++) {
        int w_start, w_end;
        window_bounds(w, length, window_len, step, type, &w_start, &w_end);
        int w_len = w_end - w_start;
        /* Extract window data */
        double* xw = malloc((size_t)w_len * sizeof(double));
        double* yw = malloc((size_t)w_len * sizeof(double));
        for (int i = 0; i < w_len; i++) {
            xw[i] = x[w_start + i];
            yw[i] = y[w_start + i];
        }
        int best_lag = gt_optimal_lag(xw, yw, w_len, max_lag);
        GrangerTestResult gt;
        gt_test(&gt, xw, yw, w_len, best_lag, alpha);
        r->f_statistics[w] = gt.f_statistic;
        r->p_values[w] = gt.p_value;
        r->lag_orders[w] = best_lag;
        r->time_points[w] = (w_start + w_end) / 2.0;
        /* Compute AIC/BIC for this window's unrestricted model */
        int n = w_len - best_lag, k = 1 + 2 * best_lag;
        double rss = gt.rss_unrestricted;
        if (rss < 1e-15) rss = 1e-15;
        r->aic_values[w] = n * log(rss / n) + 2.0 * k;
        r->bic_values[w] = n * log(rss / n) + log((double)n) * k;
        free(xw); free(yw);
    }
    return 0;
}

void tvg_print(const TVGrangerResult* r) {
    if (!r) return;
    printf("Time-Varying Granger: %d windows (len=%d, step=%d)\n",
        r->n_windows, r->window_length, r->step_size);
    for (int w = 0; w < r->n_windows && w < 20; w++)
        printf("  t=%.0f: F=%.3f p=%.4f lag=%d AIC=%.1f\n",
            r->time_points[w], r->f_statistics[w], r->p_values[w],
            r->lag_orders[w], r->aic_values[w]);
    if (r->n_windows > 20) printf("  ... (%d more windows)\n", r->n_windows - 20);
}

/* CUSUM-based changepoint detection.
 * Computes cumulative sum of deviations from the mean F-statistic.
 * Changepoints occur where the CUSUM statistic is maximized.
 * Returns indices of detected changepoints, terminated by -1. */
int* tvg_detect_changepoints(const TVGrangerResult* r, double threshold) {
    if (!r || r->n_windows < 3) return NULL;
    int n = r->n_windows;
    /* Compute mean F-statistic */
    double mean_F = 0.0;
    for (int i = 0; i < n; i++) mean_F += r->f_statistics[i];
    mean_F /= n;
    /* Compute CUSUM */
    double* cusum = malloc((size_t)n * sizeof(double));
    double running = 0.0, max_cusum = 0.0;
    for (int i = 0; i < n; i++) {
        running += r->f_statistics[i] - mean_F;
        cusum[i] = running;
        if (fabs(running) > max_cusum) max_cusum = fabs(running);
    }
    /* Detect peaks exceeding threshold * max_cusum */
    int* cps = malloc((size_t)(n + 1) * sizeof(int));
    int ncp = 0;
    for (int i = 1; i < n - 1; i++)
        if (fabs(cusum[i]) > threshold * max_cusum &&
            fabs(cusum[i]) > fabs(cusum[i - 1]) &&
            fabs(cusum[i]) >= fabs(cusum[i + 1]))
            cps[ncp++] = i;
    cps[ncp] = -1;
    free(cusum);
    return cps;
}

/* Find the longest contiguous interval where p < alpha.
 * This identifies the dominant causal regime. */
void tvg_significant_interval(const TVGrangerResult* r,
    int* start_idx, int* end_idx) {
    *start_idx = -1; *end_idx = -1;
    if (!r) return;
    int best_len = 0, cur_start = -1;
    for (int i = 0; i < r->n_windows; i++) {
        if (r->p_values[i] < 0.05) {
            if (cur_start < 0) cur_start = i;
            int cur_len = i - cur_start + 1;
            if (cur_len > best_len) {
                best_len = cur_len;
                *start_idx = cur_start;
                *end_idx = i;
            }
        } else {
            cur_start = -1;
        }
    }
}