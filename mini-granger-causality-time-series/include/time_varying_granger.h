#ifndef TIME_VARYING_GRANGER_H
#define TIME_VARYING_GRANGER_H
#include "ts_core.h"

/* Time-Varying Granger Causality via Sliding Windows
 *
 * Definition: Granger causality computed over rolling/sliding windows
 * to detect temporal changes in causal relationships.
 *
 * Window types supported:
 *  - Sliding window (fixed length, step by 1)
 *  - Recursive expanding window
 *  - Rolling window with fixed step
 *
 * Reference: Hesse et al. (2003), J. Neurosci. Methods, 124(1):89-95.
 *            Seth (2010), J. Neurosci. Methods, 186(2):262-273.
 */

typedef struct {
    int n_windows;
    int window_length;
    int step_size;
    double* time_points;     /* [n_windows] center time of each window */
    double* f_statistics;    /* [n_windows] F-stat per window */
    double* p_values;        /* [n_windows] p-value per window */
    int* lag_orders;         /* [n_windows] selected lag order per window */
    double* aic_values;      /* [n_windows] AIC of unrestricted model */
    double* bic_values;      /* [n_windows] BIC of unrestricted model */
} TVGrangerResult;

typedef enum { WINDOW_SLIDING, WINDOW_EXPANDING, WINDOW_ROLLING } WindowType;

/* Compute time-varying Granger causality.
 * window_len: length of each window
 * step: distance between consecutive window starts
 * type: sliding, expanding, or rolling window
 * Complexity: O(n_windows * max_lag * window_len). */
TVGrangerResult* tvg_create(void);
void tvg_free(TVGrangerResult* r);
int tvg_compute(TVGrangerResult* r, const double* x, const double* y,
    int length, int window_len, int step, int max_lag,
    double alpha, WindowType type);
void tvg_print(const TVGrangerResult* r);

/* Detect change points in Granger causality (sudden shifts in F-stat).
 * Uses CUSUM-like statistic: cumulative sum of F-stat deviations.
 * Returns indices of detected change points, terminated by -1.
 * Complexity: O(n_windows). */
int* tvg_detect_changepoints(const TVGrangerResult* r, double threshold);

/* Compute the dominance interval: contiguous time region where
 * X→Y Granger causality is significant at level alpha.
 * Returns [start_idx, end_idx] or [-1,-1] if no dominance. */
void tvg_significant_interval(const TVGrangerResult* r,
    int* start_idx, int* end_idx);

#endif
