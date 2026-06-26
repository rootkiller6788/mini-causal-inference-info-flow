#ifndef TE_SCHREIBER_H
#define TE_SCHREIBER_H
#include "te_core.h"

/* ================================================================
 * Schreiber Transfer Entropy (2000)
 * TE_{X->Y} = Σ p(y_{n+1}, y_n^{(k)}, x_n^{(l)}) log [p(y_{n+1}|y_n^{(k)},x_n^{(l)}) / p(y_{n+1}|y_n^{(k)})]
 * The original binned-histogram TE estimator.
 * ================================================================ */

/* ---- Core Schreiber Estimators ---- */
TEResult te_schreiber_compute(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int delay);
TEResult te_schreiber_binning(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_bins);

/* ---- Information-Theoretic Building Blocks ---- */
double te_schreiber_entropy_rate(const TEBinnedSeries *bs, int k);
double te_schreiber_mutual_info(const TEBinnedSeries *x, const TEBinnedSeries *y, int lag);

/* ---- Lag Scanning ---- */
void te_schreiber_scan_lag(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int max_lag, double *te_values);

/* ---- Multiscale TE (coarse-graining over τ scales) ---- */
void te_schreiber_multiscale(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_scales, double *te_scales);

/* ---- Time-Varying TE (sliding window) ---- */
int te_schreiber_sliding_window(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int window, int step, double *te_tv, int max_out);

/* ---- Embedding Optimization ---- */
void te_schreiber_optimal_embedding(const TETimeSeries *x, const TETimeSeries *y, int max_k, int max_l, int *best_k, int *best_l);
TEResult te_schreiber_variable_embedding(const TETimeSeries *x, const TETimeSeries *y, int k_min, int k_max, int l_min, int l_max);

/* ---- Directionality and Normalization ---- */
double te_schreiber_directionality(const TEResult *r);
double te_schreiber_uncertainty_coefficient(const TEResult *r);
double te_schreiber_normalized(const TEResult *r);

/* ---- Robust Estimation Methods ---- */
TEResult te_schreiber_robust_binning(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int min_bins, int max_bins, int n_bin_sizes);
TEResult te_schreiber_theiler(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_bins, int theiler_w);

/* ---- Integrated TE (sum over embedding dimensions) ---- */
double te_schreiber_integrated(const TETimeSeries *x, const TETimeSeries *y, int max_k, int max_l);

#endif
