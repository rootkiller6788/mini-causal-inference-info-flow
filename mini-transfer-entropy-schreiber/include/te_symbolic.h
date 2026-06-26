#ifndef TE_SYMBOLIC_H
#define TE_SYMBOLIC_H
#include "te_core.h"

/* ================================================================
 * Symbolic Transfer Entropy (Staniek & Lehnertz, 2008)
 * Uses ordinal (permutation) patterns instead of binning.
 * Robust to non-stationarity and amplitude variations.
 * ================================================================ */

/* ---- Core Symbolic TE ---- */
TEResult te_symbolic_compute(const TETimeSeries *x, const TETimeSeries *y, int order, int delay);
int      te_symbolic_permutation_encode(const double *window, int order, int *pattern);
void     te_symbolic_pattern_count(const TETimeSeries *ts, int order, int delay, int *counts, int n_patterns);

/* ---- Permutation Entropy (Bandt & Pompe, 2002) ---- */
double te_symbolic_permutation_entropy(const TETimeSeries *ts, int order, int delay);
double te_symbolic_weighted_permutation_entropy(const TETimeSeries *ts, int order, int delay);
double te_symbolic_conditional_pe(const TETimeSeries *ts, int order, int delay);

/* ---- Symbolic Transfer Entropy ---- */
double te_symbolic_transfer_entropy(const TETimeSeries *x, const TETimeSeries *y, int order, int delay);
double te_symbolic_bias_corrected(const TETimeSeries *x, const TETimeSeries *y, int order, int delay, int n_shuffles);

/* ---- Multiscale Symbolic Analysis ---- */
void te_symbolic_multiscale_pe(const TETimeSeries *ts, int order, int delay, int n_scales, double *pe_scales);

/* ---- Transition Matrix & AAPE ---- */
int    te_symbolic_transition_matrix(const TETimeSeries *ts, int order, int delay, double *trans_matrix);
double te_symbolic_aape(const TETimeSeries *ts, int order, int delay, double A);

#endif
