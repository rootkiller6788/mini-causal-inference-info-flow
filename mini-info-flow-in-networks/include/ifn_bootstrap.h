#ifndef IFN_BOOTSTRAP_H
#define IFN_BOOTSTRAP_H
#include "ifn_core.h"

/* Bootstrap resampling and surrogate data methods for causal inference.
 *
 * Key methods:
 *   1. Block bootstrap — preserves temporal dependence structure
 *   2. IAAFT surrogates — preserves amplitude spectrum, destroys phase coupling
 *   3. Permutation testing — non-parametric significance for TE/Granger
 *   4. Confidence intervals — percentile and BCa methods
 *
 * Ref: Efron & Tibshirani (1993) An Introduction to the Bootstrap
 *      Schreiber & Schmitz (2000) Phys. Rev. E 62, 1952
 */

/* Block bootstrap: resample time series in contiguous blocks
 * to preserve serial dependence structure.
 * block_len: typical length 2-5 for weak dependence, 10-20 for strong
 * n_bootstrap: number of bootstrap replicates
 * Returns array of n_bootstrap * n bootstrap samples (row-major) */
double* ifn_block_bootstrap(const double* series, int n, int block_len,
    int n_bootstrap);

/* Stationary bootstrap: random block lengths from geometric distribution
 * with mean block length p. Better for unknown dependence structure.
 * Ref: Politis & Romano (1994) JASA 89, 1303 */
double* ifn_stationary_bootstrap(const double* series, int n, double p,
    int n_bootstrap);

/* IAAFT (Iterative Amplitude Adjusted Fourier Transform) surrogates.
 * Generates surrogate data with same amplitude spectrum but randomized phases.
 * Used for testing nonlinear causality (TE significance).
 * n_surrogates: typically 100-500 */
double** ifn_iaaft_surrogates(const double* signal, int n, int n_surrogates);

/* Permutation-based significance for transfer entropy.
 * Shuffles source series to destroy causal relationship, recomputes TE.
 * Returns p-value from empirical null distribution. */
double ifn_te_permutation_test(const IFN_TimeSeries* source,
    const IFN_TimeSeries* target, int k, int l,
    int n_permutations, const IFN_TransferResult* observed);

/* Percentile bootstrap confidence interval.
 * Returns [lower, upper] at given confidence level. */
void ifn_bootstrap_ci_percentile(const double* bootstrap_vals,
    int n_bootstrap, double alpha, double* lower, double* upper);

/* BCa (Bias-Corrected and accelerated) bootstrap confidence interval.
 * More accurate than percentile for skewed distributions.
 * Ref: Efron (1987) JASA 82, 171 */
void ifn_bootstrap_ci_bca(const double* bootstrap_vals, int n_bootstrap,
    const double* jackknife_vals, int n_original, double alpha,
    double* lower, double* upper);

/* Cross-validation for Granger causality lag selection.
 * Uses rolling-origin evaluation to find optimal max_lag.
 * Returns best lag that minimizes prediction error. */
int ifn_granger_cv_lag(const IFN_TimeSeries* x, const IFN_TimeSeries* y,
    int max_lag, int n_folds);

#endif
