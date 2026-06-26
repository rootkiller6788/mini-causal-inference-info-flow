#ifndef TE_EFFECTIVE_H
#define TE_EFFECTIVE_H
#include "te_core.h"

/* ================================================================
 * Effective Transfer Entropy (Marschinski & Kantz, 2002)
 * ETE = TE(X->Y) - mean(TE(X_surrogate->Y))
 * ================================================================ */
TEResult te_effective_compute(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_surrogates);
TEResult te_effective_bias_corrected(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_shuffles);
double   te_effective_shuffle_bias(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_shuffles);
double   te_jackknife_te(const TETimeSeries *x, const TETimeSeries *y, int k, int l);

/* ---- Surrogate Generation Methods ---- */
void te_generate_iid_surrogates(const double *data, int n, double *surrogate);
void te_generate_ft_surrogates(const double *data, int n, double *surrogate);
void te_generate_iaaft_surrogates(const double *data, int n, double *surrogate, int max_iter);
void te_generate_block_bootstrap(const double *data, int n, int block_len, double *surrogate);
void te_generate_stationary_bootstrap(const double *data, int n, double p, double *surrogate);
void te_generate_circular_shift(const double *data, int n, double *surrogate);

/* ---- Statistical Testing ---- */
double te_zscore_to_pvalue(double z);
int    te_kolmogorov_smirnov(const double *a, int na, const double *b, int nb, double *d_stat, double *p_value);
int    te_permutation_test(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_perms, double *p_value, double *null_dist);
void   te_bootstrap_ci(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_boot, double alpha, double *ci_lower, double *ci_upper, double *te_bootstrap);

/* ---- Multiple Comparison Corrections ---- */
void te_holm_correction(double *p_values, int n, double *corrected);
void te_benjamini_hochberg(double *p_values, int n, double alpha, int *significant);

#endif
