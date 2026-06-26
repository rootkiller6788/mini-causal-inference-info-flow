#ifndef TE_CORE_H
#define TE_CORE_H
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#define TE_EPSILON  1e-12
#define TE_PI       3.14159265358979323846
#define TE_E        2.71828182845904523536
#define TE_MAX_BINS 256
#define TE_MAX_DIM  32
#define TE_MAX_EMBED 16

/* ================================================================
 * Core Data Types
 * ================================================================ */

typedef struct { double *data; int length; char *name; } TETimeSeries;

/* ---- Time Series Create/Free/Transform ---- */
TETimeSeries* te_ts_create(const double *data, int len, const char *name);
void         te_ts_free(TETimeSeries *ts);
TETimeSeries* te_ts_clone(const TETimeSeries *ts);
void         te_ts_normalize(TETimeSeries *ts);
void         te_ts_standardize(TETimeSeries *ts);
void         te_ts_diff(TETimeSeries *ts);
void         te_ts_shuffle(TETimeSeries *ts);
void         te_ts_resample(TETimeSeries *ts, int new_length);
int          te_ts_subsample(const TETimeSeries *ts, int step, TETimeSeries **out);

/* ---- Time Series Statistics ---- */
double te_ts_mean(const TETimeSeries *ts);
double te_ts_variance(const TETimeSeries *ts);
double te_ts_std(const TETimeSeries *ts);
double te_ts_min(const TETimeSeries *ts);
double te_ts_max(const TETimeSeries *ts);
double te_ts_skewness(const TETimeSeries *ts);
double te_ts_kurtosis(const TETimeSeries *ts);
double te_ts_rms(const TETimeSeries *ts);
double te_ts_quantile(const TETimeSeries *ts, double q);
double te_ts_percentile(const TETimeSeries *ts, double pct);
double te_ts_empirical_cdf(const TETimeSeries *ts, double x);
int    te_ts_zero_crossings(const TETimeSeries *ts);
void   te_ts_summary(const TETimeSeries *ts, char *buf, int buf_size);
void   te_ts_plot_ascii(const TETimeSeries *ts, int width, int height);

/* ---- Time Series Signal Processing ---- */
void te_ts_detrend_linear(TETimeSeries *ts);
void te_ts_moving_average(TETimeSeries *ts, int window);
void te_ts_exponential_smooth(TETimeSeries *ts, double alpha);
void te_ts_add_noise(TETimeSeries *ts, double sigma);
void te_ts_hamming_window(TETimeSeries *ts);
void te_ts_hann_window(TETimeSeries *ts);
void te_ts_blackman_window(TETimeSeries *ts);
void te_ts_welch_psd(const TETimeSeries *ts, int nfft, int overlap, double *psd, double *freqs);

/* ---- Autocorrelation & Cross-Correlation ---- */
double te_ts_autocorrelation(const TETimeSeries *ts, int lag);
void   te_ts_autocorrelation_function(const TETimeSeries *ts, int max_lag, double *acf);
double te_ts_cross_correlation(const TETimeSeries *x, const TETimeSeries *y, int lag);
double te_ts_durbin_watson(const TETimeSeries *ts);

/* ---- Time Series Stationarity Tests ---- */
double te_ts_adf_statistic(const TETimeSeries *ts, int max_lag);
double te_ts_kpss_stat(const TETimeSeries *ts);
double te_ts_hurst_rs(const TETimeSeries *ts, int min_window, int max_window, int n_windows);

/* ---- Time Series Generation ---- */
TETimeSeries* te_ts_generate_sine(int n, double freq, double amp, double phase);
TETimeSeries* te_ts_generate_henon(int n, double a, double b);
TETimeSeries* te_ts_generate_logistic(int n, double r, double x0);
TETimeSeries* te_ts_generate_ar1(int n, double phi, double sigma);
TETimeSeries* te_ts_generate_coupled_ar1(int n, double c, double phi, double sigma);

/* ---- Time Series Window Analysis ---- */
void   te_ts_sliding_stats(const TETimeSeries *ts, int window, double *mean_out, double *var_out, double *skew_out);
void   te_ts_segment_stats(const TETimeSeries *ts, int n_segments, double *means, double *vars, double *mins, double *maxs);
TETimeSeries* te_ts_aggregate(const TETimeSeries *ts, int factor, int mode);
int    te_ts_histogram(const TETimeSeries *ts, int n_bins, double *edges, double *counts);
void   te_ts_poisson_disc_sample(const TETimeSeries *ts, double min_dist, int *indices, int *n_indices, int max_indices);

/* ---- Hjorth Parameters (Activity / Mobility / Complexity) ---- */
void te_ts_hjorth(const TETimeSeries *ts, double *activity, double *mobility, double *complexity);

/* ---- Phase Randomization Surrogate ---- */
void te_ts_phase_randomize(const TETimeSeries *ts, double *surrogate);

/* ---- Q-Q Plot Data ---- */
int te_ts_qq_data(const TETimeSeries *a, const TETimeSeries *b, double *qa, double *qb, int n_quantiles);

/* ---- Time Series I/O ---- */
int te_ts_load_csv(const char *filename, TETimeSeries **ts_out, int *n_series);
int te_ts_save_csv(const TETimeSeries *ts, const char *filename);
int te_ts_save_binary(const TETimeSeries *ts, const char *filename);
int te_ts_load_binary(const char *filename, TETimeSeries **ts);

/* ================================================================
 * Binning
 * ================================================================ */
typedef enum { TE_BIN_EQUAL_WIDTH=0, TE_BIN_EQUAL_FREQ=1, TE_BIN_KD_TREE=2, TE_BIN_SYMBOLIC=3 } TEBinningMethod;
typedef struct { int *bins; int n_bins; int length; double *edges; TEBinningMethod method; } TEBinnedSeries;

TEBinnedSeries* te_bin_create(const TETimeSeries *ts, int n_bins, TEBinningMethod method);
TEBinnedSeries* te_bin_create_adaptive(const TETimeSeries *ts, int min_bins, int max_bins);
void   te_bin_free(TEBinnedSeries *bs);
int    te_bin_get(const TEBinnedSeries *bs, int idx);
void   te_bin_count(const TEBinnedSeries *bs, double *hist);
double te_bin_entropy(const TEBinnedSeries *bs);
int    te_bin_index_2d(const TEBinnedSeries *bx, const TEBinnedSeries *by, int idx);
int    te_estimate_optimal_bins(const TETimeSeries *ts);

/* ================================================================
 * Embedding
 * ================================================================ */
typedef struct { int **states; int n_states; int embed_dim; int delay; } TEEmbedding;
TEEmbedding* te_embed_create(const TEBinnedSeries *bs, int dim, int delay);
void te_embed_free(TEEmbedding *emb);

/* ================================================================
 * Joint Probability & Entropy
 * ================================================================ */
typedef struct { double *joint; double *marg_x; double *marg_y; int nx, ny; int n_samples; } TEJointProb;
TEJointProb* te_joint_estimate(const TEBinnedSeries *x, const TEBinnedSeries *y, int x_lag);
void   te_joint_free(TEJointProb *jp);
double te_joint_entropy(const TEJointProb *jp);
double te_conditional_entropy(const TEJointProb *jp);
double te_joint_entropy_2d(const TEBinnedSeries *bx, const TEBinnedSeries *by);
double te_conditional_entropy_2d(const TEBinnedSeries *bx, const TEBinnedSeries *by);

/* ================================================================
 * Transfer Entropy Result
 * ================================================================ */
typedef struct { double te_xy; double te_yx; double net_te; double h_x_future; double h_x_past; double h_x_past_y; int k_embed; int l_embed; double z_score; double p_value; int n_surrogates; bool significant; } TEResult;
void te_result_print(const TEResult *r);
int  te_result_is_significant(const TEResult *r, double alpha);
void te_result_to_string(const TEResult *r, char *buf, int buf_size);

/* ================================================================
 * Mutual Information & Information-Theoretic Measures
 * ================================================================ */
double te_mutual_information(const TETimeSeries *x, const TETimeSeries *y, int n_bins, int lag);
double te_active_information_storage(const TETimeSeries *x, int k, int n_bins);
double te_predictive_information(const TETimeSeries *x, int k, int h, int n_bins);

/* ---- Generalized Entropies ---- */
double te_renyi_entropy(const TEBinnedSeries *bs, double alpha);
double te_tsallis_entropy(const TEBinnedSeries *bs, double q);
double te_js_divergence(const double *p, const double *q, int n);

/* ---- Complexity Measures ---- */
double te_lz_complexity(const TEBinnedSeries *bs);
int    te_lz76_complexity(const int *seq, int n);
double te_effective_measure_complexity(const TETimeSeries *ts, int max_k, int n_bins);
double te_excess_entropy(const TETimeSeries *ts, int max_k, int n_bins);
double te_approximate_entropy(const TETimeSeries *ts, int m, double r);
double te_sample_entropy(const TETimeSeries *ts, int m, double r);

/* ---- Kernel Density Estimation ---- */
double te_kde_entropy(const double *data, int n, double bw);
double te_kde_mutual_information(const double *x, const double *y, int n, double bw_x, double bw_y);

/* ================================================================
 * Statistical Tests
 * ================================================================ */
double te_mann_whitney_u(const double *a, int na, const double *b, int nb, double *p_value);
double te_bds_statistic(const TETimeSeries *ts, int embed_dim, double epsilon);

#endif
