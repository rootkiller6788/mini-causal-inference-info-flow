#ifndef TE_MULTIVARIATE_H
#define TE_MULTIVARIATE_H
#include "te_core.h"

/* ================================================================
 * Multivariate / Conditional Transfer Entropy
 * Extends bivariate TE to handle confounders (Z) and networks.
 * ================================================================ */

/* ---- Conditional & Multivariate TE ---- */
TEResult te_multivariate_compute(const TETimeSeries *target, const TETimeSeries *source, const TETimeSeries **conditions, int n_cond, int k, int l);
TEResult te_conditional_transfer_entropy(const TETimeSeries *x, const TETimeSeries *y, const TETimeSeries *z, int k, int l, int m);
double   te_partial_transfer_entropy(const TETimeSeries *x, const TETimeSeries *y, const TETimeSeries *z, int k, int l, int m);
double   te_corrected_te(const TETimeSeries *x, const TETimeSeries *y, int k, int l, int n_shuffles, double *raw_te, double *shuffle_mean);
void     te_multivariate_all_pairs(int n_vars, TETimeSeries **vars, int k, int l, double **te_matrix, double **p_values, int n_surr);

/* ---- Causal Network Analysis ---- */
void te_causal_network(int n_vars, TETimeSeries **vars, int k, int l, double **te_matrix);
void te_network_density(const double *te_matrix, int n, double *density, double *reciprocity, double *transitivity);
void te_threshold_matrix(double *te_matrix, int n, double threshold, int *adjacency);
int  te_network_degree_centrality(const int *adjacency, int n, int node, int in_or_out);
void te_network_all_degrees(const int *adjacency, int n, int *in_deg, int *out_deg);

/* ---- Network Graph Algorithms ---- */
void te_pagerank(const double *te_matrix, int n, double damping, double *ranks, int max_iter);
int  te_label_propagation(const double *te_matrix, int n, int *labels, int max_iter);

/* ---- Partial Information Decomposition (Williams & Beer, 2010) ---- */
void te_partial_information_decomposition(const TETimeSeries *target, const TETimeSeries *s1, const TETimeSeries *s2, int k, double *redundant, double *unique1, double *unique2, double *synergy);

/* ---- Granger Causality (linear alternative) ---- */
double te_granger_f_test(const TETimeSeries *x, const TETimeSeries *y, int max_lag);

/* ---- Embedding Selection ---- */
int  te_estimate_optimal_delay(const TETimeSeries *ts, int max_delay, int n_bins);
int  te_estimate_embedding_dimension(const TETimeSeries *ts, int max_dim, int delay, double rtol, double atol);
void te_embedding_spectrum(const TETimeSeries *x, const TETimeSeries *y, int max_dim, double *te_values);

#endif
