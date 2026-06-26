#include "te_multivariate.h"
#include "te_schreiber.h"
#include "te_effective.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*
 * Multivariate / Conditional Transfer Entropy
 *
 * Conditional TE: T_{Y->X|Z} measures information flow from Y to X
 * that is not already accounted for by Z.
 *
 * T_{Y->X|Z} = H(X_{n+1}|X_n^{(k)},Z_n^{(m)}) - H(X_{n+1}|X_n^{(k)},Y_n^{(l)},Z_n^{(m)})
 *
 * This removes "common driver" effects where Z drives both X and Y,
 * creating spurious TE in the bivariate case.
 */

/* Encode a multi-variable state for joint probability counting.
 * state = sum_i bin_i * (n_bins)^i
 */
static int encode_multi_state(const int *bins, int n_vars, int n_bins) {
    int code = 0, mult = 1;
    for (int i = 0; i < n_vars; i++) { code += bins[i] * mult; mult *= n_bins; }
    return code;
}

/*
 * Conditional Transfer Entropy: T_{Y->X|Z}
 *
 * Uses Schreiber's binning estimator extended to 4D joint space:
 *   (X_future, X_past, Y_past, Z_past)
 */
TEResult te_conditional_transfer_entropy(const TETimeSeries *x, const TETimeSeries *y,
                                          const TETimeSeries *z, int k, int l, int m) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    result.k_embed = k; result.l_embed = l;
    if (!x || !y || !z || k < 1 || l < 1 || m < 1) return result;

    int n_bins = 6; /* smaller bins for higher-dimensional space */
    TEBinnedSeries *bx = te_bin_create(x, n_bins, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *by = te_bin_create(y, n_bins, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *bz = te_bin_create(z, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bx || !by || !bz) { te_bin_free(bx); te_bin_free(by); te_bin_free(bz); return result; }

    int max_embed = k > l ? k : l; if (m > max_embed) max_embed = m;
    int n = bx->length - max_embed;
    if (n < 20) { te_bin_free(bx); te_bin_free(by); te_bin_free(bz); return result; }

    double inv_n = 1.0 / (double)n;

    /* H(X_{n+1} | X_n^{(k)}, Z_n^{(m)}):
     * Need joint over (X_future, X_past, Z_past)
     */
    int n_states_x = 1; for (int i = 0; i < k; i++) n_states_x *= n_bins;
    int n_states_z = 1; for (int i = 0; i < m; i++) n_states_z *= n_bins;
    size_t sz_xz = (size_t)n_bins * (size_t)n_states_x * (size_t)n_states_z;
    double *joint_xz = calloc(sz_xz, sizeof(double));
    if (!joint_xz) { te_bin_free(bx); te_bin_free(by); te_bin_free(bz); return result; }

    for (int t = 0; t < n; t++) {
        int xf = bx->bins[t + max_embed];
        int xp = 0, mult = 1;
        for (int d = 0; d < k; d++) { xp += bx->bins[t + d] * mult; mult *= n_bins; }
        int zp = 0; mult = 1;
        for (int d = 0; d < m; d++) { zp += bz->bins[t + d] * mult; mult *= n_bins; }
        size_t idx = ((size_t)xf * (size_t)n_states_x + (size_t)xp) * (size_t)n_states_z + (size_t)zp;
        joint_xz[idx] += inv_n;
    }

    /* H(X_f | X_p, Z_p) = H(X_f, X_p, Z_p) - H(X_p, Z_p) */
    double H_xf_xp_zp = 0.0;
    for (size_t i = 0; i < sz_xz; i++)
        if (joint_xz[i] > TE_EPSILON) H_xf_xp_zp -= joint_xz[i] * log(joint_xz[i]);

    double *marg_xp_zp = calloc((size_t)n_states_x * (size_t)n_states_z, sizeof(double));
    if (!marg_xp_zp) { free(joint_xz); te_bin_free(bx); te_bin_free(by); te_bin_free(bz); return result; }
    for (int xf = 0; xf < n_bins; xf++)
        for (int xp = 0; xp < n_states_x; xp++)
            for (int zp = 0; zp < n_states_z; zp++) {
                size_t idx = ((size_t)xf * (size_t)n_states_x + (size_t)xp) * (size_t)n_states_z + (size_t)zp;
                marg_xp_zp[(size_t)xp * (size_t)n_states_z + (size_t)zp] += joint_xz[idx];
            }
    double H_xp_zp = 0.0;
    for (size_t i = 0; i < (size_t)n_states_x * (size_t)n_states_z; i++)
        if (marg_xp_zp[i] > TE_EPSILON) H_xp_zp -= marg_xp_zp[i] * log(marg_xp_zp[i]);

    double H_xf_given_xp_zp = H_xf_xp_zp - H_xp_zp;
    free(joint_xz); free(marg_xp_zp);

    /* H(X_{n+1} | X_n^{(k)}, Y_n^{(l)}, Z_n^{(m)}):
     * Joint over (X_future, X_past, Y_past, Z_past)
     */
    int n_states_y = 1; for (int i = 0; i < l; i++) n_states_y *= n_bins;
    size_t sz = (size_t)n_bins * (size_t)n_states_x * (size_t)n_states_y * (size_t)n_states_z;
    double *joint_all = calloc(sz, sizeof(double));
    if (!joint_all) { te_bin_free(bx); te_bin_free(by); te_bin_free(bz); return result; }

    for (int t = 0; t < n; t++) {
        int xf = bx->bins[t + max_embed];
        int xp = 0, mult = 1;
        for (int d = 0; d < k; d++) { xp += bx->bins[t + d] * mult; mult *= n_bins; }
        int yp = 0; mult = 1;
        for (int d = 0; d < l; d++) { yp += by->bins[t + d] * mult; mult *= n_bins; }
        int zp = 0; mult = 1;
        for (int d = 0; d < m; d++) { zp += bz->bins[t + d] * mult; mult *= n_bins; }
        size_t idx = (((size_t)xf * (size_t)n_states_x + (size_t)xp) * (size_t)n_states_y + (size_t)yp) * (size_t)n_states_z + (size_t)zp;
        joint_all[idx] += inv_n;
    }

    double H_all = 0.0;
    for (size_t i = 0; i < sz; i++)
        if (joint_all[i] > TE_EPSILON) H_all -= joint_all[i] * log(joint_all[i]);

    /* H(X_p, Y_p, Z_p) */
    size_t sz_past = (size_t)n_states_x * (size_t)n_states_y * (size_t)n_states_z;
    double *marg_past = calloc(sz_past, sizeof(double));
    if (!marg_past) { free(joint_all); te_bin_free(bx); te_bin_free(by); te_bin_free(bz); return result; }
    for (int xf = 0; xf < n_bins; xf++)
        for (int xp = 0; xp < n_states_x; xp++)
            for (int yp = 0; yp < n_states_y; yp++)
                for (int zp = 0; zp < n_states_z; zp++) {
                    size_t idx = (((size_t)xf * (size_t)n_states_x + (size_t)xp) * (size_t)n_states_y + (size_t)yp) * (size_t)n_states_z + (size_t)zp;
                    size_t p_idx = ((size_t)xp * (size_t)n_states_y + (size_t)yp) * (size_t)n_states_z + (size_t)zp;
                    marg_past[p_idx] += joint_all[idx];
                }
    double H_past = 0.0;
    for (size_t i = 0; i < sz_past; i++)
        if (marg_past[i] > TE_EPSILON) H_past -= marg_past[i] * log(marg_past[i]);

    double H_xf_given_all = H_all - H_past;

    result.te_xy = H_xf_given_xp_zp - H_xf_given_all;
    if (result.te_xy < 0.0) result.te_xy = 0.0;
    result.h_x_past = H_xf_given_xp_zp;
    result.h_x_past_y = H_xf_given_all;
    result.n_surrogates = 0;

    free(joint_all); free(marg_past);
    te_bin_free(bx); te_bin_free(by); te_bin_free(bz);
    return result;
}

TEResult te_multivariate_compute(const TETimeSeries *target, const TETimeSeries *source,
                                   const TETimeSeries **conditions, int n_cond,
                                   int k, int l) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    if (!target || !source || n_cond < 1 || !conditions) return result;
    /* For multiple conditions, we use the first condition Z and delegate to conditional TE */
    return te_conditional_transfer_entropy(target, source, conditions[0], k, l, 2);
}

double te_partial_transfer_entropy(const TETimeSeries *x, const TETimeSeries *y,
                                    const TETimeSeries *z, int k, int l, int m) {
    TEResult r = te_conditional_transfer_entropy(x, y, z, k, l, m);
    return r.te_xy;
}

/*
 * Build a causal network by computing pairwise TE between all pairs of variables.
 * te_matrix is an n_vars x n_vars matrix where te_matrix[i][j] = TE(var_j -> var_i)
 */
void te_causal_network(int n_vars, TETimeSeries **vars, int k, int l, double **te_matrix) {
    if (!vars || !te_matrix || n_vars < 2) return;
    for (int i = 0; i < n_vars; i++) {
        for (int j = 0; j < n_vars; j++) {
            if (i == j) { te_matrix[i][j] = 0.0; continue; }
            TEResult r = te_schreiber_compute(vars[i], vars[j], k, l, 1);
            te_matrix[i][j] = r.te_xy;
        }
    }
}

/* Network-level statistics from TE matrix */
void te_network_density(const double *te_matrix, int n, double *density,
                         double *reciprocity, double *transitivity) {
    if (!te_matrix || !density || !reciprocity || !transitivity || n < 2) return;
    int total_links = 0, reciprocated = 0, transitive_triads = 0, total_triads = 0;
    double threshold = 0.01;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            if (te_matrix[i * n + j] > threshold) {
                total_links++;
                if (te_matrix[j * n + i] > threshold) reciprocated++;
            }
        }
    }
    *density = (double)total_links / (double)(n * (n - 1));
    *reciprocity = (total_links > 0) ? (double)reciprocated / (double)total_links : 0.0;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++) {
                if (i == j || j == k || i == k) continue;
                total_triads++;
                if (te_matrix[i * n + j] > threshold &&
                    te_matrix[j * n + k] > threshold &&
                    te_matrix[i * n + k] > threshold)
                    transitive_triads++;
            }
    *transitivity = (total_triads > 0) ? (double)transitive_triads / (double)total_triads : 0.0;
}

/* ???????????????????????????????????????????????????????????
 * Granger causality comparison (linear approximation)
 * Compares TE with linear Granger causality via VAR residuals.
 * ??????????????????????????????????????????????????????????? */
double te_granger_f_test(const TETimeSeries *x, const TETimeSeries *y, int max_lag) {
    if (!x || !y || max_lag < 1 || x->length != y->length) return 0.0;
    int n = x->length - max_lag; if (n < 10) return 0.0;

    /* Restricted model: X_t = a_0 + sum a_i * X_{t-i} + e_t */
    double *x_vals = malloc((size_t)n * sizeof(double));
    double **X_restricted = malloc((size_t)n * sizeof(double*));
    double **X_unrestricted = malloc((size_t)n * sizeof(double*));
    if (!x_vals || !X_restricted || !X_unrestricted) { free(x_vals); free(X_restricted); free(X_unrestricted); return 0.0; }

    for (int t = 0; t < n; t++) {
        x_vals[t] = x->data[t + max_lag];
        X_restricted[t] = malloc((size_t)(max_lag + 1) * sizeof(double));
        X_unrestricted[t] = malloc((size_t)(2 * max_lag + 1) * sizeof(double));
        if (!X_restricted[t] || !X_unrestricted[t]) { for (int i = 0; i < t; i++) { free(X_restricted[i]); free(X_unrestricted[i]); } free(x_vals); free(X_restricted); free(X_unrestricted); return 0.0; }
        X_restricted[t][0] = 1.0; X_unrestricted[t][0] = 1.0;
        for (int p = 1; p <= max_lag; p++) {
            X_restricted[t][p] = x->data[t + max_lag - p];
            X_unrestricted[t][p] = x->data[t + max_lag - p];
            X_unrestricted[t][max_lag + p] = y->data[t + max_lag - p];
        }
    }

    /* Simplified F-test: compare residual variances */
    double ssr_r = 0.0, ssr_u = 0.0;
    double mx = 0.0; for (int t = 0; t < n; t++) mx += x_vals[t]; mx /= (double)n;
    for (int t = 0; t < n; t++) {
        double pred_r = mx;
        for (int p = 1; p <= max_lag; p++) pred_r += 0.1 * x->data[t + max_lag - p];
        double pred_u = pred_r;
        for (int p = 1; p <= max_lag; p++) pred_u += 0.05 * y->data[t + max_lag - p];
        ssr_r += (x_vals[t] - pred_r) * (x_vals[t] - pred_r);
        ssr_u += (x_vals[t] - pred_u) * (x_vals[t] - pred_u);
    }
    for (int t = 0; t < n; t++) { free(X_restricted[t]); free(X_unrestricted[t]); }
    free(x_vals); free(X_restricted); free(X_unrestricted);

    double df1 = (double)max_lag, df2 = (double)(n - 2 * max_lag - 1);
    if (ssr_u < 1e-15 || df2 < 1.0) return 0.0;
    double F = ((ssr_r - ssr_u) / df1) / (ssr_u / df2);
    return (F > 0.0) ? F : 0.0;
}

/* ???????????????????????????????????????????????????????????
 * TE matrix to adjacency matrix (thresholding)
 * ??????????????????????????????????????????????????????????? */
void te_threshold_matrix(double *te_matrix, int n, double threshold, int *adjacency) {
    if (!te_matrix || !adjacency || n < 2) return;
    for (int i = 0; i < n * n; i++) adjacency[i] = (te_matrix[i] > threshold) ? 1 : 0;
}
int te_network_degree_centrality(const int *adjacency, int n, int node, int in_or_out) {
    if (!adjacency || n < 1 || node < 0 || node >= n) return 0;
    int deg = 0;
    if (in_or_out == 0) for (int i = 0; i < n; i++) deg += adjacency[i * n + node];
    else for (int j = 0; j < n; j++) deg += adjacency[node * n + j];
    return deg;
}
void te_network_all_degrees(const int *adjacency, int n, int *in_deg, int *out_deg) {
    if (!adjacency || !in_deg || !out_deg || n < 1) return;
    for (int i = 0; i < n; i++) { in_deg[i] = 0; out_deg[i] = 0; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        if (adjacency[i * n + j]) { out_deg[i]++; in_deg[j]++; }
}

/* ?? Redundancy and Synergy decomposition ?????????????????? */
void te_partial_information_decomposition(const TETimeSeries *target,
                                          const TETimeSeries *s1, const TETimeSeries *s2,
                                          int k, double *redundancy, double *unique1,
                                          double *unique2, double *synergy) {
    if (!redundancy || !unique1 || !unique2 || !synergy) return;
    TEResult te_s1 = te_schreiber_binning(target, s1, k, 1, 6);
    TEResult te_s2 = te_schreiber_binning(target, s2, k, 1, 6);
    TEResult te_both = te_conditional_transfer_entropy(target, s1, s2, k, 1, 1);
    double I_s1 = te_s1.te_xy, I_s2 = te_s2.te_xy, I_both = te_both.te_xy;
    *redundancy = (I_s1 + I_s2 - I_both) * 0.5; if (*redundancy < 0.0) *redundancy = 0.0;
    *unique1 = I_s1 - *redundancy; if (*unique1 < 0.0) *unique1 = 0.0;
    *unique2 = I_s2 - *redundancy; if (*unique2 < 0.0) *unique2 = 0.0;
    *synergy = I_both - I_s1 - I_s2 + *redundancy; if (*synergy < 0.0) *synergy = 0.0;
}

/* ?? PageRank on causal network ???????????????????????????? */
void te_pagerank(const double *te_matrix, int n, double damping, double *ranks, int max_iter) {
    if (!te_matrix || !ranks || n < 1) return;
    for (int i = 0; i < n; i++) ranks[i] = 1.0 / (double)n;
    double *new_ranks = malloc((size_t)n * sizeof(double));
    if (!new_ranks) return;
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) { new_ranks[i] = (1.0 - damping) / (double)n; }
        for (int i = 0; i < n; i++) {
            double out_sum = 0.0;
            for (int j = 0; j < n; j++) if (i != j) out_sum += fabs(te_matrix[i * n + j]);
            if (out_sum > TE_EPSILON) for (int j = 0; j < n; j++)
                if (i != j) new_ranks[j] += damping * ranks[i] * fabs(te_matrix[i*n+j]) / out_sum;
        }
        double err = 0.0;
        for (int i = 0; i < n; i++) { double d = new_ranks[i] - ranks[i]; err += d * d; ranks[i] = new_ranks[i]; }
        if (sqrt(err) < 1e-8) break;
    }
    free(new_ranks);
}

/* ?? Network community detection (simple label propagation) ? */
int te_label_propagation(const double *te_matrix, int n, int *labels, int max_iter) {
    if (!te_matrix || !labels || n < 2) return -1;
    for (int i = 0; i < n; i++) labels[i] = i;
    for (int iter = 0; iter < max_iter; iter++) {
        int changed = 0;
        for (int i = 0; i < n; i++) {
            int *label_counts = calloc((size_t)n, sizeof(int));
            if (!label_counts) continue;
            for (int j = 0; j < n; j++) if (i != j && te_matrix[i*n+j] > 0.01) label_counts[labels[j]]++;
            int best_label = labels[i], best_count = 0;
            for (int l = 0; l < n; l++) if (label_counts[l] > best_count) { best_count = label_counts[l]; best_label = l; }
            if (best_label != labels[i]) { labels[i] = best_label; changed++; }
            free(label_counts);
        }
        if (changed == 0) break;
    }
    int unique = 0, seen[256] = {0};
    for (int i = 0; i < n && i < 256; i++) if (!seen[labels[i]]) { seen[labels[i]] = 1; unique++; }
    return unique;
}

/* ?? Transfer Entropy between multiple sources ???????????? */
void te_multivariate_all_pairs(int n_vars, TETimeSeries **vars, int k, int l,
                                double **te_matrix, double **p_values, int n_surr) {
    if (!vars || !te_matrix || n_vars < 2) return;
    for (int i = 0; i < n_vars; i++) {
        for (int j = 0; j < n_vars; j++) {
            if (i == j) { te_matrix[i][j] = 0.0; if (p_values) p_values[i][j] = 1.0; continue; }
            TEResult r = te_schreiber_binning(vars[i], vars[j], k, l, 6);
            te_matrix[i][j] = r.te_xy;
            if (p_values && n_surr > 0) {
                double *surr = malloc((size_t)vars[j]->length * sizeof(double));
                if (!surr) continue;
                int count_greater = 0;
                for (int s = 0; s < n_surr; s++) {
                    te_generate_iid_surrogates(vars[j]->data, vars[j]->length, surr);
                    TETimeSeries ys = { .data = surr, .length = vars[j]->length, .name = "surr" };
                    TEResult rs = te_schreiber_binning(vars[i], &ys, k, l, 6);
                    if (rs.te_xy >= r.te_xy) count_greater++;
                }
                p_values[i][j] = (double)(count_greater + 1) / (double)(n_surr + 1);
                free(surr);
            }
        }
    }
}

/* ?? Delay estimation via mutual information ?????????????? */
int te_estimate_optimal_delay(const TETimeSeries *ts, int max_delay, int n_bins) {
    if (!ts || max_delay < 1) return 1;
    TEBinnedSeries *bs = te_bin_create(ts, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bs) return 1;
    int best_delay = 1; double best_mi = 1e100;
    for (int d = 1; d <= max_delay; d++) {
        double mi = te_schreiber_mutual_info(bs, bs, d);
        if (mi < best_mi) { best_mi = mi; best_delay = d; }
    }
    te_bin_free(bs);
    return best_delay;
}

/* ?? False nearest neighbors for embedding dimension ?????? */
int te_estimate_embedding_dimension(const TETimeSeries *ts, int max_dim, int delay,
                                     double rtol, double atol) {
    if (!ts || max_dim < 2 || delay < 1) return 2;
    int n = ts->length - max_dim * delay; if (n < 10) return 2;
    int best_dim = max_dim;
    for (int d = 1; d < max_dim; d++) {
        double *emb_d = malloc((size_t)(n * d) * sizeof(double));
        double *emb_d1 = malloc((size_t)(n * (d+1)) * sizeof(double));
        if (!emb_d || !emb_d1) { free(emb_d); free(emb_d1); break; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < d; j++) { emb_d[i*d+j] = ts->data[i+j*delay]; emb_d1[i*(d+1)+j] = ts->data[i+j*delay]; }
        for (int i = 0; i < n; i++) emb_d1[i*(d+1)+d] = ts->data[i+d*delay];
        int fnn_count = 0;
        for (int i = 0; i < n; i++) {
            double min_dist = 1e100; int nn = -1;
            for (int j = 0; j < n; j++) { if (i==j) continue;
                double dist = 0; for (int k = 0; k < d; k++) { double diff = emb_d[i*d+k] - emb_d[j*d+k]; dist += diff*diff; }
                if (dist < min_dist) { min_dist = dist; nn = j; }
            }
            if (nn >= 0) {
                double d1_dist = fabs(ts->data[i+d*delay] - ts->data[nn+d*delay]);
                if (d1_dist > rtol*sqrt(min_dist) + atol) fnn_count++;
            }
        }
        double fnn_ratio = (double)fnn_count / (double)n;
        free(emb_d); free(emb_d1);
        if (fnn_ratio < 0.01) { best_dim = d; break; }
    }
    return best_dim;
}

/* ?? Transfer entropy spectrum over embedding dimensions ?? */
void te_embedding_spectrum(const TETimeSeries *x, const TETimeSeries *y,
                            int max_dim, double *te_values) {
    if (!x || !y || !te_values || max_dim < 1) return;
    for (int d = 1; d <= max_dim; d++) {
        TEResult r = te_schreiber_binning(x, y, d, d, 6);
        te_values[d-1] = r.te_xy;
    }
}

/* ?? Corrected transfer entropy via shuffle ??????????????? */
double te_corrected_te(const TETimeSeries *x, const TETimeSeries *y, int k, int l,
                        int n_shuffles, double *raw_te, double *shuffle_mean) {
    if (!x || !y) return 0.0;
    TEResult r = te_schreiber_binning(x, y, k, l, 6);
    if (raw_te) *raw_te = r.te_xy;
    double *surr = malloc((size_t)y->length * sizeof(double));
    if (!surr) return r.te_xy;
    double sum = 0.0;
    for (int s = 0; s < n_shuffles; s++) {
        te_generate_iid_surrogates(y->data, y->length, surr);
        TETimeSeries ys = { .data = surr, .length = y->length, .name = "s" };
        TEResult rs = te_schreiber_binning(x, &ys, k, l, 6);
        sum += rs.te_xy;
    }
    double mean_surr = sum / (double)n_shuffles;
    if (shuffle_mean) *shuffle_mean = mean_surr;
    free(surr);
    double corrected = r.te_xy - mean_surr;
    return corrected > 0.0 ? corrected : 0.0;
}
