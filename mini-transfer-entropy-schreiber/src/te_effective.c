#include "te_effective.h"
#include "te_schreiber.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/*
 * Effective Transfer Entropy with bias correction.
 *
 * The naive TE estimator has a positive bias, especially for small samples.
 * Effective TE subtracts the bias estimated from shuffled surrogates:
 *
 *   ETE(Y->X) = TE(Y->X) - TE(Y_shuffled->X)
 *
 * Additionally, statistical significance is assessed via surrogate data
 * testing with multiple (typically 100-1000) random shuffles.
 */

/*
 * Fisher-Yates shuffle for generating IID surrogates.
 * Destroys temporal correlations while preserving the marginal distribution.
 */
void te_generate_iid_surrogates(const double *data, int n, double *surrogate) {
    if (!data || !surrogate || n < 1) return;
    memcpy(surrogate, data, (size_t)n * sizeof(double));
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        double t = surrogate[i];
        surrogate[i] = surrogate[j];
        surrogate[j] = t;
    }
}

/*
 * Fourier Transform (FT) surrogates:
 * Randomize phases in frequency domain while preserving power spectrum.
 * This preserves linear correlations (autocorrelation) while destroying
 * nonlinear / information-theoretic structure.
 */
void te_generate_ft_surrogates(const double *data, int n, double *surrogate) {
    if (!data || !surrogate || n < 4) return;
    int n2 = n / 2;
    double *real = malloc((size_t)n * sizeof(double));
    double *imag = calloc((size_t)n, sizeof(double));
    if (!real || !imag) { free(real); free(imag); return; }
    memcpy(real, data, (size_t)n * sizeof(double));

    /* Simple DFT */
    for (int k = 0; k < n; k++) {
        double re = 0.0, im = 0.0;
        for (int j = 0; j < n; j++) {
            double ang = -2.0 * TE_PI * (double)(k * j) / (double)n;
            re += real[j] * cos(ang);
            im += real[j] * sin(ang);
        }
        real[k] = re; imag[k] = im;
    }

    /* Randomize phases (preserve conjugate symmetry) */
    for (int k = 1; k < n2; k++) {
        double phase = 2.0 * TE_PI * (double)rand() / (double)RAND_MAX;
        double mag = sqrt(real[k] * real[k] + imag[k] * imag[k]);
        real[k] = mag * cos(phase);
        imag[k] = mag * sin(phase);
        real[n - k] = real[k];
        imag[n - k] = -imag[k];
    }
    /* DC and Nyquist remain real */
    imag[0] = 0.0;
    if (n % 2 == 0) imag[n2] = 0.0;

    /* Inverse DFT */
    for (int j = 0; j < n; j++) {
        double re = 0.0;
        for (int k = 0; k < n; k++) {
            double ang = 2.0 * TE_PI * (double)(k * j) / (double)n;
            re += real[k] * cos(ang) - imag[k] * sin(ang);
        }
        surrogate[j] = re / (double)n;
    }

    /* Amplitude adjust: map back to original sorted values to preserve distribution */
    double *orig_sorted = malloc((size_t)n * sizeof(double));
    double *surr_sorted = malloc((size_t)n * sizeof(double));
    int *ranks = malloc((size_t)n * sizeof(int));
    if (!orig_sorted || !surr_sorted || !ranks) {
        free(orig_sorted); free(surr_sorted); free(ranks);
        free(real); free(imag); return;
    }
    memcpy(orig_sorted, data, (size_t)n * sizeof(double));
    memcpy(surr_sorted, surrogate, (size_t)n * sizeof(double));
    /* Simple sort */
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++) {
            if (orig_sorted[i] > orig_sorted[j]) { double t = orig_sorted[i]; orig_sorted[i] = orig_sorted[j]; orig_sorted[j] = t; }
            if (surr_sorted[i] > surr_sorted[j]) { double t = surr_sorted[i]; surr_sorted[i] = surr_sorted[j]; surr_sorted[j] = t; }
        }
    /* Rank surrogate and map to original amplitudes */
    for (int i = 0; i < n; i++) {
        int rank = 0;
        for (int j = 0; j < n; j++) if (surrogate[j] < surrogate[i] || (surrogate[j] == surrogate[i] && j < i)) rank++;
        surrogate[i] = orig_sorted[rank];
    }

    free(orig_sorted); free(surr_sorted); free(ranks);
    free(real); free(imag);
}

/*
 * IAAFT (Iterative Amplitude Adjusted Fourier Transform) surrogates.
 * Iteratively corrects both power spectrum and amplitude distribution.
 */
void te_generate_iaaft_surrogates(const double *data, int n, double *surrogate, int max_iter) {
    if (!data || !surrogate || n < 4) return;
    double *orig_sorted = malloc((size_t)n * sizeof(double));
    if (!orig_sorted) return;
    memcpy(orig_sorted, data, (size_t)n * sizeof(double));
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if (orig_sorted[i] > orig_sorted[j]) { double t = orig_sorted[i]; orig_sorted[i] = orig_sorted[j]; orig_sorted[j] = t; }

    /* Start with a random shuffle */
    te_generate_iid_surrogates(data, n, surrogate);

    double *freq_re = malloc((size_t)n * sizeof(double));
    double *freq_im = calloc((size_t)n, sizeof(double));
    if (!freq_re || !freq_im) { free(orig_sorted); free(freq_re); free(freq_im); return; }

    for (int iter = 0; iter < max_iter; iter++) {
        /* Step 1: Adjust to match Fourier amplitudes of original data */
        for (int k = 0; k < n; k++) {
            double re = 0.0, im = 0.0;
            for (int j = 0; j < n; j++) {
                double ang = -2.0 * TE_PI * (double)(k * j) / (double)n;
                re += surrogate[j] * cos(ang);
                im += surrogate[j] * sin(ang);
            }
            freq_re[k] = re; freq_im[k] = im;
        }
        /* Replace with original amplitudes */
        double *orig_re = malloc((size_t)n * sizeof(double));
        double *orig_im = calloc((size_t)n, sizeof(double));
        if (!orig_re || !orig_im) { free(orig_re); free(orig_im); continue; }
        for (int j = 0; j < n; j++) {
            double re2 = 0.0;
            for (int k = 0; k < n; k++) {
                double ang = -2.0 * TE_PI * (double)(k * j) / (double)n;
                re2 += data[k] * cos(ang);
            }
            orig_re[j] = re2;
        }
        /* Keep original magnitude, surrogate phase */
        for (int k = 0; k < n; k++) {
            double mag_o = sqrt(orig_re[k] * orig_re[k] + orig_im[k] * orig_im[k]);
            double mag_s = sqrt(freq_re[k] * freq_re[k] + freq_im[k] * freq_im[k]);
            if (mag_s > TE_EPSILON) {
                double scale = mag_o / mag_s;
                freq_re[k] *= scale; freq_im[k] *= scale;
            }
        }
        free(orig_re); free(orig_im);

        /* Inverse DFT */
        for (int j = 0; j < n; j++) {
            double re = 0.0;
            for (int k = 0; k < n; k++) {
                double ang = 2.0 * TE_PI * (double)(k * j) / (double)n;
                re += freq_re[k] * cos(ang) - freq_im[k] * sin(ang);
            }
            surrogate[j] = re / (double)n;
        }

        /* Step 2: Amplitude adjust (rank-order remapping) */
        double *surr_sorted = malloc((size_t)n * sizeof(double));
        if (!surr_sorted) continue;
        memcpy(surr_sorted, surrogate, (size_t)n * sizeof(double));
        for (int i = 0; i < n - 1; i++)
            for (int j = i + 1; j < n; j++)
                if (surr_sorted[i] > surr_sorted[j]) { double t = surr_sorted[i]; surr_sorted[i] = surr_sorted[j]; surr_sorted[j] = t; }
        for (int i = 0; i < n; i++) {
            int rank = 0;
            for (int j = 0; j < n; j++) if (surrogate[j] < surrogate[i] || (surrogate[j] == surrogate[i] && j < i)) rank++;
            surrogate[i] = orig_sorted[rank];
        }
        free(surr_sorted);
    }
    free(orig_sorted); free(freq_re); free(freq_im);
}

/* Estimate shuffle bias by averaging TE over n_shuffles random permutations */
double te_effective_shuffle_bias(const TETimeSeries *x, const TETimeSeries *y,
                                   int k, int l, int n_shuffles) {
    if (!x || !y || n_shuffles < 1) return 0.0;
    double *y_shuffled = malloc((size_t)y->length * sizeof(double));
    if (!y_shuffled) return 0.0;
    /* Fill initial data to avoid uninitialized warning */
    memcpy(y_shuffled, y->data, (size_t)y->length * sizeof(double));
    TETimeSeries *ys = te_ts_create(y_shuffled, y->length, "shuffled");
    if (!ys) { free(y_shuffled); return 0.0; }

    double te_sum = 0.0;
    for (int s = 0; s < n_shuffles; s++) {
        te_generate_iid_surrogates(y->data, y->length, y_shuffled);
        memcpy(ys->data, y_shuffled, (size_t)y->length * sizeof(double));
        TEResult r = te_schreiber_binning(x, ys, k, l, 8);
        te_sum += r.te_xy;
    }
    te_ts_free(ys);
    free(y_shuffled);
    return te_sum / (double)n_shuffles;
}

/*
 * Effective Transfer Entropy with surrogate-based bias correction.
 */
TEResult te_effective_compute(const TETimeSeries *x, const TETimeSeries *y,
                                int k, int l, int n_surrogates) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    result.k_embed = k; result.l_embed = l;
    result.n_surrogates = n_surrogates;
    if (!x || !y || k < 1 || l < 1) return result;

    /* Raw TE */
    TEResult raw = te_schreiber_binning(x, y, k, l, 8);
    result.te_xy = raw.te_xy;
    result.te_yx = raw.te_yx;

    /* Surrogate distribution */
    double *surr_te = malloc((size_t)n_surrogates * sizeof(double));
    double *y_shuffled = malloc((size_t)y->length * sizeof(double));
    if (!surr_te || !y_shuffled) { free(surr_te); free(y_shuffled); return result; }

    TETimeSeries *ys = te_ts_create(y_shuffled, y->length, "surr");
    if (!ys) { free(surr_te); free(y_shuffled); return result; }

    double te_mean = 0.0, te_std = 0.0;
    for (int s = 0; s < n_surrogates; s++) {
        te_generate_iid_surrogates(y->data, y->length, y_shuffled);
        memcpy(ys->data, y_shuffled, (size_t)y->length * sizeof(double));
        TEResult r = te_schreiber_binning(x, ys, k, l, 8);
        surr_te[s] = r.te_xy;
        te_mean += surr_te[s];
    }
    te_mean /= (double)n_surrogates;

    for (int s = 0; s < n_surrogates; s++) {
        double d = surr_te[s] - te_mean;
        te_std += d * d;
    }
    te_std = sqrt(te_std / (double)(n_surrogates - 1));

    /* Effective TE = raw - mean(surrogate) */
    double eff_te = result.te_xy - te_mean;
    if (eff_te < 0.0) eff_te = 0.0;

    /* Z-score */
    if (te_std > TE_EPSILON) {
        result.z_score = (result.te_xy - te_mean) / te_std;
        result.p_value = te_zscore_to_pvalue(result.z_score);
    }
    result.significant = (result.p_value < 0.05);

    /* Store effective TE in the result */
    result.te_xy = eff_te;
    result.net_te = eff_te - result.te_yx;

    te_ts_free(ys);
    free(surr_te); free(y_shuffled);
    return result;
}

TEResult te_effective_bias_corrected(const TETimeSeries *x, const TETimeSeries *y,
                                       int k, int l, int n_shuffles) {
    TEResult result = te_schreiber_binning(x, y, k, l, 8);
    if (n_shuffles > 0) {
        double bias = te_effective_shuffle_bias(x, y, k, l, n_shuffles);
        result.te_xy -= bias;
        if (result.te_xy < 0.0) result.te_xy = 0.0;
        result.net_te = result.te_xy - result.te_yx;
    }
    return result;
}

/* Z-score to p-value using error function approximation */
double te_zscore_to_pvalue(double z) {
    /* Two-tailed test */
    double abs_z = fabs(z);
    double x = abs_z / sqrt(2.0);
    /* ERF approximation (Abramowitz & Stegun 7.1.26) */
    double t = 1.0 / (1.0 + 0.3275911 * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    double a4 = -1.453152027, a5 = 1.061405429;
    double erf = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * exp(-x * x);
    return 1.0 - erf;
}

/*
 * Kolmogorov-Smirnov two-sample test.
 * Tests whether two samples come from the same distribution.
 */
int te_kolmogorov_smirnov(const double *a, int na, const double *b, int nb,
                           double *d_stat, double *p_value) {
    if (!a || !b || na < 2 || nb < 2 || !d_stat || !p_value) return -1;
    double *combined = malloc((size_t)(na + nb) * sizeof(double));
    if (!combined) return -1;
    memcpy(combined, a, (size_t)na * sizeof(double));
    memcpy(combined + na, b, (size_t)nb * sizeof(double));
    for (int i = 0; i < na + nb - 1; i++)
        for (int j = i + 1; j < na + nb; j++)
            if (combined[i] > combined[j]) { double t = combined[i]; combined[i] = combined[j]; combined[j] = t; }

    double max_diff = 0.0;
    double cdf_a = 0.0, cdf_b = 0.0;
    int ai = 0, bi = 0;
    for (int k = 0; k < na + nb; k++) {
        while (ai < na && a[ai] <= combined[k]) { cdf_a += 1.0 / (double)na; ai++; }
        while (bi < nb && b[bi] <= combined[k]) { cdf_b += 1.0 / (double)nb; bi++; }
        double diff = fabs(cdf_a - cdf_b);
        if (diff > max_diff) max_diff = diff;
    }
    *d_stat = max_diff;
    double lambda = sqrt((double)na * (double)nb / (double)(na + nb)) * max_diff;
    *p_value = 2.0 * exp(-2.0 * lambda * lambda);
    free(combined);
    return 0;
}

/* ???????????????????????????????????????????????????????????
 * Block Bootstrap for dependent time series
 * Resamples blocks of length L to preserve temporal structure.
 * ??????????????????????????????????????????????????????????? */
void te_generate_block_bootstrap(const double *data, int n, int block_len, double *surrogate) {
    if (!data || !surrogate || n < 1 || block_len < 1) return;
    for (int pos = 0; pos < n; ) {
        int block_start = rand() % (n - block_len + 1);
        for (int j = 0; j < block_len && pos < n; j++, pos++)
            surrogate[pos] = data[block_start + j];
    }
}

/* ???????????????????????????????????????????????????????????
 * Stationary Bootstrap (Politis & Romano, 1994)
 * Random block lengths with geometric distribution.
 * ??????????????????????????????????????????????????????????? */
void te_generate_stationary_bootstrap(const double *data, int n, double p, double *surrogate) {
    if (!data || !surrogate || n < 1 || p <= 0.0 || p > 1.0) return;
    int pos = 0;
    while (pos < n) {
        int start = rand() % n;
        int block_len = 1;
        while ((double)rand() / (double)RAND_MAX > p) block_len++;
        for (int j = 0; j < block_len && pos < n; j++, pos++)
            surrogate[pos] = data[(start + j) % n];
    }
}

/* ???????????????????????????????????????????????????????????
 * Circular shift surrogate (preserves autocorrelation)
 * ??????????????????????????????????????????????????????????? */
void te_generate_circular_shift(const double *data, int n, double *surrogate) {
    if (!data || !surrogate || n < 2) return;
    int shift = rand() % (n - 1) + 1;
    memcpy(surrogate, data + shift, (size_t)(n - shift) * sizeof(double));
    memcpy(surrogate + n - shift, data, (size_t)shift * sizeof(double));
}

/* ???????????????????????????????????????????????????????????
 * Multiple comparison correction (Bonferroni-Holm)
 * Corrects p-values for multiple TE tests in a network.
 * ??????????????????????????????????????????????????????????? */
void te_holm_correction(double *p_values, int n, double *corrected) {
    if (!p_values || !corrected || n < 1) return;
    int *order = malloc((size_t)n * sizeof(int));
    if (!order) return;
    for (int i = 0; i < n; i++) order[i] = i;
    for (int i = 0; i < n - 1; i++) for (int j = i + 1; j < n; j++)
        if (p_values[order[i]] > p_values[order[j]]) { int t = order[i]; order[i] = order[j]; order[j] = t; }
    for (int i = 0; i < n; i++) {
        double holm = p_values[order[i]] * (double)(n - i);
        if (holm > 1.0) holm = 1.0;
        corrected[order[i]] = holm;
    }
    free(order);
}

/* ?? Permutation test for TE significance ????????????????? */
int te_permutation_test(const TETimeSeries *x, const TETimeSeries *y,
                         int k, int l, int n_perms, double *p_value, double *null_dist) {
    if (!x || !y || !p_value || n_perms < 10) return -1;
    TEResult obs = te_schreiber_binning(x, y, k, l, 6);
    double *surr = malloc((size_t)y->length * sizeof(double));
    if (!surr) return -1;
    if (null_dist) null_dist[0] = obs.te_xy;
    int count_greater = 0;
    for (int p = 0; p < n_perms; p++) {
        te_generate_iid_surrogates(y->data, y->length, surr);
        TETimeSeries ys = { .data = surr, .length = y->length, .name = "perm" };
        TEResult r = te_schreiber_binning(x, &ys, k, l, 6);
        if (null_dist) null_dist[p+1] = r.te_xy;
        if (r.te_xy >= obs.te_xy) count_greater++;
    }
    *p_value = (double)(count_greater + 1) / (double)(n_perms + 1);
    free(surr);
    return 0;
}

/* ?? Bootstrap confidence interval for TE ?????????????????? */
void te_bootstrap_ci(const TETimeSeries *x, const TETimeSeries *y,
                      int k, int l, int n_boot, double alpha,
                      double *ci_lower, double *ci_upper, double *te_bootstrap) {
    if (!x || !y || !ci_lower || !ci_upper || n_boot < 10) return;
    double *te_vals = malloc((size_t)n_boot * sizeof(double));
    if (!te_vals) return;
    for (int b = 0; b < n_boot; b++) {
        int n = x->length;
        double *x_boot = malloc((size_t)n * sizeof(double));
        double *y_boot = malloc((size_t)n * sizeof(double));
        if (!x_boot || !y_boot) { free(x_boot); free(y_boot); continue; }
        for (int i = 0; i < n; i++) { int idx = rand() % n; x_boot[i]=x->data[idx]; y_boot[i]=y->data[idx]; }
        TETimeSeries xs={.data=x_boot,.length=n,.name="xb"}, ys={.data=y_boot,.length=n,.name="yb"};
        TEResult r = te_schreiber_binning(&xs, &ys, k, l, 6);
        te_vals[b] = r.te_xy; free(x_boot); free(y_boot);
    }
    for (int i=0;i<n_boot-1;i++) for (int j=i+1;j<n_boot;j++) if (te_vals[i]>te_vals[j]){double t=te_vals[i];te_vals[i]=te_vals[j];te_vals[j]=t;}
    int lo_idx = (int)(alpha/2.0 * (double)n_boot), hi_idx = (int)((1.0-alpha/2.0)*(double)n_boot);
    if (lo_idx<0) lo_idx=0; if (hi_idx>=n_boot) hi_idx=n_boot-1;
    *ci_lower=te_vals[lo_idx]; *ci_upper=te_vals[hi_idx];
    if (te_bootstrap) memcpy(te_bootstrap, te_vals, (size_t)n_boot*sizeof(double));
    free(te_vals);
}

/* ?? Jackknife bias correction ???????????????????????????? */
double te_jackknife_te(const TETimeSeries *x, const TETimeSeries *y, int k, int l) {
    if (!x || !y || x->length < 10) return 0.0;
    int n = x->length;
    TEResult full = te_schreiber_binning(x, y, k, l, 6);
    double sum = 0.0; int count = 0;
    int block = n / 20; if (block < 1) block = 1;
    for (int start = 0; start < n; start += block) {
        int len = n - block;
        double *xj = malloc((size_t)len*sizeof(double)), *yj = malloc((size_t)len*sizeof(double));
        if (!xj || !yj) { free(xj); free(yj); continue; }
        int pos = 0;
        for (int i=0;i<start;i++) { xj[pos]=x->data[i]; yj[pos]=y->data[i]; pos++; }
        for (int i=start+block;i<n;i++) { xj[pos]=x->data[i]; yj[pos]=y->data[i]; pos++; }
        TETimeSeries xs={.data=xj,.length=len,.name="xj"}, ys={.data=yj,.length=len,.name="yj"};
        TEResult r = te_schreiber_binning(&xs, &ys, k, l, 6);
        sum += r.te_xy; count++; free(xj); free(yj);
    }
    if (count < 2) return full.te_xy;
    double mean_jk = sum / (double)count;
    double bias = ((double)(count-1)) * (mean_jk - full.te_xy);
    return full.te_xy - bias;
}

/* ?? False Discovery Rate (Benjamini-Hochberg) ???????????? */
void te_benjamini_hochberg(double *p_values, int n, double alpha, int *significant) {
    if (!p_values || !significant || n < 1) return;
    int *order = malloc((size_t)n * sizeof(int));
    if (!order) return;
    for (int i=0;i<n;i++) order[i]=i;
    for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) if (p_values[order[i]]>p_values[order[j]]){int t=order[i];order[i]=order[j];order[j]=t;}
    for (int i=0;i<n;i++) significant[i]=0;
    for (int i=n-1;i>=0;i--) {
        double thresh = alpha * (double)(i+1) / (double)n;
        if (p_values[order[i]] <= thresh) { for (int j=0;j<=i;j++) significant[order[j]]=1; break; }
    }
    free(order);
}
