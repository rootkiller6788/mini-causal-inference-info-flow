#include "te_symbolic.h"
#include "te_effective.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*
 * Symbolic Transfer Entropy (Staniek & Lehnertz, PRL 2008)
 *
 * Uses ordinal pattern (permutation) encoding instead of binning.
 * For an embedding of order m with delay tau, each window of m values
 * is mapped to one of m! possible ordinal patterns based on the
 * rank ordering of the values.
 *
 * Advantages:
 *   - Parameter-free (no bin width selection)
 *   - Robust to monotonic nonlinear transformations
 *   - Computationally efficient
 */

/* Encode a window of 'order' values into a permutation pattern.
 * pattern[0..order-1] gets the ranks (0 = smallest, order-1 = largest).
 * Returns the Lehmer code (lexicographic index 0..order!-1).
 */
int te_symbolic_permutation_encode(const double *window, int order, int *pattern) {
    if (!window || order < 2 || order > 8) return -1;
    /* Compute the ordinal pattern: for each position i, pattern[i] = number of
     * elements j < i such that window[j] <= window[i].
     * This gives the inversion vector (Lehmer code) of the permutation. */
    for (int i = 0; i < order; i++) {
        pattern[i] = 0;
        for (int j = 0; j < order; j++) {
            if (window[j] < window[i]) pattern[i]++;
            /* Tie-breaking: use index order for equal values */
            if (fabs(window[j] - window[i]) < TE_EPSILON && j < i) pattern[i]++;
        }
    }
    /* Convert pattern to integer index via factorial number system.
     * The pattern[i] is the rank among ALL elements (including those
     * already accounted for). We need the rank among REMAINING elements.
     *
     * Correct Lehmer code: for each i, count how many UNUSED elements
     * are smaller than the element at position i. Then multiply by
     * (order-1-i)! and sum.
     *
     * Algorithm: scan positions 0..order-1, finding the rank-th unused
     * element. The number of smaller unused elements is exactly
     * pattern[i] minus the number of already-used smaller elements.
     */
    int index = 0;
    int used_positions[8] = {0};  /* which original positions have been consumed */
    /* Precompute factorials */
    int fact[8];
    fact[0] = 1;
    for (int i = 1; i < order; i++) fact[i] = fact[i-1] * i;
    for (int i = 0; i < order; i++) {
        int rank = pattern[i];
        /* Count how many SMALLER elements are still unused */
        int smaller_unused = 0;
        for (int j = 0; j < order; j++) {
            if (pattern[j] < rank && !used_positions[j]) smaller_unused++;
        }
        index += smaller_unused * fact[order - 1 - i];
        /* Mark the current element as used */
        used_positions[i] = 1;
    }
    return index;
}

/* Count occurrences of each ordinal pattern in a time series.
 * n_patterns should be order! (factorial of order).
 */
void te_symbolic_pattern_count(const TETimeSeries *ts, int order, int delay,
                                int *counts, int n_patterns) {
    if (!ts || !counts || order < 2 || delay < 1) return;
    int n = ts->length - (order - 1) * delay;
    if (n < 1) return;
    memset(counts, 0, (size_t)n_patterns * sizeof(int));
    double window[8];
    int pattern[8];
    for (int t = 0; t < n; t++) {
        for (int d = 0; d < order; d++) window[d] = ts->data[t + d * delay];
        int idx = te_symbolic_permutation_encode(window, order, pattern);
        if (idx >= 0 && idx < n_patterns) counts[idx]++;
    }
}

/* Permutation entropy: H = -sum p(pi) * log(p(pi))
 * Normalized: H_norm = H / log(order!)
 */
double te_symbolic_permutation_entropy(const TETimeSeries *ts, int order, int delay) {
    if (!ts || order < 2 || order > 8 || delay < 1) return 0.0;
    int n_patterns = 1;
    for (int i = 2; i <= order; i++) n_patterns *= i;
    int *counts = calloc((size_t)n_patterns, sizeof(int));
    if (!counts) return 0.0;
    te_symbolic_pattern_count(ts, order, delay, counts, n_patterns);
    int n = ts->length - (order - 1) * delay;
    double H = 0.0, inv_n = 1.0 / (double)n;
    for (int i = 0; i < n_patterns; i++) {
        if (counts[i] > 0) { double p = (double)counts[i] * inv_n; H -= p * log(p); }
    }
    free(counts);
    return H;
}

/*
 * Symbolic Transfer Entropy:
 *   TE_{Y->X} = sum p(x_{t+1}, x_t^{(m)}, y_t^{(m)}) *
 *               log( p(x_{t+1} | x_t^{(m)}, y_t^{(m)}) /
 *                    p(x_{t+1} | x_t^{(m)}) )
 *
 * where probabilities are estimated from ordinal pattern counts.
 */
double te_symbolic_transfer_entropy(const TETimeSeries *x, const TETimeSeries *y,
                                     int order, int delay) {
    if (!x || !y || order < 2 || delay < 1) return 0.0;
    if (x->length != y->length) return 0.0;

    int n_patterns = 1;
    for (int i = 2; i <= order; i++) n_patterns *= i;
    int n = x->length - order * delay;
    if (n < 50) return 0.0;

    /* Build joint probability table:
     * Joint = p(x_future_pattern, x_past_pattern, y_past_pattern)
     */
    size_t total = (size_t)n_patterns * (size_t)n_patterns * (size_t)n_patterns;
    double *joint = calloc(total, sizeof(double));
    double *marg_xx = calloc((size_t)n_patterns * (size_t)n_patterns, sizeof(double));
    double *marg_x = calloc((size_t)n_patterns, sizeof(double));
    double *marg_xy = calloc((size_t)n_patterns * (size_t)n_patterns, sizeof(double));
    if (!joint || !marg_xx || !marg_x || !marg_xy) {
        free(joint); free(marg_xx); free(marg_x); free(marg_xy); return 0.0;
    }

    double window[8]; int pat[8];
    double inv_n = 1.0 / (double)n;

    for (int t = 0; t < n; t++) {
        /* X future: 1 step ahead */
        for (int d = 0; d < order; d++) window[d] = x->data[t + 1 + d * delay];
        int xf = te_symbolic_permutation_encode(window, order, pat);

        /* X past: current state */
        for (int d = 0; d < order; d++) window[d] = x->data[t + d * delay];
        int xp = te_symbolic_permutation_encode(window, order, pat);

        /* Y past: current state */
        for (int d = 0; d < order; d++) window[d] = y->data[t + d * delay];
        int yp = te_symbolic_permutation_encode(window, order, pat);

        if (xf < 0 || xp < 0 || yp < 0) continue;

        size_t idx = ((size_t)xf * (size_t)n_patterns + (size_t)xp) * (size_t)n_patterns + (size_t)yp;
        joint[idx] += inv_n;
        marg_xx[(size_t)xf * (size_t)n_patterns + (size_t)xp] += inv_n;
        marg_x[xp] += inv_n;
        marg_xy[(size_t)xp * (size_t)n_patterns + (size_t)yp] += inv_n;
    }

    double TE = 0.0;
    for (int xf = 0; xf < n_patterns; xf++) {
        for (int xp = 0; xp < n_patterns; xp++) {
            for (int yp = 0; yp < n_patterns; yp++) {
                size_t idx = ((size_t)xf * (size_t)n_patterns + (size_t)xp) * (size_t)n_patterns + (size_t)yp;
                double p_xyz = joint[idx];
                if (p_xyz < TE_EPSILON) continue;
                double p_xx = marg_xx[(size_t)xf * (size_t)n_patterns + (size_t)xp];
                double p_xp = marg_x[xp];
                double p_xpyp = marg_xy[(size_t)xp * (size_t)n_patterns + (size_t)yp];
                if (p_xx < TE_EPSILON || p_xp < TE_EPSILON || p_xpyp < TE_EPSILON) continue;
                TE += p_xyz * log(p_xyz * p_xp / (p_xx * p_xpyp));
            }
        }
    }

    free(joint); free(marg_xx); free(marg_x); free(marg_xy);
    return (TE > 0.0) ? TE : 0.0;
}

TEResult te_symbolic_compute(const TETimeSeries *x, const TETimeSeries *y,
                              int order, int delay) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    result.k_embed = order; result.l_embed = order;
    if (!x || !y || order < 2 || delay < 1) return result;
    result.te_xy = te_symbolic_transfer_entropy(x, y, order, delay);
    result.te_yx = te_symbolic_transfer_entropy(y, x, order, delay);
    result.net_te = result.te_xy - result.te_yx;
    result.h_x_future = te_symbolic_permutation_entropy(x, order, delay);
    result.n_surrogates = 0;
    return result;
}

/* ???????????????????????????????????????????????????????????
 * Weighted Permutation Entropy (Fadlallah et al., 2013)
 * Each pattern is weighted by the variance of the embedding window.
 * ??????????????????????????????????????????????????????????? */
double te_symbolic_weighted_permutation_entropy(const TETimeSeries *ts, int order, int delay) {
    if (!ts || order < 2 || order > 8 || delay < 1) return 0.0;
    int n_patterns = 1; for (int i = 2; i <= order; i++) n_patterns *= i;
    int n = ts->length - (order - 1) * delay; if (n < 1) return 0.0;
    double *weights = calloc((size_t)n_patterns, sizeof(double));
    int *counts = calloc((size_t)n_patterns, sizeof(int));
    if (!weights || !counts) { free(weights); free(counts); return 0.0; }
    double total_weight = 0.0, window[8]; int pattern[8];
    for (int t = 0; t < n; t++) {
        double mean = 0.0, var = 0.0;
        for (int d = 0; d < order; d++) { window[d] = ts->data[t + d * delay]; mean += window[d]; }
        mean /= (double)order;
        for (int d = 0; d < order; d++) { double dev = window[d] - mean; var += dev * dev; }
        double w = var / (double)order;
        int idx = te_symbolic_permutation_encode(window, order, pattern);
        if (idx >= 0 && idx < n_patterns) { weights[idx] += w; counts[idx]++; total_weight += w; }
    }
    double H = 0.0;
    for (int i = 0; i < n_patterns; i++) {
        if (weights[i] > TE_EPSILON && total_weight > TE_EPSILON) {
            double p = weights[i] / total_weight; H -= p * log(p);
        }
    }
    free(weights); free(counts); return H;
}

/* ???????????????????????????????????????????????????????????
 * Multi-scale Permutation Entropy
 * Coarse-graining followed by permutation entropy at each scale.
 * ??????????????????????????????????????????????????????????? */
void te_symbolic_multiscale_pe(const TETimeSeries *ts, int order, int delay,
                                int max_scale, double *pe_scales) {
    if (!ts || !pe_scales || max_scale < 1) return;
    TETimeSeries *tss = te_ts_create(ts->data, ts->length, "temp");
    if (!tss) return;
    for (int scale = 1; scale <= max_scale; scale++) {
        int n_new = ts->length / scale; if (n_new < order) { pe_scales[scale-1] = 0.0; continue; }
        double *cg = malloc((size_t)n_new * sizeof(double));
        if (!cg) continue;
        memset(cg, 0, (size_t)n_new * sizeof(double));
        for (int i = 0; i < n_new * scale; i++) cg[i / scale] += ts->data[i];
        for (int i = 0; i < n_new; i++) cg[i] /= (double)scale;
        free(tss->data); tss->data = cg; tss->length = n_new;
        pe_scales[scale-1] = te_symbolic_permutation_entropy(tss, order, delay);
    }
    te_ts_free(tss);
}

/* ???????????????????????????????????????????????????????????
 * Ordinal pattern transition matrix
 * Computes transition probabilities between ordinal patterns.
 * ??????????????????????????????????????????????????????????? */
int te_symbolic_transition_matrix(const TETimeSeries *ts, int order, int delay,
                                   double *trans_matrix) {
    if (!ts || !trans_matrix || order < 2) return -1;
    int np = 1; for (int i = 2; i <= order; i++) np *= i;
    int n = ts->length - order * delay; if (n < 2) return -1;
    memset(trans_matrix, 0, (size_t)(np * np) * sizeof(double));
    int *counts = calloc((size_t)np, sizeof(int));
    double window[8]; int pattern[8];
    int prev = -1;
    for (int t = 0; t < n; t++) {
        for (int d = 0; d < order; d++) window[d] = ts->data[t + d * delay];
        int curr = te_symbolic_permutation_encode(window, order, pattern);
        if (curr < 0) continue;
        if (prev >= 0) { trans_matrix[prev * np + curr] += 1.0; counts[prev]++; }
        prev = curr;
    }
    for (int i = 0; i < np; i++) if (counts[i] > 0) for (int j = 0; j < np; j++) trans_matrix[i * np + j] /= (double)counts[i];
    free(counts); return 0;
}

/* ?? Amplitude-aware permutation entropy ?????????????????? */
double te_symbolic_aape(const TETimeSeries *ts, int order, int delay, double A) {
    if (!ts || order<2 || order>8 || delay<1 || A<0.0) return 0.0;
    int np=1; for(int i=2;i<=order;i++) np*=i;
    int n=ts->length-(order-1)*delay; if(n<1) return 0.0;
    double *weights=calloc((size_t)np,sizeof(double)); if(!weights)return 0.0;
    double window[8]; int pattern[8]; double total_w=0.0;
    for(int t=0;t<n;t++) {
        for(int d=0;d<order;d++) window[d]=ts->data[t+d*delay];
        int idx=te_symbolic_permutation_encode(window,order,pattern);
        if(idx<0||idx>=np) continue;
        double mean=0; for(int d=0;d<order;d++) mean+=fabs(window[d]); mean/=(double)order;
        double w=A/(A+mean+1e-15);
        weights[idx]+=w; total_w+=w;
    }
    double H=0.0;
    for(int i=0;i<np;i++) if(weights[i]>TE_EPSILON&&total_w>TE_EPSILON) { double p=weights[i]/total_w; H-=p*log(p); }
    free(weights); return H;
}

/* ?? Conditional permutation entropy ?????????????????????? */
double te_symbolic_conditional_pe(const TETimeSeries *ts, int order, int delay) {
    if (!ts||order<2||order>6||delay<1) return 0.0;
    int np=1; for(int i=2;i<=order;i++) np*=i;
    int n=ts->length-order*delay; if(n<50) return 0.0;
    size_t sz=(size_t)np*(size_t)np;
    double *trans=calloc(sz,sizeof(double)), *marg=calloc((size_t)np,sizeof(double));
    if(!trans||!marg) { free(trans);free(marg); return 0.0; }
    double window[8]; int pattern[8]; int prev=-1;
    for(int t=0;t<n;t++) {
        for(int d=0;d<order;d++) window[d]=ts->data[t+d*delay];
        int curr=te_symbolic_permutation_encode(window,order,pattern);
        if(curr<0||curr>=np) continue;
        if(prev>=0) { trans[prev*np+curr]+=1.0; marg[prev]+=1.0; }
        prev=curr;
    }
    double H=0.0;
    for(int i=0;i<np;i++) if(marg[i]>TE_EPSILON) for(int j=0;j<np;j++) if(trans[i*np+j]>TE_EPSILON)
        H -= (trans[i*np+j]/(double)n) * log(trans[i*np+j]/marg[i]+1e-15);
    free(trans);free(marg); return H;
}

/* ?? Symbolic TE with bias correction ????????????????????? */
double te_symbolic_bias_corrected(const TETimeSeries *x, const TETimeSeries *y, int order, int delay, int n_shuffles) {
    if (!x||!y) return 0.0;
    double raw=te_symbolic_transfer_entropy(x,y,order,delay);
    if (n_shuffles<1) return raw;
    double *shuf=malloc((size_t)y->length*sizeof(double));
    if(!shuf) return raw;
    double bias=0.0;
    for(int s=0;s<n_shuffles;s++) {
        te_generate_iid_surrogates(y->data,y->length,shuf);
        TETimeSeries ys={.data=shuf,.length=y->length,.name="shuf"};
        bias+=te_symbolic_transfer_entropy(x,&ys,order,delay);
    }
    free(shuf);
    bias/=n_shuffles;
    double corr=raw-bias;
    return corr>0.0?corr:0.0;
}
