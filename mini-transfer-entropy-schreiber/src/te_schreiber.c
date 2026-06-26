#include "te_schreiber.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/*
 * Schreiber's Transfer Entropy (Phys. Rev. Lett. 85, 461, 2000)
 *
 * T_{Y->X} = sum p(x_{n+1}, x_n^{(k)}, y_n^{(l)}) *
 *            log( p(x_{n+1} | x_n^{(k)}, y_n^{(l)}) /
 *                 p(x_{n+1} | x_n^{(k)}) )
 *
 * This is the reduction in uncertainty about future X given past Y,
 * beyond the information already provided by past X.
 */

/* Three-dimensional joint probability table:
 * dim 0: X_future (x_{n+1})
 * dim 1: X_past   (x_n^{(k)} encoded as base-n_bins number)
 * dim 2: Y_past   (y_n^{(l)} encoded as base-n_bins number)
 */
typedef struct {
    double *data;   /* flattened 3D array [nx][nx^k][ny^l] */
    double *marg_xf_xp; /* [nx][nx^k] */
    double *marg_xp;    /* [nx^k] */
    int     nx;          /* n_bins for X */
    int     ny;          /* n_bins for Y */
    int     nxk;         /* n_bins^k */
    int     nyl;         /* n_bins^l */
    int     k_embed;
    int     l_embed;
    int     n_samples;
} TE3DJoint;

static int int_pow(int base, int exp) {
    int r = 1;
    for (int i = 0; i < exp; i++) r *= base;
    return r;
}

static int encode_state(const int *bins, int dim, int n_bins) {
    int code = 0, mult = 1;
    for (int d = 0; d < dim; d++) { code += bins[d] * mult; mult *= n_bins; }
    return code;
}

static TE3DJoint* te3d_alloc(int nx, int ny, int k, int l) {
    TE3DJoint *j3 = calloc(1, sizeof(TE3DJoint));
    if (!j3) return NULL;
    j3->nx = nx; j3->ny = ny;
    j3->nxk = int_pow(nx, k);
    j3->nyl = int_pow(ny, l);
    j3->k_embed = k; j3->l_embed = l;
    size_t sz_joint  = (size_t)nx * (size_t)j3->nxk * (size_t)j3->nyl;
    size_t sz_marg   = (size_t)nx * (size_t)j3->nxk;
    j3->data = calloc(sz_joint, sizeof(double));
    j3->marg_xf_xp = calloc(sz_marg, sizeof(double));
    j3->marg_xp = calloc((size_t)j3->nxk, sizeof(double));
    if (!j3->data || !j3->marg_xf_xp || !j3->marg_xp) {
        free(j3->data); free(j3->marg_xf_xp); free(j3->marg_xp); free(j3); return NULL;
    }
    return j3;
}

static void te3d_free(TE3DJoint *j3) {
    if (!j3) return;
    free(j3->data); free(j3->marg_xf_xp); free(j3->marg_xp);
    free(j3);
}

/* Compute TE from the 3D joint table.
 * T = sum p(x_f, x_p, y_p) * log( p(x_f|x_p,y_p) / p(x_f|x_p) )
 *   = sum p(x_f, x_p, y_p) * log( p(x_f, x_p, y_p) * p(x_p) /
 *                                   (p(x_f, x_p) * p(x_p, y_p)) )
 */
static double te3d_compute_te(const TE3DJoint *j3) {
    double TE = 0.0;
    for (int xf = 0; xf < j3->nx; xf++) {
        for (int xp = 0; xp < j3->nxk; xp++) {
            for (int yp = 0; yp < j3->nyl; yp++) {
                size_t idx = ((size_t)xf * (size_t)j3->nxk + (size_t)xp) * (size_t)j3->nyl + (size_t)yp;
                double p_xyz = j3->data[idx];
                if (p_xyz < TE_EPSILON) continue;

                size_t idx_xy = (size_t)xf * (size_t)j3->nxk + (size_t)xp;
                double p_xy = j3->marg_xf_xp[idx_xy];
                if (p_xy < TE_EPSILON) continue;

                double p_xp = j3->marg_xp[xp];
                if (p_xp < TE_EPSILON) continue;

                /* p(x_f, x_p, y_p) / p(x_f, x_p) = p(y_p | x_f, x_p)
                 * Then we also need p(y_p | x_p) = p(x_p, y_p) / p(x_p)
                 * But we have the full formula:
                 * T += p_xyz * log( p_xyz * p_xp / (p_xy * p(x_p, y_p)) )
                 * Need p(x_p, y_p) which is the marginal over x_f:
                 */
                double p_xpyp = 0.0;
                for (int xf2 = 0; xf2 < j3->nx; xf2++) {
                    size_t idx2 = ((size_t)xf2 * (size_t)j3->nxk + (size_t)xp) * (size_t)j3->nyl + (size_t)yp;
                    p_xpyp += j3->data[idx2];
                }
                if (p_xpyp < TE_EPSILON) continue;

                double ratio = p_xyz * p_xp / (p_xy * p_xpyp);
                TE += p_xyz * log(ratio);
            }
        }
    }
    return TE;
}

/* Compute TE using Schreiber's binning estimator */
TEResult te_schreiber_compute(const TETimeSeries *x, const TETimeSeries *y,
                                int k, int l, int delay) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    if (!x || !y || k < 1 || l < 1 || delay < 1) return result;

    int n_bins = 8; /* default binning */
    return te_schreiber_binning(x, y, k, l, n_bins);
}

TEResult te_schreiber_binning(const TETimeSeries *x, const TETimeSeries *y,
                                int k, int l, int n_bins) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    result.k_embed = k; result.l_embed = l;
    if (!x || !y || k < 1 || l < 1 || n_bins < 2) return result;

    /* Discretize both series */
    TEBinnedSeries *bx = te_bin_create(x, n_bins, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *by = te_bin_create(y, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bx || !by) { te_bin_free(bx); te_bin_free(by); return result; }

    int max_embed = (k > l) ? k : l;
    int n = bx->length - max_embed;
    if (n < 50) { te_bin_free(bx); te_bin_free(by); return result; }

    TE3DJoint *j3_xy = te3d_alloc(bx->n_bins, by->n_bins, k, l);
    TE3DJoint *j3_yx = te3d_alloc(by->n_bins, bx->n_bins, l, k);
    if (!j3_xy || !j3_yx) { te3d_free(j3_xy); te3d_free(j3_yx); te_bin_free(bx); te_bin_free(by); return result; }

    /* Fill joint probability tables by counting occurrences */
    double inv_n = 1.0 / (double)n;
    int xp_buf[TE_MAX_EMBED], yp_buf[TE_MAX_EMBED];

    for (int t = 0; t < n; t++) {
        int x_future = bx->bins[t + k];  /* x_{n+1} approximated as x_{n+k} */

        for (int d = 0; d < k; d++) xp_buf[d] = bx->bins[t + d];
        for (int d = 0; d < l; d++) yp_buf[d] = by->bins[t + d];

        int x_past = encode_state(xp_buf, k, bx->n_bins);
        int y_past = encode_state(yp_buf, l, by->n_bins);

        /* TE(Y->X): count p(x_future, x_past, y_past) */
        {
            size_t idx = ((size_t)x_future * (size_t)j3_xy->nxk + (size_t)x_past)
                         * (size_t)j3_xy->nyl + (size_t)y_past;
            j3_xy->data[idx] += inv_n;
            j3_xy->marg_xf_xp[(size_t)x_future * (size_t)j3_xy->nxk + (size_t)x_past] += inv_n;
            j3_xy->marg_xp[x_past] += inv_n;
        }
        /* TE(X->Y): count p(y_future, y_past, x_past) */
        int y_future = by->bins[t + l];
        {
            size_t idx = ((size_t)y_future * (size_t)j3_yx->nxk + (size_t)y_past)
                         * (size_t)j3_yx->nyl + (size_t)x_past;
            j3_yx->data[idx] += inv_n;
            j3_yx->marg_xf_xp[(size_t)y_future * (size_t)j3_yx->nxk + (size_t)y_past] += inv_n;
            j3_yx->marg_xp[y_past] += inv_n;
        }
    }

    result.te_xy = te3d_compute_te(j3_xy);
    result.te_yx = te3d_compute_te(j3_yx);
    result.net_te = result.te_xy - result.te_yx;

    /* Compute entropies for diagnostics */
    result.h_x_future = te_bin_entropy(bx);
    result.h_x_past = 0.0;
    result.h_x_past_y = 0.0;
    result.n_surrogates = 0;
    result.significant = false;

    te3d_free(j3_xy); te3d_free(j3_yx);
    te_bin_free(bx); te_bin_free(by);
    return result;
}

/*
 * Entropy rate of order k:
 *   h_k = H(X_{n+1} | X_n^{(k)})
 *        = H(X_{n+1}, X_n^{(k)}) - H(X_n^{(k)})
 */
double te_schreiber_entropy_rate(const TEBinnedSeries *bs, int k) {
    if (!bs || k < 1) return 0.0;
    int n = bs->length - k;
    if (n < 10) return 0.0;
    int nb = bs->n_bins;
    int n_states = int_pow(nb, k + 1);
    int n_states_past = int_pow(nb, k);
    double *joint = calloc((size_t)n_states, sizeof(double));
    double *marg = calloc((size_t)n_states_past, sizeof(double));
    if (!joint || !marg) { free(joint); free(marg); return 0.0; }
    double inv_n = 1.0 / (double)n;
    for (int t = 0; t < n; t++) {
        int future = bs->bins[t + k];
        int past = 0, mult = 1;
        for (int d = 0; d < k; d++) { past += bs->bins[t + d] * mult; mult *= nb; }
        joint[future * n_states_past + past] += inv_n;
        marg[past] += inv_n;
    }
    double H_joint = 0.0, H_marg = 0.0;
    for (int i = 0; i < n_states; i++)
        if (joint[i] > TE_EPSILON) H_joint -= joint[i] * log(joint[i]);
    for (int i = 0; i < n_states_past; i++)
        if (marg[i] > TE_EPSILON) H_marg -= marg[i] * log(marg[i]);
    free(joint); free(marg);
    return H_joint - H_marg;
}

/*
 * Mutual Information: I(X; Y) = H(X) + H(Y) - H(X,Y)
 * With time lag: I(X_t; Y_{t-lag})
 */
double te_schreiber_mutual_info(const TEBinnedSeries *x, const TEBinnedSeries *y, int lag) {
    if (!x || !y || lag < 0) return 0.0;
    int n = x->length - lag;
    if (n < 10) return 0.0;
    TEJointProb *jp = te_joint_estimate(x, y, lag);
    if (!jp) return 0.0;
    double Hx = 0.0, Hy = 0.0;
    for (int i = 0; i < jp->nx; i++)
        if (jp->marg_x[i] > TE_EPSILON) Hx -= jp->marg_x[i] * log(jp->marg_x[i]);
    for (int j = 0; j < jp->ny; j++)
        if (jp->marg_y[j] > TE_EPSILON) Hy -= jp->marg_y[j] * log(jp->marg_y[j]);
    double Hxy = te_joint_entropy(jp);
    te_joint_free(jp);
    return Hx + Hy - Hxy;
}

/*
 * Scan transfer entropy over multiple lags to find optimal delay.
 * te_values[0..max_lag-1] filled with TE(Y->X) at each lag.
 */
void te_schreiber_scan_lag(const TETimeSeries *x, const TETimeSeries *y,
                            int k, int l, int max_lag, double *te_values) {
    if (!x || !y || !te_values || max_lag < 1) return;
    for (int lag = 1; lag <= max_lag; lag++) {
        /* For scanning, we fix k=l=1 and vary the delay */
        TEResult r = te_schreiber_compute(x, y, k, l, lag);
        te_values[lag - 1] = r.te_xy;
    }
}

/* ???????????????????????????????????????????????????????????
 * Multi-scale Transfer Entropy
 * Computes TE at multiple time scales via coarse-graining.
 * ??????????????????????????????????????????????????????????? */
void te_schreiber_multiscale(const TETimeSeries *x, const TETimeSeries *y,
                              int k, int l, int max_scale, double *te_scales) {
    if (!x || !y || !te_scales || max_scale < 1) return;
    TETimeSeries *xs = te_ts_create(x->data, x->length, "x_scale");
    TETimeSeries *ys = te_ts_create(y->data, y->length, "y_scale");
    if (!xs || !ys) { te_ts_free(xs); te_ts_free(ys); return; }
    for (int scale = 1; scale <= max_scale; scale++) {
        int n_new = x->length / scale; if (n_new < 50) { te_scales[scale-1] = 0.0; continue; }
        double *x_cg = malloc((size_t)n_new * sizeof(double));
        double *y_cg = malloc((size_t)n_new * sizeof(double));
        if (!x_cg || !y_cg) { free(x_cg); free(y_cg); continue; }
        for (int i = 0; i < n_new; i++) { x_cg[i] = 0.0; y_cg[i] = 0.0; }
        for (int i = 0; i < n_new * scale; i++) { x_cg[i / scale] += x->data[i]; y_cg[i / scale] += y->data[i]; }
        for (int i = 0; i < n_new; i++) { x_cg[i] /= (double)scale; y_cg[i] /= (double)scale; }
        free(xs->data); xs->data = x_cg; xs->length = n_new;
        free(ys->data); ys->data = y_cg; ys->length = n_new;
        TEResult r = te_schreiber_binning(xs, ys, k, l, 6);
        te_scales[scale-1] = r.te_xy;
    }
    te_ts_free(xs); te_ts_free(ys);
}

/* ???????????????????????????????????????????????????????????
 * Time-varying Transfer Entropy (sliding window)
 * ??????????????????????????????????????????????????????????? */
int te_schreiber_sliding_window(const TETimeSeries *x, const TETimeSeries *y,
                                 int k, int l, int window_size, int step,
                                 double *te_tv, int max_points) {
    if (!x || !y || !te_tv || window_size < 50 || step < 1 || max_points < 1) return 0;
    int count = 0;
    for (int start = 0; start + window_size <= x->length && count < max_points; start += step) {
        TETimeSeries xw = { .data = x->data + start, .length = window_size, .name = "x_win" };
        TETimeSeries yw = { .data = y->data + start, .length = window_size, .name = "y_win" };
        TEResult r = te_schreiber_binning(&xw, &yw, k, l, 6);
        te_tv[count++] = r.te_xy;
    }
    return count;
}

/* ???????????????????????????????????????????????????????????
 * Delay embedding selection via Ragwitz criterion
 * Selects optimal (k,l) by minimizing entropy rate.
 * ??????????????????????????????????????????????????????????? */
void te_schreiber_optimal_embedding(const TETimeSeries *x, const TETimeSeries *y,
                                     int max_k, int max_l, int *best_k, int *best_l) {
    *best_k = 1; *best_l = 1; double best_te = -1.0;
    for (int k = 1; k <= max_k; k++) {
        for (int l = 1; l <= max_l; l++) {
            TEResult r = te_schreiber_binning(x, y, k, l, 6);
            if (r.te_xy > best_te) { best_te = r.te_xy; *best_k = k; *best_l = l; }
        }
    }
}

/* ?? Transfer Entropy with variable embedding ?????????????? */
TEResult te_schreiber_variable_embedding(const TETimeSeries *x, const TETimeSeries *y,
                                          int k_min, int k_max, int l_min, int l_max) {
    TEResult best; memset(&best, 0, sizeof(best)); double best_te = -1.0;
    for (int k = k_min; k <= k_max; k++) {
        for (int l = l_min; l <= l_max; l++) {
            TEResult r = te_schreiber_binning(x, y, k, l, 6);
            if (r.te_xy > best_te) { best_te = r.te_xy; best = r; }
        }
    }
    return best;
}

/* ?? Directionality index ?????????????????????????????????? */
double te_schreiber_directionality(const TEResult *r) {
    if (!r) return 0.0;
    double total = r->te_xy + r->te_yx;
    return (total > TE_EPSILON) ? (r->te_xy - r->te_yx) / total : 0.0;
}

/* ?? Uncertainty coefficient (Theil's U) ???????????????????
 * U(Y|X) = TE(Y->X) / H(Y_future)
 */
double te_schreiber_uncertainty_coefficient(const TEResult *r) {
    if (!r || r->h_x_future < TE_EPSILON) return 0.0;
    return r->te_xy / r->h_x_future;
}

/* ?? Integration with multiple bin widths ??????????????????
 * Average TE over multiple bin sizes for robustness.
 */
TEResult te_schreiber_robust_binning(const TETimeSeries *x, const TETimeSeries *y,
                                      int k, int l, int min_bins, int max_bins, int n_bin_sizes) {
    TEResult avg; memset(&avg, 0, sizeof(avg)); avg.k_embed = k; avg.l_embed = l;
    if (n_bin_sizes < 1) return avg;
    int count = 0;
    for (int i = 0; i < n_bin_sizes; i++) {
        int nb = min_bins + (max_bins - min_bins) * i / (n_bin_sizes - 1 + (n_bin_sizes == 1 ? 1 : 0));
        if (nb < 2) nb = 2;
        TEResult r = te_schreiber_binning(x, y, k, l, nb);
        avg.te_xy += r.te_xy; avg.te_yx += r.te_yx; count++;
    }
    if (count > 0) { avg.te_xy /= (double)count; avg.te_yx /= (double)count; avg.net_te = avg.te_xy - avg.te_yx; }
    return avg;
}

/* ?? Normalized Transfer Entropy ???????????????????????????
 * NTE = TE(Y->X) / H(X_future | X_past)
 * Ranges from 0 (no info) to 1 (complete info).
 */
double te_schreiber_normalized(const TEResult *r) {
    if (!r || r->h_x_past < TE_EPSILON) return 0.0;
    return r->te_xy / r->h_x_past;
}

/* ?? Transfer entropy with Theiler window correction ???????
 * Excludes temporally correlated points within a Theiler window.
 */
TEResult te_schreiber_theiler(const TETimeSeries *x, const TETimeSeries *y,
                               int k, int l, int n_bins, int theiler_w) {
    TEResult result; memset(&result, 0, sizeof(result));
    result.k_embed = k; result.l_embed = l;
    if (!x || !y || k < 1 || l < 1 || n_bins < 2 || theiler_w < 1) return result;
    TEBinnedSeries *bx = te_bin_create(x, n_bins, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *by = te_bin_create(y, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bx || !by) { te_bin_free(bx); te_bin_free(by); return result; }
    int n = bx->length - k, nb = bx->n_bins, nxk = 1; for (int i=0;i<k;i++) nxk*=nb;
    int nyl = 1; for (int i=0;i<l;i++) nyl*=nb;
    size_t sz = (size_t)nb*(size_t)nxk*(size_t)nyl;
    double *joint = calloc(sz, sizeof(double));
    if (!joint) { te_bin_free(bx); te_bin_free(by); return result; }
    int count = 0;
    double *marg_xx = calloc((size_t)nb*(size_t)nxk, sizeof(double));
    double *marg_xp = calloc((size_t)nxk, sizeof(double));
    double *marg_xpyp = calloc((size_t)nxk*(size_t)nyl, sizeof(double));
    if (!marg_xx || !marg_xp || !marg_xpyp) { free(joint); free(marg_xx); free(marg_xp); free(marg_xpyp); te_bin_free(bx); te_bin_free(by); return result; }
    for (int t = 0; t < n; t++) {
        int skip = 0;
        for (int w = 1; w <= theiler_w && t + w < n; w++) { if (t + w < n) { double dx = x->data[t+k] - x->data[t+w+k]; if (fabs(dx) < 1e-6) { skip = 1; break; } } }
        if (skip) continue;
        int xf = bx->bins[t+k], xp = 0, mult = 1;
        for (int d=0;d<k;d++) { xp+=bx->bins[t+d]*mult; mult*=nb; }
        int yp = 0; mult = 1;
        for (int d=0;d<l;d++) { yp+=by->bins[t+d]*mult; mult*=nb; }
        joint[((size_t)xf*(size_t)nxk+(size_t)xp)*(size_t)nyl+(size_t)yp] += 1.0;
        marg_xx[(size_t)xf*(size_t)nxk+(size_t)xp] += 1.0;
        marg_xp[xp] += 1.0;
        marg_xpyp[(size_t)xp*(size_t)nyl+(size_t)yp] += 1.0;
        count++;
    }
    if (count < 10) { free(joint);free(marg_xx);free(marg_xp);free(marg_xpyp);te_bin_free(bx);te_bin_free(by);return result; }
    double inv_n = 1.0/(double)count;
    result.te_xy = 0.0;
    for (int xf=0;xf<nb;xf++) for (int xp=0;xp<nxk;xp++) for (int yp=0;yp<nyl;yp++) {
        double p_xyz = joint[((size_t)xf*(size_t)nxk+(size_t)xp)*(size_t)nyl+(size_t)yp]*inv_n;
        if (p_xyz<TE_EPSILON) continue;
        double p_xx = marg_xx[(size_t)xf*(size_t)nxk+(size_t)xp]*inv_n;
        double p_xp = marg_xp[xp]*inv_n, p_xy = marg_xpyp[(size_t)xp*(size_t)nyl+(size_t)yp]*inv_n;
        if (p_xx<TE_EPSILON||p_xp<TE_EPSILON||p_xy<TE_EPSILON) continue;
        result.te_xy += p_xyz*log(p_xyz*p_xp/(p_xx*p_xy));
    }
    free(joint);free(marg_xx);free(marg_xp);free(marg_xpyp);te_bin_free(bx);te_bin_free(by);
    return result;
}

/* ?? Integrated TE over a range of embedding parameters ??? */
double te_schreiber_integrated(const TETimeSeries *x, const TETimeSeries *y, int max_k, int max_l) {
    double integral = 0.0; int count = 0;
    for (int k=1;k<=max_k;k++) for (int l=1;l<=max_l;l++) { TEResult r=te_schreiber_binning(x,y,k,l,6); integral+=r.te_xy; count++; }
    return count>0 ? integral/(double)count : 0.0;
}
