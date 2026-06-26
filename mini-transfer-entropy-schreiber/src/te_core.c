#include "te_core.h"
#include "te_schreiber.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

TETimeSeries* te_ts_create(const double *data, int len, const char *name) {
    if (!data || len < 1) return NULL;
    TETimeSeries *ts = calloc(1, sizeof(TETimeSeries));
    if (!ts) return NULL;
    ts->data = malloc((size_t)len * sizeof(double));
    if (!ts->data) { free(ts); return NULL; }
    memcpy(ts->data, data, (size_t)len * sizeof(double));
    ts->length = len;
    ts->name = name ? strdup(name) : strdup("unnamed");
    return ts;
}
void te_ts_free(TETimeSeries *ts) { if (!ts) return; free(ts->data); free(ts->name); free(ts); }

void te_ts_normalize(TETimeSeries *ts) {
    if (!ts || ts->length < 1) return;
    double mn = te_ts_min(ts), mx = te_ts_max(ts), rng = mx - mn;
    if (rng < TE_EPSILON) { for (int i = 0; i < ts->length; i++) ts->data[i] = 0.5; return; }
    for (int i = 0; i < ts->length; i++) ts->data[i] = (ts->data[i] - mn) / rng;
}
void te_ts_standardize(TETimeSeries *ts) {
    if (!ts || ts->length < 2) return;
    double m = te_ts_mean(ts), s = te_ts_std(ts);
    if (s < TE_EPSILON) { for (int i = 0; i < ts->length; i++) ts->data[i] = 0.0; return; }
    for (int i = 0; i < ts->length; i++) ts->data[i] = (ts->data[i] - m) / s;
}
double te_ts_mean(const TETimeSeries *ts) {
    if (!ts || ts->length < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < ts->length; i++) sum += ts->data[i];
    return sum / (double)ts->length;
}
double te_ts_variance(const TETimeSeries *ts) {
    if (!ts || ts->length < 2) return 0.0;
    double m = te_ts_mean(ts), sum = 0.0;
    for (int i = 0; i < ts->length; i++) { double d = ts->data[i] - m; sum += d * d; }
    return sum / (double)(ts->length - 1);
}
double te_ts_std(const TETimeSeries *ts) { return sqrt(te_ts_variance(ts)); }
double te_ts_min(const TETimeSeries *ts) {
    if (!ts || ts->length < 1) return 0.0;
    double mn = ts->data[0];
    for (int i = 1; i < ts->length; i++) if (ts->data[i] < mn) mn = ts->data[i];
    return mn;
}
double te_ts_max(const TETimeSeries *ts) {
    if (!ts || ts->length < 1) return 0.0;
    double mx = ts->data[0];
    for (int i = 1; i < ts->length; i++) if (ts->data[i] > mx) mx = ts->data[i];
    return mx;
}
int te_ts_subsample(const TETimeSeries *ts, int step, TETimeSeries **out) {
    if (!ts || !out || step < 1) return -1;
    int n_out = ts->length / step;
    if (n_out < 2) return -1;
    double *buf = malloc((size_t)n_out * sizeof(double));
    if (!buf) return -1;
    for (int i = 0; i < n_out; i++) buf[i] = ts->data[i * step];
    *out = te_ts_create(buf, n_out, "subsampled");
    free(buf);
    return (*out) ? 0 : -1;
}
void te_ts_diff(TETimeSeries *ts) {
    if (!ts || ts->length < 2) return;
    for (int i = 0; i < ts->length - 1; i++) ts->data[i] = ts->data[i + 1] - ts->data[i];
    ts->length--;
}

TEBinnedSeries* te_bin_create(const TETimeSeries *ts, int n_bins, TEBinningMethod method) {
    if (!ts || n_bins < 2 || n_bins > TE_MAX_BINS) return NULL;
    TEBinnedSeries *bs = calloc(1, sizeof(TEBinnedSeries));
    if (!bs) return NULL;
    bs->bins = malloc((size_t)ts->length * sizeof(int));
    bs->edges = malloc((size_t)(n_bins + 1) * sizeof(double));
    if (!bs->bins || !bs->edges) { te_bin_free(bs); return NULL; }
    bs->n_bins = n_bins; bs->length = ts->length; bs->method = method;
    double mn = te_ts_min(ts), mx = te_ts_max(ts), rng = mx - mn;
    if (rng < TE_EPSILON) rng = 1.0;
    if (method == TE_BIN_EQUAL_WIDTH || method == TE_BIN_KD_TREE || method == TE_BIN_SYMBOLIC) {
        for (int i = 0; i <= n_bins; i++) bs->edges[i] = mn + (double)i * rng / (double)n_bins;
        for (int i = 0; i < ts->length; i++) {
            int b = (int)((double)n_bins * (ts->data[i] - mn) / rng);
            if (b < 0) b = 0; if (b >= n_bins) b = n_bins - 1;
            bs->bins[i] = b;
        }
    } else if (method == TE_BIN_EQUAL_FREQ) {
        double *sorted = malloc((size_t)ts->length * sizeof(double));
        if (!sorted) { te_bin_free(bs); return NULL; }
        memcpy(sorted, ts->data, (size_t)ts->length * sizeof(double));
        for (int i = 0; i < ts->length - 1; i++)
            for (int j = i + 1; j < ts->length; j++)
                if (sorted[i] > sorted[j]) { double t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }
        for (int b = 0; b <= n_bins; b++) { int idx = (b * ts->length) / n_bins; if (idx >= ts->length) idx = ts->length - 1; bs->edges[b] = sorted[idx]; }
        for (int i = 0; i < ts->length; i++) { int b = 0; while (b < n_bins - 1 && ts->data[i] > bs->edges[b + 1]) b++; bs->bins[i] = b; }
        free(sorted);
    }
    return bs;
}
void te_bin_free(TEBinnedSeries *bs) { if (bs) { free(bs->bins); free(bs->edges); free(bs); } }
int te_bin_get(const TEBinnedSeries *bs, int idx) { return (bs && idx >= 0 && idx < bs->length) ? bs->bins[idx] : -1; }

void te_bin_count(const TEBinnedSeries *bs, double *hist) {
    if (!bs || !hist) return;
    memset(hist, 0, (size_t)bs->n_bins * sizeof(double));
    for (int i = 0; i < bs->length; i++) hist[bs->bins[i]] += 1.0;
    double inv_n = 1.0 / (double)bs->length;
    for (int b = 0; b < bs->n_bins; b++) hist[b] *= inv_n;
}
double te_bin_entropy(const TEBinnedSeries *bs) {
    if (!bs) return 0.0;
    double *hist = calloc((size_t)bs->n_bins, sizeof(double));
    if (!hist) return 0.0;
    for (int i = 0; i < bs->length; i++) hist[bs->bins[i]] += 1.0;
    double H = 0.0, inv_n = 1.0 / (double)bs->length;
    for (int b = 0; b < bs->n_bins; b++) if (hist[b] > 0.0) { double p = hist[b] * inv_n; H -= p * log(p); }
    free(hist);
    return H;
}

TEEmbedding* te_embed_create(const TEBinnedSeries *bs, int dim, int delay) {
    if (!bs || dim < 1 || dim > TE_MAX_EMBED || delay < 1) return NULL;
    int n_states = bs->length - (dim - 1) * delay;
    if (n_states < 10) return NULL;
    TEEmbedding *emb = calloc(1, sizeof(TEEmbedding));
    if (!emb) return NULL;
    emb->states = malloc((size_t)n_states * sizeof(int*));
    if (!emb->states) { free(emb); return NULL; }
    for (int i = 0; i < n_states; i++) {
        emb->states[i] = malloc((size_t)dim * sizeof(int));
        if (!emb->states[i]) { for (int j = 0; j < i; j++) free(emb->states[j]); free(emb->states); free(emb); return NULL; }
        for (int d = 0; d < dim; d++) emb->states[i][d] = bs->bins[i + d * delay];
    }
    emb->n_states = n_states; emb->embed_dim = dim; emb->delay = delay;
    return emb;
}
void te_embed_free(TEEmbedding *emb) {
    if (!emb) return;
    for (int i = 0; i < emb->n_states; i++) free(emb->states[i]);
    free(emb->states); free(emb);
}

TEJointProb* te_joint_estimate(const TEBinnedSeries *x, const TEBinnedSeries *y, int x_lag) {
    if (!x || !y || x_lag < 0) return NULL;
    int start = x_lag, n = x->length - start;
    if (n < 10) return NULL;
    TEJointProb *jp = calloc(1, sizeof(TEJointProb));
    if (!jp) return NULL;
    jp->nx = x->n_bins; jp->ny = y->n_bins; jp->n_samples = n;
    jp->joint = calloc((size_t)(jp->nx * jp->ny), sizeof(double));
    jp->marg_x = calloc((size_t)jp->nx, sizeof(double));
    jp->marg_y = calloc((size_t)jp->ny, sizeof(double));
    if (!jp->joint || !jp->marg_x || !jp->marg_y) { te_joint_free(jp); return NULL; }
    for (int i = 0; i < n; i++) {
        int xi = x->bins[i + start], yi = y->bins[i];
        jp->joint[xi * jp->ny + yi] += 1.0;
        jp->marg_x[xi] += 1.0; jp->marg_y[yi] += 1.0;
    }
    double inv_n = 1.0 / (double)n;
    for (int i = 0; i < jp->nx * jp->ny; i++) jp->joint[i] *= inv_n;
    for (int i = 0; i < jp->nx; i++) jp->marg_x[i] *= inv_n;
    for (int i = 0; i < jp->ny; i++) jp->marg_y[i] *= inv_n;
    return jp;
}
void te_joint_free(TEJointProb *jp) { if (jp) { free(jp->joint); free(jp->marg_x); free(jp->marg_y); free(jp); } }
double te_joint_entropy(const TEJointProb *jp) {
    if (!jp) return 0.0; double H = 0.0;
    for (int i = 0; i < jp->nx * jp->ny; i++) if (jp->joint[i] > TE_EPSILON) H -= jp->joint[i] * log(jp->joint[i]);
    return H;
}
double te_conditional_entropy(const TEJointProb *jp) {
    if (!jp) return 0.0; double H_xy = te_joint_entropy(jp), H_y = 0.0;
    for (int j = 0; j < jp->ny; j++) if (jp->marg_y[j] > TE_EPSILON) H_y -= jp->marg_y[j] * log(jp->marg_y[j]);
    return H_xy - H_y;
}

void te_result_print(const TEResult *r) {
    if (!r) { printf("TEResult: NULL\n"); return; }
    printf("=== Transfer Entropy Result ===\n");
    printf("  TE(Y->X):        %.6f nats\n", r->te_xy);
    printf("  TE(X->Y):        %.6f nats\n", r->te_yx);
    printf("  Net TE (Y-X):    %.6f nats\n", r->net_te);
    printf("  H(X_future):     %.6f\n", r->h_x_future);
    printf("  H(X_fut|X_past): %.6f\n", r->h_x_past);
    printf("  H(X_fut|X_past,Y_past): %.6f\n", r->h_x_past_y);
    printf("  Embed (k,l):     (%d,%d)\n", r->k_embed, r->l_embed);
    if (r->n_surrogates > 0) { printf("  Z-score: %.4f  P-value: %.6f  Significant: %s\n", r->z_score, r->p_value, r->significant ? "YES" : "no"); }
}
int te_result_is_significant(const TEResult *r, double alpha) { return (r && r->n_surrogates > 0 && r->p_value < alpha) ? 1 : 0; }

/* ?? CSV I/O ?????????????????????????????????????????????? */
int te_ts_load_csv(const char *filename, TETimeSeries **ts_out, int *n_series) {
    if (!filename || !ts_out || !n_series) return -1;
    FILE *fp = fopen(filename, "r"); if (!fp) return -1;
    char line[4096]; int col_count = 0, row_cap = 1024;
    double **cols = NULL; int *row_counts = NULL;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        if (col_count == 0) {
            char *tok = strtok(line, ",");
            while (tok) { col_count++; tok = strtok(NULL, ",\n\r"); }
            cols = calloc((size_t)col_count, sizeof(double*));
            row_counts = calloc((size_t)col_count, sizeof(int));
            if (!cols || !row_counts) { fclose(fp); return -1; }
            for (int c = 0; c < col_count; c++) { cols[c] = malloc((size_t)row_cap * sizeof(double)); row_counts[c] = 0; }
            continue;
        }
        char *tok = strtok(line, ","); int c = 0;
        while (tok && c < col_count) { double v = atof(tok); if (row_counts[c] < row_cap) cols[c][row_counts[c]++] = v; tok = strtok(NULL, ",\n\r"); c++; }
    }
    fclose(fp); *n_series = col_count;
    *ts_out = malloc((size_t)col_count * sizeof(TETimeSeries));
    for (int c = 0; c < col_count; c++) { char nm[32]; snprintf(nm, 32, "col_%d", c);
        (*ts_out)[c].data = cols[c]; (*ts_out)[c].length = row_counts[c]; (*ts_out)[c].name = strdup(nm); }
    free(cols); free(row_counts); return 0;
}
int te_ts_save_csv(const TETimeSeries *ts, const char *filename) {
    if (!ts || !filename) return -1; FILE *fp = fopen(filename, "w"); if (!fp) return -1;
    for (int i = 0; i < ts->length; i++) fprintf(fp, "%.10g\n", ts->data[i]);
    fclose(fp); return 0;
}

/* ?? Advanced Statistics ??????????????????????????????????? */
double te_ts_skewness(const TETimeSeries *ts) {
    if (!ts || ts->length < 3) return 0.0;
    double m = te_ts_mean(ts), s = te_ts_std(ts);
    if (s < TE_EPSILON) return 0.0; double sum = 0.0;
    for (int i = 0; i < ts->length; i++) { double z = (ts->data[i] - m) / s; sum += z * z * z; }
    return sum / (double)ts->length;
}
double te_ts_kurtosis(const TETimeSeries *ts) {
    if (!ts || ts->length < 4) return 0.0;
    double m = te_ts_mean(ts), s = te_ts_std(ts);
    if (s < TE_EPSILON) return 0.0; double sum = 0.0;
    for (int i = 0; i < ts->length; i++) { double z = (ts->data[i] - m) / s; sum += z * z * z * z; }
    return sum / (double)ts->length - 3.0;
}
double te_ts_autocorrelation(const TETimeSeries *ts, int lag) {
    if (!ts || lag < 0 || lag >= ts->length) return 0.0;
    double m = te_ts_mean(ts), num = 0.0, den = 0.0;
    for (int i = 0; i < ts->length - lag; i++) num += (ts->data[i] - m) * (ts->data[i + lag] - m);
    for (int i = 0; i < ts->length; i++) { double d = ts->data[i] - m; den += d * d; }
    return den > TE_EPSILON ? num / den : 0.0;
}
void te_ts_autocorrelation_function(const TETimeSeries *ts, int max_lag, double *acf) {
    if (!ts || !acf || max_lag < 1) return;
    for (int lag = 0; lag < max_lag && lag < ts->length; lag++) acf[lag] = te_ts_autocorrelation(ts, lag);
}
double te_ts_cross_correlation(const TETimeSeries *x, const TETimeSeries *y, int lag) {
    if (!x || !y || x->length != y->length) return 0.0;
    int n = x->length; double mx = te_ts_mean(x), my = te_ts_mean(y), sxy = 0.0, sx = 0.0, sy = 0.0;
    for (int i = 0; i < n; i++) { double dx = x->data[i] - mx, dy = y->data[i] - my; sx += dx * dx; sy += dy * dy; }
    if (sx < TE_EPSILON || sy < TE_EPSILON) return 0.0;
    for (int i = 0; i < n - lag; i++) { double dx = x->data[i + lag] - mx, dy = y->data[i] - my; sxy += dx * dy; }
    return sxy / sqrt(sx * sy);
}

/* ?? Optimal bin width (Freedman-Diaconis) ????????????????? */
int te_estimate_optimal_bins(const TETimeSeries *ts) {
    if (!ts || ts->length < 2) return 2;
    int n = ts->length; double *sorted = malloc((size_t)n * sizeof(double));
    if (!sorted) return 8; memcpy(sorted, ts->data, (size_t)n * sizeof(double));
    for (int i = 0; i < n - 1; i++) for (int j = i + 1; j < n; j++)
        if (sorted[i] > sorted[j]) { double t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }
    double q1 = sorted[n / 4], q3 = sorted[3 * n / 4], iqr = q3 - q1;
    if (iqr < TE_EPSILON) iqr = te_ts_std(ts);
    double bw = 2.0 * iqr / pow((double)n, 1.0 / 3.0);
    double rng = te_ts_max(ts) - te_ts_min(ts);
    int nb = (int)(rng / bw) + 1; if (nb < 2) nb = 2; if (nb > TE_MAX_BINS) nb = TE_MAX_BINS;
    free(sorted); return nb;
}

/* ?? Mutual Information ???????????????????????????????????? */
double te_mutual_information(const TETimeSeries *x, const TETimeSeries *y, int n_bins, int lag) {
    if (!x || !y || n_bins < 2) return 0.0;
    TEBinnedSeries *bx = te_bin_create(x, n_bins, TE_BIN_EQUAL_WIDTH);
    TEBinnedSeries *by = te_bin_create(y, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bx || !by) { te_bin_free(bx); te_bin_free(by); return 0.0; }
    TEJointProb *jp = te_joint_estimate(bx, by, lag);
    if (!jp) { te_bin_free(bx); te_bin_free(by); return 0.0; }
    double Hx = 0.0, Hy = 0.0;
    for (int i = 0; i < jp->nx; i++) if (jp->marg_x[i] > TE_EPSILON) Hx -= jp->marg_x[i] * log(jp->marg_x[i]);
    for (int j = 0; j < jp->ny; j++) if (jp->marg_y[j] > TE_EPSILON) Hy -= jp->marg_y[j] * log(jp->marg_y[j]);
    double Hxy = te_joint_entropy(jp);
    te_joint_free(jp); te_bin_free(bx); te_bin_free(by);
    return Hx + Hy - Hxy;
}

/* ?? Active Information Storage (Lizier et al., 2012) ??????
 * AIS = I(X_{n+1}; X_n^{(k)}) = H(X_{n+1}) - H(X_{n+1}|X_n^{(k)})
 */
double te_active_information_storage(const TETimeSeries *x, int k, int n_bins) {
    if (!x || k < 1) return 0.0;
    TEBinnedSeries *bx = te_bin_create(x, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bx) return 0.0;
    double H_xf = te_bin_entropy(bx);
    int nb = bx->n_bins, n_states = 1; for (int i = 0; i < k; i++) n_states *= nb;
    int n = bx->length - k; if (n < 10) { te_bin_free(bx); return 0.0; }
    size_t sz = (size_t)nb * (size_t)n_states;
    double *joint = calloc(sz, sizeof(double));
    double *marg = calloc((size_t)n_states, sizeof(double));
    if (!joint || !marg) { free(joint); free(marg); te_bin_free(bx); return 0.0; }
    double inv_n = 1.0 / (double)n;
    for (int t = 0; t < n; t++) {
        int xf = bx->bins[t + k], xp = 0, mult = 1;
        for (int d = 0; d < k; d++) { xp += bx->bins[t + d] * mult; mult *= nb; }
        joint[(size_t)xf * (size_t)n_states + (size_t)xp] += inv_n;
        marg[xp] += inv_n;
    }
    double H_xf_xp = 0.0, H_xp = 0.0;
    for (size_t i = 0; i < sz; i++) if (joint[i] > TE_EPSILON) H_xf_xp -= joint[i] * log(joint[i]);
    for (int i = 0; i < n_states; i++) if (marg[i] > TE_EPSILON) H_xp -= marg[i] * log(marg[i]);
    free(joint); free(marg); te_bin_free(bx);
    return H_xf + H_xp - H_xf_xp;
}

/* ?? Predictive Information (Bialek et al., 2001) ??????????
 * PI = I(X_{future}^{h}; X_{past}^{k})
 */
double te_predictive_information(const TETimeSeries *x, int k, int h, int n_bins) {
    if (!x || k < 1 || h < 1) return 0.0;
    TEBinnedSeries *bx = te_bin_create(x, n_bins, TE_BIN_EQUAL_WIDTH);
    if (!bx) return 0.0;
    int n = bx->length - k - h, nb = bx->n_bins;
    if (n < 10) { te_bin_free(bx); return 0.0; }
    int nsf = 1; for (int i = 0; i < h; i++) nsf *= nb;
    int nsp = 1; for (int i = 0; i < k; i++) nsp *= nb;
    size_t sz = (size_t)nsf * (size_t)nsp;
    double *joint = calloc(sz, sizeof(double)), *mf = calloc((size_t)nsf, sizeof(double)), *mp = calloc((size_t)nsp, sizeof(double));
    if (!joint || !mf || !mp) { free(joint); free(mf); free(mp); te_bin_free(bx); return 0.0; }
    double inv_n = 1.0 / (double)n;
    for (int t = 0; t < n; t++) {
        int cf = 0, mult = 1; for (int d = 0; d < h; d++) { cf += bx->bins[t + k + d] * mult; mult *= nb; }
        int cp = 0; mult = 1; for (int d = 0; d < k; d++) { cp += bx->bins[t + d] * mult; mult *= nb; }
        joint[(size_t)cf * (size_t)nsp + (size_t)cp] += inv_n; mf[cf] += inv_n; mp[cp] += inv_n;
    }
    double Hf = 0, Hp = 0, Hfp = 0;
    for (int i = 0; i < nsf; i++) if (mf[i] > TE_EPSILON) Hf -= mf[i] * log(mf[i]);
    for (int i = 0; i < nsp; i++) if (mp[i] > TE_EPSILON) Hp -= mp[i] * log(mp[i]);
    for (size_t i = 0; i < sz; i++) if (joint[i] > TE_EPSILON) Hfp -= joint[i] * log(joint[i]);
    free(joint); free(mf); free(mp); te_bin_free(bx);
    return Hf + Hp - Hfp;
}

/* ???????????????????????????????????????????????????????????
 * Extended Time Series Operations
 * ??????????????????????????????????????????????????????????? */
void te_ts_detrend_linear(TETimeSeries *ts) {
    if (!ts || ts->length < 2) return;
    double mx = (double)(ts->length - 1) / 2.0, my = te_ts_mean(ts);
    double sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < ts->length; i++) { double dx = (double)i - mx; sxx += dx * dx; sxy += dx * (ts->data[i] - my); }
    double slope = sxy / sxx;
    for (int i = 0; i < ts->length; i++) ts->data[i] -= slope * (double)i;
}
void te_ts_moving_average(TETimeSeries *ts, int window) {
    if (!ts || window < 2 || window >= ts->length) return;
    double *smoothed = malloc((size_t)ts->length * sizeof(double));
    if (!smoothed) return;
    for (int i = 0; i < ts->length; i++) {
        int lo = i - window / 2, hi = i + window / 2;
        if (lo < 0) lo = 0; if (hi >= ts->length) hi = ts->length - 1;
        double sum = 0.0; for (int j = lo; j <= hi; j++) sum += ts->data[j];
        smoothed[i] = sum / (double)(hi - lo + 1);
    }
    memcpy(ts->data, smoothed, (size_t)ts->length * sizeof(double));
    free(smoothed);
}
void te_ts_exponential_smooth(TETimeSeries *ts, double alpha) {
    if (!ts || ts->length < 2 || alpha <= 0.0 || alpha > 1.0) return;
    double s = ts->data[0];
    for (int i = 1; i < ts->length; i++) { s = alpha * ts->data[i] + (1.0 - alpha) * s; ts->data[i] = s; }
}
void te_ts_add_noise(TETimeSeries *ts, double sigma) {
    if (!ts || sigma <= 0.0) return;
    for (int i = 0; i < ts->length; i++) {
        double u1 = (double)rand() / (double)RAND_MAX, u2 = (double)rand() / (double)RAND_MAX;
        double z = sqrt(-2.0 * log(u1 + 1e-15)) * cos(6.283185307179586 * u2);
        ts->data[i] += sigma * z;
    }
}
void te_ts_resample(TETimeSeries *ts, int new_length) {
    if (!ts || new_length < 2 || new_length == ts->length) return;
    double *new_data = malloc((size_t)new_length * sizeof(double));
    if (!new_data) return;
    double ratio = (double)(ts->length - 1) / (double)(new_length - 1);
    for (int i = 0; i < new_length; i++) {
        double pos = (double)i * ratio; int lo = (int)pos; double frac = pos - (double)lo;
        if (lo >= ts->length - 1) new_data[i] = ts->data[ts->length - 1];
        else new_data[i] = ts->data[lo] * (1.0 - frac) + ts->data[lo + 1] * frac;
    }
    free(ts->data); ts->data = new_data; ts->length = new_length;
}
double te_ts_quantile(const TETimeSeries *ts, double q) {
    if (!ts || ts->length < 1 || q < 0.0 || q > 1.0) return 0.0;
    double *sorted = malloc((size_t)ts->length * sizeof(double));
    if (!sorted) return 0.0;
    memcpy(sorted, ts->data, (size_t)ts->length * sizeof(double));
    for (int i = 0; i < ts->length - 1; i++) for (int j = i + 1; j < ts->length; j++)
        if (sorted[i] > sorted[j]) { double t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }
    int idx = (int)(q * (double)(ts->length - 1) + 0.5);
    double val = sorted[idx]; free(sorted); return val;
}

/* ???????????????????????????????????????????????????????????
 * Multi-variate binning for joint entropy estimation
 * ??????????????????????????????????????????????????????????? */
TEBinnedSeries* te_bin_create_adaptive(const TETimeSeries *ts, int min_bins, int max_bins) {
    int optimal = te_estimate_optimal_bins(ts);
    if (optimal < min_bins) optimal = min_bins; if (optimal > max_bins) optimal = max_bins;
    return te_bin_create(ts, optimal, TE_BIN_EQUAL_WIDTH);
}
int te_bin_index_2d(const TEBinnedSeries *bx, const TEBinnedSeries *by, int idx) {
    if (!bx || !by || idx < 0 || idx >= bx->length || idx >= by->length) return -1;
    return bx->bins[idx] * by->n_bins + by->bins[idx];
}
double te_joint_entropy_2d(const TEBinnedSeries *bx, const TEBinnedSeries *by) {
    if (!bx || !by || bx->length != by->length) return 0.0;
    int nb = bx->n_bins * by->n_bins, n = bx->length;
    int *counts = calloc((size_t)nb, sizeof(int));
    if (!counts) return 0.0;
    for (int i = 0; i < n; i++) counts[bx->bins[i] * by->n_bins + by->bins[i]]++;
    double H = 0.0, inv_n = 1.0 / (double)n;
    for (int i = 0; i < nb; i++) if (counts[i] > 0) { double p = (double)counts[i] * inv_n; H -= p * log(p); }
    free(counts); return H;
}
double te_conditional_entropy_2d(const TEBinnedSeries *bx, const TEBinnedSeries *by) {
    double H_xy = te_joint_entropy_2d(bx, by), H_y = te_bin_entropy(by);
    return H_xy - H_y;
}

/* Serialization to/from binary format for efficient storage */
int te_ts_save_binary(const TETimeSeries *ts, const char *filename) {
    if (!ts || !filename) return -1;
    FILE *fp = fopen(filename, "wb"); if (!fp) return -1;
    fwrite(&ts->length, sizeof(int), 1, fp);
    fwrite(ts->data, sizeof(double), (size_t)ts->length, fp);
    fclose(fp); return 0;
}
int te_ts_load_binary(const char *filename, TETimeSeries **ts) {
    if (!filename || !ts) return -1;
    FILE *fp = fopen(filename, "rb"); if (!fp) return -1;
    int len; if (fread(&len, sizeof(int), 1, fp) != 1) { fclose(fp); return -1; }
    if (len < 1 || len > 10000000) { fclose(fp); return -1; }
    double *buf = malloc((size_t)len * sizeof(double));
    if (!buf) { fclose(fp); return -1; }
    if (fread(buf, sizeof(double), (size_t)len, fp) != (size_t)len) { free(buf); fclose(fp); return -1; }
    fclose(fp);
    *ts = te_ts_create(buf, len, filename); free(buf);
    return (*ts) ? 0 : -1;
}

/* Generate test signals */
TETimeSeries* te_ts_generate_sine(int n, double freq, double amp, double phase) {
    double *buf = malloc((size_t)n * sizeof(double));
    if (!buf) return NULL;
    for (int i = 0; i < n; i++) buf[i] = amp * sin(6.283185307179586 * freq * (double)i / (double)n + phase);
    TETimeSeries *ts = te_ts_create(buf, n, "sine"); free(buf); return ts;
}
TETimeSeries* te_ts_generate_henon(int n, double a, double b) {
    double *buf = malloc((size_t)n * sizeof(double));
    if (!buf) return NULL;
    double x = 0.0, y = 0.0;
    for (int i = 0; i < n; i++) {
        buf[i] = x;
        double xn = 1.0 - a * x * x + y;
        y = b * x; x = xn;
    }
    TETimeSeries *ts = te_ts_create(buf, n, "henon"); free(buf); return ts;
}
TETimeSeries* te_ts_generate_logistic(int n, double r, double x0) {
    double *buf = malloc((size_t)n * sizeof(double));
    if (!buf) return NULL;
    double x = x0;
    for (int i = 0; i < n; i++) { buf[i] = x; x = r * x * (1.0 - x); }
    TETimeSeries *ts = te_ts_create(buf, n, "logistic"); free(buf); return ts;
}
TETimeSeries* te_ts_generate_ar1(int n, double phi, double sigma) {
    double *buf = malloc((size_t)n * sizeof(double));
    if (!buf) return NULL;
    buf[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double u1 = (double)rand() / (double)RAND_MAX, u2 = (double)rand() / (double)RAND_MAX;
        buf[i] = phi * buf[i-1] + sigma * sqrt(-2.0*log(u1+1e-15)) * cos(6.283185307179586*u2);
    }
    TETimeSeries *ts = te_ts_create(buf, n, "ar1"); free(buf); return ts;
}
TETimeSeries* te_ts_generate_coupled_ar1(int n, double c, double phi, double sigma) {
    double *buf_x = malloc((size_t)n * sizeof(double));
    double *buf_y = malloc((size_t)n * sizeof(double));
    if (!buf_x || !buf_y) { free(buf_x); free(buf_y); return NULL; }
    double x = 0.0, y = 0.0;
    for (int i = 0; i < n; i++) {
        buf_x[i] = x; buf_y[i] = y;
        double u1 = (double)rand()/(double)RAND_MAX, u2 = (double)rand()/(double)RAND_MAX;
        double ex = sigma*sqrt(-2.0*log(u1+1e-15))*cos(6.283185307179586*u2);
        u1 = (double)rand()/(double)RAND_MAX; u2 = (double)rand()/(double)RAND_MAX;
        double ey = sigma*sqrt(-2.0*log(u1+1e-15))*cos(6.283185307179586*u2);
        double xn = phi*x + c*y + ex;
        y = phi*y + ey; x = xn;
    }
    TETimeSeries *ts = te_ts_create(buf_x, n, "coupled_x");
    free(buf_x); free(buf_y); return ts;
}

/* Convert TEResult to JSON-like string */
void te_result_to_string(const TEResult *r, char *buf, int buf_size) {
    if (!r || !buf || buf_size < 1) return;
    snprintf(buf, (size_t)buf_size,
        "{\"te_xy\":%.6f,\"te_yx\":%.6f,\"net_te\":%.6f,\"k\":%d,\"l\":%d,\"p_value\":%.6f,\"significant\":%s}",
        r->te_xy, r->te_yx, r->net_te, r->k_embed, r->l_embed, r->p_value, r->significant ? "true" : "false");
}

/* Lempel-Ziv complexity (normalized) */
double te_lz_complexity(const TEBinnedSeries *bs) {
    if (!bs || bs->length < 2) return 0.0;
    int n = bs->length, *s = bs->bins, c = 1, l = 1, i = 0, j = 1, k = 1;
    while (i + l < n) {
        if (s[i + k - 1] == s[j + k - 1]) {
            k++;
            if (j + k > n) { c++; break; }
            if (i + k > j) { c++; i = j; j++; if (j >= n) break; k = 1; l = 1; }
        } else { k = 1; l++; j++; if (j + l > n) { c++; break; } }
    }
    double norm = (double)n / log2((double)n + 1e-15);
    return (double)c / norm;
}

/* ?? Kernel Density Estimation for continuous entropy ?????? */
double te_kde_entropy(const double *data, int n, double bw) {
    if (!data || n < 2 || bw <= 0.0) return 0.0;
    double H = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) { double d = (data[i] - data[j]) / bw; sum += exp(-0.5 * d * d); }
        double p = sum / ((double)n * bw * 2.5066282746310002);
        if (p > TE_EPSILON) H -= log(p);
    }
    return H / (double)n;
}
double te_kde_mutual_information(const double *x, const double *y, int n, double bw_x, double bw_y) {
    if (!x || !y || n < 2) return 0.0;
    double Hx = te_kde_entropy(x, n, bw_x), Hy = te_kde_entropy(y, n, bw_y);
    double *joint = malloc((size_t)(2 * n) * sizeof(double));
    if (!joint) return 0.0;
    for (int i = 0; i < n; i++) { joint[i*2] = x[i]; joint[i*2+1] = y[i]; }
    double Hxy = te_kde_entropy(joint, n, sqrt(bw_x * bw_y));
    free(joint);
    double mi = Hx + Hy - Hxy;
    return mi > 0.0 ? mi : 0.0;
}

/* ?? Time series downsampling and aggregation ?????????????? */
TETimeSeries* te_ts_aggregate(const TETimeSeries *ts, int factor, int mode) {
    if (!ts || factor < 2 || ts->length < factor) return NULL;
    int n_out = ts->length / factor; if (n_out < 2) return NULL;
    double *buf = malloc((size_t)n_out * sizeof(double)); if (!buf) return NULL;
    for (int i = 0; i < n_out; i++) {
        double sum = 0.0, mn = 1e100, mx = -1e100;
        for (int j = 0; j < factor; j++) { double v = ts->data[i*factor + j]; sum += v; if (v < mn) mn = v; if (v > mx) mx = v; }
        if (mode == 0) buf[i] = sum / (double)factor;
        else if (mode == 1) buf[i] = mn;
        else buf[i] = mx;
    }
    TETimeSeries *out = te_ts_create(buf, n_out, "aggregated"); free(buf); return out;
}

/* ?? Sliding window statistics ????????????????????????????? */
void te_ts_sliding_stats(const TETimeSeries *ts, int window, double *mean_out, double *var_out, double *skew_out) {
    if (!ts || window < 2 || window > ts->length) return;
    int n_out = ts->length - window + 1;
    for (int i = 0; i < n_out; i++) {
        double sum = 0.0, sum2 = 0.0, sum3 = 0.0;
        for (int j = 0; j < window; j++) { double v = ts->data[i+j]; sum += v; sum2 += v*v; sum3 += v*v*v; }
        double m = sum / (double)window, v = sum2/(double)window - m*m;
        if (mean_out) mean_out[i] = m;
        if (var_out) var_out[i] = v;
        if (skew_out && v > TE_EPSILON) { double s3 = sum3/(double)window - 3*m*sum2/(double)window + 2*m*m*m; skew_out[i] = s3 / pow(v, 1.5); }
    }
}

/* ?? Histogram-based probability density ??????????????????? */
int te_ts_histogram(const TETimeSeries *ts, int n_bins, double *edges, double *counts) {
    if (!ts || !edges || !counts || n_bins < 2) return -1;
    double mn = te_ts_min(ts), mx = te_ts_max(ts), rng = mx - mn;
    if (rng < TE_EPSILON) rng = 1.0;
    for (int i = 0; i <= n_bins; i++) edges[i] = mn + (double)i * rng / (double)n_bins;
    memset(counts, 0, (size_t)n_bins * sizeof(double));
    for (int i = 0; i < ts->length; i++) {
        int b = (int)((double)n_bins * (ts->data[i] - mn) / rng);
        if (b < 0) b = 0; if (b >= n_bins) b = n_bins - 1;
        counts[b] += 1.0;
    }
    return 0;
}

/* ?? Surrogate generation for time series ?????????????????? */
void te_ts_shuffle(TETimeSeries *ts) { if (!ts) return; for (int i = ts->length-1; i>0; i--) { int j = rand()%(i+1); double t = ts->data[i]; ts->data[i]=ts->data[j]; ts->data[j]=t; } }
TETimeSeries* te_ts_clone(const TETimeSeries *ts) { return ts ? te_ts_create(ts->data, ts->length, ts->name) : NULL; }

/* ?? Cumulative distribution function ?????????????????????? */
double te_ts_empirical_cdf(const TETimeSeries *ts, double x) {
    if (!ts || ts->length < 1) return 0.0;
    int count = 0;
    for (int i = 0; i < ts->length; i++) if (ts->data[i] <= x) count++;
    return (double)count / (double)ts->length;
}

/* ?? Percentile computation ???????????????????????????????? */
double te_ts_percentile(const TETimeSeries *ts, double pct) {
    if (!ts || pct < 0.0 || pct > 100.0) return 0.0;
    return te_ts_quantile(ts, pct / 100.0);
}

/* ?? Root mean square ?????????????????????????????????????? */
double te_ts_rms(const TETimeSeries *ts) {
    if (!ts || ts->length < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < ts->length; i++) sum += ts->data[i] * ts->data[i];
    return sqrt(sum / (double)ts->length);
}

/* ?? Zero-crossing rate ???????????????????????????????????? */
int te_ts_zero_crossings(const TETimeSeries *ts) {
    if (!ts || ts->length < 2) return 0;
    int count = 0;
    for (int i = 1; i < ts->length; i++) if (ts->data[i-1] * ts->data[i] < 0.0) count++;
    return count;
}

/* ?? Hjorth parameters: activity, mobility, complexity ????? */
void te_ts_hjorth(const TETimeSeries *ts, double *activity, double *mobility, double *complexity) {
    if (!ts || !activity || !mobility || !complexity) return;
    *activity = te_ts_variance(ts);
    TETimeSeries *diff = te_ts_create(ts->data, ts->length, "diff"); te_ts_diff(diff);
    double var_d1 = te_ts_variance(diff);
    *mobility = (*activity > TE_EPSILON) ? sqrt(var_d1 / *activity) : 0.0;
    TETimeSeries *diff2 = te_ts_create(diff->data, diff->length, "diff2"); te_ts_diff(diff2);
    double var_d2 = te_ts_variance(diff2);
    *complexity = (*mobility > TE_EPSILON) ? sqrt(var_d2/var_d1) / *mobility : 0.0;
    te_ts_free(diff); te_ts_free(diff2);
}

/* ?? Approximate Entropy (ApEn) ???????????????????????????? */
double te_approximate_entropy(const TETimeSeries *ts, int m, double r) {
    if (!ts || m < 1 || r <= 0.0 || ts->length < m + 2) return 0.0;
    int n = ts->length;
    double apen = 0.0;
    for (int k = 0; k < 2; k++) {
        int mm = m + k, nn = n - mm + 1; if (nn < 2) continue;
        double *c = calloc((size_t)nn, sizeof(double)); if (!c) return 0.0;
        for (int i = 0; i < nn; i++) {
            for (int j = 0; j < nn; j++) {
                if (i == j) continue;
                double max_diff = 0.0;
                for (int d = 0; d < mm; d++) { double diff = fabs(ts->data[i+d] - ts->data[j+d]); if (diff > max_diff) max_diff = diff; }
                if (max_diff <= r) c[i] += 1.0;
            }
            c[i] /= (double)(nn - 1);
        }
        double phi = 0.0;
        for (int i = 0; i < nn; i++) if (c[i] > TE_EPSILON) phi += log(c[i]);
        phi /= (double)nn;
        if (k == 0) apen = phi; else apen -= phi;
        free(c);
    }
    return apen;
}

/* ?? Sample Entropy (SampEn) ??????????????????????????????? */
double te_sample_entropy(const TETimeSeries *ts, int m, double r) {
    if (!ts || m < 1 || r <= 0.0 || ts->length < m + 2) return 0.0;
    int n = ts->length, nn = n - m;
    if (nn < 2) return 0.0;
    int A = 0, B = 0;
    for (int i = 0; i < nn; i++) {
        for (int j = i + 1; j < nn; j++) {
            double max_d = 0.0; int match = 1;
            for (int d = 0; d < m; d++) { double diff = fabs(ts->data[i+d]-ts->data[j+d]); if (diff > r) { match = 0; break; } if (diff > max_d) max_d = diff; }
            if (match) { B++;
                if (fabs(ts->data[i+m]-ts->data[j+m]) <= r) A++; }
        }
    }
    if (B == 0) return 0.0;
    return -log((double)A / (double)B);
}

/* ?? Durbin-Watson statistic for autocorrelation ?????????? */
double te_ts_durbin_watson(const TETimeSeries *ts) {
    if (!ts || ts->length < 3) return 2.0;
    double num=0.0, den=0.0;
    for (int i=1;i<ts->length;i++) { double d=ts->data[i]-ts->data[i-1]; num+=d*d; }
    for (int i=0;i<ts->length;i++) { double d=ts->data[i]-te_ts_mean(ts); den+=d*d; }
    return den>TE_EPSILON ? num/den : 2.0;
}

/* ?? Augmented Dickey-Fuller test statistic ??????????????? */
double te_ts_adf_statistic(const TETimeSeries *ts, int max_lag) {
    if (!ts || ts->length < max_lag+3) return 0.0;
    int n = ts->length - 1;
    double *dy = malloc((size_t)n*sizeof(double));
    if (!dy) return 0.0;
    for (int i=0;i<n;i++) dy[i]=ts->data[i+1]-ts->data[i];
    double sum_xy=0, sum_xx=0;
    for (int i=0;i<n;i++) { sum_xy+=dy[i]*(ts->data[i]-te_ts_mean(ts)); sum_xx+=(ts->data[i]-te_ts_mean(ts))*(ts->data[i]-te_ts_mean(ts)); }
    double rho = sum_xx>TE_EPSILON ? sum_xy/sum_xx : 0.0;
    double se = 0.0; for (int i=0;i<n;i++) { double e=dy[i]-rho*ts->data[i]; se+=e*e; }
    se = sqrt(se/(double)(n-2));
    free(dy);
    return se>TE_EPSILON ? rho/se : 0.0;
}

/* ?? Hurst exponent via R/S analysis ?????????????????????? */
double te_ts_hurst_rs(const TETimeSeries *ts, int min_window, int max_window, int n_windows) {
    if (!ts || min_window<4 || max_window<min_window || n_windows<2) return 0.5;
    double *log_n = malloc((size_t)n_windows*sizeof(double));
    double *log_rs = malloc((size_t)n_windows*sizeof(double));
    if (!log_n||!log_rs) { free(log_n);free(log_rs); return 0.5; }
    for (int w=0;w<n_windows;w++) {
        int win = min_window + (max_window-min_window)*w/(n_windows-1);
        if (win>=ts->length) win=ts->length-1;
        int n_seg = ts->length/win; if (n_seg<1) { log_n[w]=log((double)win); log_rs[w]=0.0; continue; }
        double rs_sum=0.0;
        for (int s=0;s<n_seg;s++) {
            double m=0; for (int i=s*win;i<(s+1)*win;i++) m+=ts->data[i]; m/=win;
            double *dev=malloc((size_t)win*sizeof(double));
            if (!dev) continue;
            dev[0]=ts->data[s*win]-m; for (int i=1;i<win;i++) dev[i]=dev[i-1]+ts->data[s*win+i]-m;
            double r=dev[0],mn_r=dev[0],mx_r=dev[0];
            for (int i=1;i<win;i++) { if (dev[i]>mx_r) mx_r=dev[i]; if (dev[i]<mn_r) mn_r=dev[i]; }
            r = mx_r-mn_r;
            double var=0; for (int i=s*win;i<(s+1)*win;i++) { double d=ts->data[i]-m; var+=d*d; }
            double s_val = sqrt(var/(double)win);
            if (s_val>TE_EPSILON) rs_sum+=r/s_val;
            free(dev);
        }
        log_n[w]=log((double)win); log_rs[w]=log(rs_sum/(double)n_seg+1e-15);
    }
    double mx_n=0; for (int i=0;i<n_windows;i++) mx_n+=log_n[i]; mx_n/=n_windows;
    double mx_rs=0; for (int i=0;i<n_windows;i++) mx_rs+=log_rs[i]; mx_rs/=n_windows;
    double num=0,den=0;
    for (int i=0;i<n_windows;i++) { double dn=log_n[i]-mx_n; num+=dn*(log_rs[i]-mx_rs); den+=dn*dn; }
    free(log_n);free(log_rs);
    return den>TE_EPSILON ? num/den : 0.5;
}

/* ?? Time series stationarity check (KPSS-like) ??????????? */
double te_ts_kpss_stat(const TETimeSeries *ts) {
    if (!ts || ts->length<10) return 0.0;
    double m=te_ts_mean(ts), *cum=malloc((size_t)ts->length*sizeof(double));
    if (!cum) return 0.0;
    cum[0]=ts->data[0]-m; for (int i=1;i<ts->length;i++) cum[i]=cum[i-1]+ts->data[i]-m;
    double s2=0; for (int i=0;i<ts->length;i++) s2+=cum[i]*cum[i];
    s2/=(double)(ts->length*ts->length);
    double var=te_ts_variance(ts);
    free(cum);
    return var>TE_EPSILON ? s2/var : 0.0;
}

/* ?? Phase randomization surrogate ???????????????????????? */
void te_ts_phase_randomize(const TETimeSeries *ts, double *surrogate) {
    if (!ts || !surrogate || ts->length<4) return;
    int n=ts->length, n2=n/2;
    double *re=malloc((size_t)n*sizeof(double)), *im=calloc((size_t)n,sizeof(double));
    if(!re||!im) { free(re);free(im); return; }
    for(int k=0;k<n;k++) { re[k]=0; for(int j=0;j<n;j++) re[k]+=ts->data[j]*cos(-6.283185307179586*(double)(k*j)/(double)n); im[k]=0; for(int j=0;j<n;j++) im[k]+=ts->data[j]*sin(-6.283185307179586*(double)(k*j)/(double)n); }
    for(int k=1;k<n2;k++) { double mag=sqrt(re[k]*re[k]+im[k]*im[k]), phase=6.283185307179586*(double)rand()/(double)RAND_MAX; re[k]=mag*cos(phase); im[k]=mag*sin(phase); re[n-k]=re[k]; im[n-k]=-im[k]; }
    im[0]=0; if(n%2==0) im[n2]=0;
    for(int j=0;j<n;j++) { surrogate[j]=0; for(int k=0;k<n;k++) surrogate[j]+=re[k]*cos(6.283185307179586*(double)(k*j)/(double)n)-im[k]*sin(6.283185307179586*(double)(k*j)/(double)n); surrogate[j]/=(double)n; }
    free(re);free(im);
}

/* ???????????????????????????????????????????????????????????
 * Extended Information-Theoretic Measures
 * ??????????????????????????????????????????????????????????? */

/* ?? Jensen-Shannon Divergence ????????????????????????????? */
double te_js_divergence(const double *p, const double *q, int n) {
    if (!p || !q || n < 1) return -1.0; double *m=malloc((size_t)n*sizeof(double));
    if (!m) return -1.0; for (int i=0;i<n;i++) m[i]=0.5*(p[i]+q[i]);
    double kl_pm=0, kl_qm=0;
    for (int i=0;i<n;i++) { if(p[i]>TE_EPSILON&&m[i]>TE_EPSILON) kl_pm+=p[i]*log(p[i]/m[i]); if(q[i]>TE_EPSILON&&m[i]>TE_EPSILON) kl_qm+=q[i]*log(q[i]/m[i]); }
    free(m); return 0.5*kl_pm+0.5*kl_qm;
}

/* ?? R?nyi entropy of order alpha ????????????????????????? */
double te_renyi_entropy(const TEBinnedSeries *bs, double alpha) {
    if (!bs) return 0.0; double *hist=calloc((size_t)bs->n_bins,sizeof(double));
    if (!hist) return 0.0; for (int i=0;i<bs->length;i++) hist[bs->bins[i]]+=1.0;
    double inv_n=1.0/(double)bs->length, sum=0.0;
    for (int b=0;b<bs->n_bins;b++) if (hist[b]>0) { double p=hist[b]*inv_n; sum+=pow(p,alpha); }
    free(hist);
    if (fabs(alpha-1.0)<TE_EPSILON) return te_bin_entropy(bs);
    return log(sum)/(1.0-alpha);
}

/* ?? Tsallis entropy (non-extensive) ?????????????????????? */
double te_tsallis_entropy(const TEBinnedSeries *bs, double q) {
    if (!bs) return 0.0; double *hist=calloc((size_t)bs->n_bins,sizeof(double));
    if (!hist) return 0.0; for (int i=0;i<bs->length;i++) hist[bs->bins[i]]+=1.0;
    double inv_n=1.0/(double)bs->length, sum=0.0;
    for (int b=0;b<bs->n_bins;b++) if (hist[b]>0) { double p=hist[b]*inv_n; sum+=pow(p,q); }
    free(hist);
    if (fabs(q-1.0)<TE_EPSILON) return te_bin_entropy(bs);
    return (1.0-sum)/(q-1.0);
}

/* ?? Effective measure complexity (Grassberger) ???????????? */
double te_effective_measure_complexity(const TETimeSeries *ts, int max_k, int n_bins) {
    if (!ts || max_k<1) return 0.0;
    TEBinnedSeries *bs=te_bin_create(ts,n_bins,TE_BIN_EQUAL_WIDTH);
    if (!bs) return 0.0; double emc=0.0, h1=te_bin_entropy(bs), h_inf=1e100;
    for (int k=1;k<=max_k;k++) {
        double hk=te_schreiber_entropy_rate(bs,k);
        emc+=h1-hk; if (hk<h_inf) h_inf=hk;
    }
    te_bin_free(bs); return emc;
}

/* ?? Excess entropy (predictive information limit) ???????? */
double te_excess_entropy(const TETimeSeries *ts, int max_k, int n_bins) {
    if (!ts||max_k<1) return 0.0;
    TEBinnedSeries *bs=te_bin_create(ts,n_bins,TE_BIN_EQUAL_WIDTH);
    if (!bs) return 0.0; double h=te_bin_entropy(bs), h_inf=1e100;
    for (int k=1;k<=max_k;k++) { double hk=te_schreiber_entropy_rate(bs,k); if (hk<h_inf) h_inf=hk; }
    te_bin_free(bs); return h-h_inf;
}

/* ?? Lempel-Ziv complexity (1976) ?????????????????????????? */
int te_lz76_complexity(const int *seq, int n) {
    if (!seq||n<1) return 0; int max_len=256, *dict=malloc((size_t)(max_len*max_len)*sizeof(int));
    if (!dict) return 0; int dict_size=0, i=0, c=1;
    while (i<n) { int len=1; while (i+len<=n) { int found=0;
        for (int d=0;d<dict_size;d++) { int match=1;
            for (int k=0;k<len&&dict[d*max_len+k]>=0;k++) if (dict[d*max_len+k]!=seq[i+k]){match=0;break;}
            if (match){found=1;break;}
        }
        if (found) len++; else break;
    }
    if (i+len<=n&&dict_size<max_len) { for (int k=0;k<len;k++) dict[dict_size*max_len+k]=seq[i+k]; for (int k=len;k<max_len;k++) dict[dict_size*max_len+k]=-1; dict_size++; }
    i+=len; c++;
    }
    free(dict); return c;
}

/* ?? Quantile-Quantile (Q-Q) plot data ???????????????????? */
int te_ts_qq_data(const TETimeSeries *a, const TETimeSeries *b, double *qa, double *qb, int n_quantiles) {
    if (!a||!b||!qa||!qb||n_quantiles<2) return -1;
    double *sa=malloc((size_t)a->length*sizeof(double)),*sb=malloc((size_t)b->length*sizeof(double));
    if(!sa||!sb){free(sa);free(sb);return -1;}
    memcpy(sa,a->data,(size_t)a->length*sizeof(double));memcpy(sb,b->data,(size_t)b->length*sizeof(double));
    for(int i=0;i<a->length-1;i++)for(int j=i+1;j<a->length;j++){if(sa[i]>sa[j]){double t=sa[i];sa[i]=sa[j];sa[j]=t;}}
    for(int i=0;i<b->length-1;i++)for(int j=i+1;j<b->length;j++){if(sb[i]>sb[j]){double t=sb[i];sb[i]=sb[j];sb[j]=t;}}
    for(int k=0;k<n_quantiles;k++){double q=(double)(k+1)/(double)(n_quantiles+1);
        int ia=(int)(q*(double)(a->length-1)),ib=(int)(q*(double)(b->length-1));
        qa[k]=sa[ia];qb[k]=sb[ib];}
    free(sa);free(sb);return 0;
}

/* ?? Time series segment statistics ??????????????????????? */
void te_ts_segment_stats(const TETimeSeries *ts, int n_segments, double *means, double *vars,
                          double *mins, double *maxs) {
    if(!ts||n_segments<1) return;
    int seg_len=ts->length/n_segments; if(seg_len<1) return;
    for(int s=0;s<n_segments;s++){int start=s*seg_len,end=(s==n_segments-1)?ts->length:start+seg_len;
        double sum=0,sq=0,mn=1e100,mx=-1e100;
        for(int i=start;i<end;i++){double v=ts->data[i];sum+=v;sq+=v*v;if(v<mn)mn=v;if(v>mx)mx=v;}
        int nseg=end-start; double m=sum/(double)nseg;
        if(means)means[s]=m;if(vars)vars[s]=sq/(double)nseg-m*m;
        if(mins)mins[s]=mn;if(maxs)maxs[s]=mx;}
}

/* ── Poisson disc sampling for k-NN validation ──────────── */
void te_ts_poisson_disc_sample(const TETimeSeries *ts, double min_dist, int *indices, int *n_indices, int max_indices) {
    if (!ts || !indices || !n_indices || max_indices < 1) return;
    *n_indices = 0; indices[0] = 0; (*n_indices)++;
    for (int i = 1; i < ts->length && *n_indices < max_indices; i++) {
        int valid = 1;
        for (int j = 0; j < *n_indices; j++) { double d = fabs(ts->data[i] - ts->data[indices[j]]); if (d < min_dist) { valid = 0; break; } }
        if (valid) indices[(*n_indices)++] = i;
    }
}

/* ── BDS test for nonlinearity ──────────────────────────── */
double te_bds_statistic(const TETimeSeries *ts, int embed_dim, double epsilon) {
    if (!ts || embed_dim < 2) return 0.0;
    int n = ts->length - embed_dim + 1; if (n < 50) return 0.0;
    double c_m = 0.0; int count = 0;
    for (int i = 0; i < n; i++) { for (int j = i + 1; j < n; j++) {
        double max_d = 0.0; for (int d = 0; d < embed_dim; d++) { double diff = fabs(ts->data[i+d]-ts->data[j+d]); if (diff > max_d) max_d = diff; }
        if (max_d < epsilon) count++;
    }}
    c_m = 2.0 * (double)count / (double)(n * (n - 1));
    double c_1 = 0.0; count = 0;
    for (int i = 0; i < ts->length; i++) for (int j = i + 1; j < ts->length; j++) if (fabs(ts->data[i]-ts->data[j]) < epsilon) count++;
    c_1 = 2.0 * (double)count / (double)(ts->length * (ts->length - 1));
    double sigma2 = 0.0, k_sum = 0.0;
    for (int k = 1; k < embed_dim; k++) { double pk = pow(c_1, (double)k); k_sum += pow(c_1, (double)(2*(embed_dim-k))); sigma2 += pk*pk; }
    sigma2 = 4.0 * (k_sum + (double)(embed_dim-1)*(double)(embed_dim-1)*pow(c_1,2.0*(double)embed_dim) - (double)embed_dim*(double)embed_dim*pow(c_1,2.0)*pow(c_1,2.0*(double)(embed_dim-1)));
    for (int j = 1; j < embed_dim; j++) sigma2 += 8.0 * (pow(c_1,2.0*(double)j) * (pow(c_1,(double)(embed_dim-j)) - pow(c_1,2.0*(double)(embed_dim-j))));
    if (sigma2 < 1e-15) return 0.0;
    return sqrt((double)n) * (c_m - pow(c_1, (double)embed_dim)) / sqrt(sigma2);
}

/* ── Mann-Whitney U test ────────────────────────────────── */
double te_mann_whitney_u(const double *a, int na, const double *b, int nb, double *p_value) {
    if (!a || !b || na < 2 || nb < 2) return 0.0;
    int *ranks_a = malloc((size_t)na * sizeof(int));
    if (!ranks_a) return 0.0;
    for (int i = 0; i < na; i++) { ranks_a[i] = 1; for (int j = 0; j < na; j++) if (a[j] < a[i]) ranks_a[i]++; for (int j = 0; j < nb; j++) if (b[j] < a[i]) ranks_a[i]++; }
    int R1 = 0; for (int i = 0; i < na; i++) R1 += ranks_a[i];
    double U1 = (double)R1 - (double)na * (double)(na + 1) / 2.0;
    double U2 = (double)na * (double)nb - U1;
    double U = U1 < U2 ? U1 : U2;
    double mu = (double)na * (double)nb / 2.0;
    double sigma = sqrt((double)na * (double)nb * (double)(na + nb + 1) / 12.0);
    double z = (sigma > 1e-15) ? (U - mu) / sigma : 0.0;
    if (p_value) { double abs_z = fabs(z), x = abs_z / 1.4142135623730951;
        double t = 1.0 / (1.0 + 0.3275911 * x);
        double erf = 1.0 - ((((1.061405429*t-1.453152027)*t+1.421413741)*t-0.284496736)*t+0.254829592)*t*exp(-x*x);
        *p_value = 2.0 * (1.0 - erf); }
    free(ranks_a); return z;
}

/* ── Time series summary ────────────────────────────────── */
void te_ts_summary(const TETimeSeries *ts, char *buf, int buf_size) {
    if (!ts || !buf || buf_size < 1) return;
    snprintf(buf, (size_t)buf_size,
        "%s: n=%d mean=%.4f std=%.4f min=%.4f max=%.4f skew=%.4f kurt=%.4f",
        ts->name ? ts->name : "?", ts->length,
        te_ts_mean(ts), te_ts_std(ts), te_ts_min(ts), te_ts_max(ts),
        te_ts_skewness(ts), te_ts_kurtosis(ts));
}
void te_ts_plot_ascii(const TETimeSeries *ts, int width, int height) {
    if (!ts || width < 10 || height < 5) return;
    double mn = te_ts_min(ts), mx = te_ts_max(ts), rng = mx - mn;
    if (rng < TE_EPSILON) rng = 1.0;
    int step = ts->length / width; if (step < 1) step = 1;
    char line[256];
    for (int row = height - 1; row >= 0; row--) {
        double thresh = mn + (double)row * rng / (double)(height - 1);
        int pos = 0; line[pos++] = '|';
        for (int i = 0; i < ts->length && pos < 254; i += step) {
            line[pos++] = (ts->data[i] >= thresh) ? '#' : ' ';
        }
        line[pos] = '\0'; printf("%s\n", line);
    }
    for (int i = 0; i < width + 2; i++) printf("-");
    printf("\n");
}

/* ── Data windowing utilities ───────────────────────────── */
void te_ts_hamming_window(TETimeSeries *ts) {
    if (!ts || ts->length < 2) return;
    for (int i = 0; i < ts->length; i++) ts->data[i] *= 0.54 - 0.46 * cos(6.283185307179586 * (double)i / (double)(ts->length - 1));
}
void te_ts_hann_window(TETimeSeries *ts) {
    if (!ts || ts->length < 2) return;
    for (int i = 0; i < ts->length; i++) ts->data[i] *= 0.5 * (1.0 - cos(6.283185307179586 * (double)i / (double)(ts->length - 1)));
}
void te_ts_blackman_window(TETimeSeries *ts) {
    if (!ts || ts->length < 2) return;
    for (int i = 0; i < ts->length; i++) ts->data[i] *= 0.42 - 0.5*cos(6.283185307179586*(double)i/(double)(ts->length-1)) + 0.08*cos(12.566370614359172*(double)i/(double)(ts->length-1));
}

/* ── Welch power spectral density ───────────────────────── */
void te_ts_welch_psd(const TETimeSeries *ts, int nfft, int overlap, double *psd, double *freqs) {
    if (!ts || !psd || nfft < 4 || overlap < 0 || overlap >= nfft) return;
    int step = nfft - overlap, n_seg = (ts->length - nfft) / step + 1;
    if (n_seg < 1) return;
    double *window = malloc((size_t)nfft * sizeof(double));
    if (!window) return;
    for (int i = 0; i < nfft; i++) window[i] = 0.54 - 0.46 * cos(6.283185307179586 * (double)i / (double)(nfft - 1));
    double *seg = malloc((size_t)nfft * sizeof(double));
    if (!seg) { free(window); return; }
    memset(psd, 0, (size_t)(nfft / 2 + 1) * sizeof(double));
    double win_pow = 0.0; for (int i = 0; i < nfft; i++) win_pow += window[i] * window[i];
    for (int s = 0; s < n_seg; s++) {
        for (int i = 0; i < nfft; i++) seg[i] = ts->data[s * step + i] * window[i];
        for (int k = 0; k <= nfft / 2; k++) { double re = 0.0; for (int j = 0; j < nfft; j++) re += seg[j] * cos(6.283185307179586 * (double)(k * j) / (double)nfft);
            psd[k] += (re * re) / ((double)nfft * win_pow); }
    }
    for (int k = 0; k <= nfft / 2; k++) { psd[k] /= (double)n_seg; if (freqs) freqs[k] = (double)k / (double)nfft; }
    free(window); free(seg);
}
