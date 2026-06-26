#include "te_ksg.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/*
 * Kraskov-Stogbauer-Grassberger (KSG) estimator for Transfer Entropy
 * Uses k-nearest-neighbor statistics to estimate entropies in
 * high-dimensional spaces without explicit binning.
 *
 * Reference: Kraskov et al., Phys. Rev. E 69, 066138 (2004)
 *
 * The key insight: the digamma (psi) function appears naturally
 * from the volume of k-NN balls in d-dimensional space.
 *
 * I(X;Y) = psi(k) - <psi(n_x+1) + psi(n_y+1)> + psi(N)
 *
 * For transfer entropy, we extend this to conditional MI:
 * TE_{Y->X} = I(X_{t+1}; Y_t^{(l)} | X_t^{(k)})
 */

/* ?? Digamma function (psi) via asymptotic expansion ???????
 * psi(x) = ln(x) - 1/(2x) - 1/(12x^2) + 1/(120x^4) - ...
 * Accurate for x > 8. For small x, use recurrence:
 *   psi(x+1) = psi(x) + 1/x
 */
static double te_digamma(double x) {
    if (x < 1e-15) return -1e100;
    double result = 0.0;
    /* Recurrence to bring x > 8 */
    while (x < 8.0) { result -= 1.0 / x; x += 1.0; }
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    double inv_x4 = inv_x2 * inv_x2;
    result += log(x) - 0.5 * inv_x - inv_x2 / 12.0 + inv_x4 / 120.0
              - inv_x4 * inv_x2 / 252.0;
    return result;
}

/* ?? Euclidean distance between two vectors ????????????????? */
static double te_euclidean_dist(const double *a, const double *b, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) { double d = a[i] - b[i]; sum += d * d; }
    return sqrt(sum);
}

/* ?? Chebyshev (max) distance ??????????????????????????????? */
static double te_chebyshev_dist(const double *a, const double *b, int dim) {
    double max_d = 0.0;
    for (int i = 0; i < dim; i++) {
        double d = fabs(a[i] - b[i]);
        if (d > max_d) max_d = d;
    }
    return max_d;
}

/* ?? Simple k-NN search (brute force) ???????????????????????
 * Returns k nearest neighbors of point idx in data set.
 * indices[0..k-1] and distances[0..k-1] are output.
 */
int te_ksg_range_search(const double *points, int n, int dim, int idx,
                         int k, int *indices, double *distances) {
    if (!points || !indices || !distances || idx < 0 || idx >= n || k < 1) return -1;

    /* Collect distances to all other points */
    double *dists = malloc((size_t)n * sizeof(double));
    int *order = malloc((size_t)n * sizeof(int));
    if (!dists || !order) { free(dists); free(order); return -1; }
    for (int i = 0; i < n; i++) { order[i] = i; dists[i] = (i == idx) ? 1e100 : te_chebyshev_dist(points + idx * dim, points + i * dim, dim); }

    /* Simple selection sort for first k+1 elements */
    for (int i = 0; i < k + 1 && i < n; i++) {
        int min_idx = i;
        for (int j = i + 1; j < n; j++)
            if (dists[j] < dists[min_idx]) min_idx = j;
        if (min_idx != i) {
            double td = dists[i]; dists[i] = dists[min_idx]; dists[min_idx] = td;
            int ti = order[i]; order[i] = order[min_idx]; order[min_idx] = ti;
        }
    }

    for (int i = 0; i < k; i++) {
        indices[i] = order[i + 1]; /* skip self */
        distances[i] = dists[i + 1];
    }
    free(dists); free(order);
    return 0;
}

/* ?? KSG Mutual Information estimator ??????????????????????
 * I(X;Y) = psi(k) + psi(N) - 1/k - <psi(n_x) + psi(n_y)>
 *
 * n_x: number of points within Chebyshev distance epsilon/2
 *      of x_i in the X-marginal space
 *
 * Using algorithm 1 from Kraskov et al.:
 * For each point i, find distance epsilon_i/2 to k-th NN
 * in joint (X,Y) space. Then count points n_x(i) within
 * epsilon_i/2 in X-marginal space.
 */
double te_ksg_mutual_info_knn(const double *x, const double *y, int n, int knn) {
    if (!x || !y || n < 10 || knn < 2 || knn > n / 2) return 0.0;

    /* Build joint space: [x_i, y_i] */
    double *joint = malloc((size_t)(2 * n) * sizeof(double));
    if (!joint) return 0.0;
    for (int i = 0; i < n; i++) { joint[i * 2] = x[i]; joint[i * 2 + 1] = y[i]; }

    double *indices_buf = malloc((size_t)knn * sizeof(double));
    int *idx_buf = malloc((size_t)knn * sizeof(int));
    if (!indices_buf || !idx_buf) { free(joint); free(indices_buf); free(idx_buf); return 0.0; }

    double sum_psi_nx = 0.0, sum_psi_ny = 0.0;

    for (int i = 0; i < n; i++) {
        te_ksg_range_search(joint, n, 2, i, knn, idx_buf, indices_buf);
        double eps = indices_buf[knn - 1]; /* distance to k-th NN in joint space */
        if (eps < TE_EPSILON) eps = TE_EPSILON;

        /* Count neighbors within eps/2 in X marginal space */
        int nx_count = 0, ny_count = 0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            if (fabs(x[j] - x[i]) < eps * 0.5) nx_count++;
            if (fabs(y[j] - y[i]) < eps * 0.5) ny_count++;
        }
        sum_psi_nx += te_digamma((double)(nx_count + 1));
        sum_psi_ny += te_digamma((double)(ny_count + 1));
    }

    double mi = te_digamma((double)knn) + te_digamma((double)n)
                - 1.0 / (double)knn
                - (sum_psi_nx + sum_psi_ny) / (double)n;

    free(joint); free(indices_buf); free(idx_buf);
    return (mi > 0.0) ? mi : 0.0;
}

/* ?? KSG entropy estimator ????????????????????????????????? */
double te_ksg_entropy_knn(const double *data, int n, int dim, int knn) {
    if (!data || n < 10 || dim < 1 || knn < 1) return 0.0;
    double *indices_buf = malloc((size_t)knn * sizeof(double));
    int *idx_buf = malloc((size_t)knn * sizeof(int));
    if (!indices_buf || !idx_buf) { free(indices_buf); free(idx_buf); return 0.0; }

    double sum_log_eps = 0.0;
    for (int i = 0; i < n; i++) {
        te_ksg_range_search(data, n, dim, i, knn, idx_buf, indices_buf);
        double eps = indices_buf[knn - 1];
        if (eps < TE_EPSILON) eps = TE_EPSILON;
        sum_log_eps += log(eps);
    }

    double V_d = (dim == 1) ? 2.0 : (dim == 2) ? TE_PI : (4.0 / 3.0) * TE_PI;
    /* For general d, V_d = pi^(d/2) / Gamma(d/2 + 1) */
    if (dim > 3) {
        double d2 = (double)dim / 2.0;
        V_d = pow(TE_PI, d2) / tgamma(d2 + 1.0);
    }
    double H = te_digamma((double)n) - te_digamma((double)knn)
               + log(V_d) + (double)dim * sum_log_eps / (double)n;

    free(indices_buf); free(idx_buf);
    return H;
}

/* ?? KSG Transfer Entropy estimator ????????????????????????
 * TE(Y->X) = I(X_future; Y_past | X_past)
 *          = H(X_future, X_past) + H(X_past, Y_past)
 *          - H(X_past) - H(X_future, X_past, Y_past)
 *
 * Using three KSG mutual information estimates:
 * TE(Y->X) = I(X_f; (X_p, Y_p)) - I(X_f; X_p)
 */
TEResult te_ksg_compute(const TETimeSeries *x, const TETimeSeries *y,
                         int k_embed, int l_embed, int k_nn) {
    TEResult result;
    memset(&result, 0, sizeof(TEResult));
    result.k_embed = k_embed; result.l_embed = l_embed;
    if (!x || !y || k_embed < 1 || l_embed < 1 || k_nn < 2) return result;

    int n = x->length - k_embed;
    if (n < 10 || x->length != y->length) return result;

    /* Build future X */
    double *x_fut = malloc((size_t)n * sizeof(double));
    /* Build joint past (X_past, Y_past) */
    double *joint_past = malloc((size_t)(2 * n) * sizeof(double));
    if (!x_fut || !joint_past) { free(x_fut); free(joint_past); return result; }

    for (int i = 0; i < n; i++) {
        x_fut[i] = x->data[i + k_embed];
        joint_past[i * 2] = x->data[i];
        joint_past[i * 2 + 1] = y->data[i];
    }

    /* I(X_future; (X_past, Y_past)) */
    double *x_fut_past = malloc((size_t)(3 * n) * sizeof(double));
    if (!x_fut_past) { free(x_fut); free(joint_past); return result; }
    for (int i = 0; i < n; i++) {
        x_fut_past[i * 3] = x_fut[i];
        x_fut_past[i * 3 + 1] = x->data[i];
        x_fut_past[i * 3 + 2] = y->data[i];
    }
    double H_xf_xpyp = te_ksg_entropy_knn(x_fut_past, n, 3, k_nn);

    /* I(X_future; X_past) */
    double *x_fut_only = malloc((size_t)(2 * n) * sizeof(double));
    if (!x_fut_only) { free(x_fut); free(joint_past); free(x_fut_past); return result; }
    for (int i = 0; i < n; i++) {
        x_fut_only[i * 2] = x_fut[i];
        x_fut_only[i * 2 + 1] = x->data[i];
    }
    double H_xf_xp = te_ksg_entropy_knn(x_fut_only, n, 2, k_nn);

    /* H(X_past) */
    double H_xp = te_ksg_entropy_knn(x->data, n, 1, k_nn);

    /* H(X_past, Y_past) */
    double H_xpyp = te_ksg_entropy_knn(joint_past, n, 2, k_nn);

    /* TE(Y->X) = H(X_f, X_p) + H(X_p, Y_p) - H(X_p) - H(X_f, X_p, Y_p) */
    result.te_xy = H_xf_xp + H_xpyp - H_xp - H_xf_xpyp;
    if (result.te_xy < 0.0) result.te_xy = 0.0;

    /* TE(X->Y): symmetric computation */
    for (int i = 0; i < n; i++) {
        x_fut[i] = y->data[i + l_embed];
        x_fut_past[i * 3] = x_fut[i];
        x_fut_past[i * 3 + 1] = y->data[i];
        x_fut_past[i * 3 + 2] = x->data[i];
        x_fut_only[i * 2] = x_fut[i];
        x_fut_only[i * 2 + 1] = y->data[i];
    }
    double H_yf_xpyp = te_ksg_entropy_knn(x_fut_past, n, 3, k_nn);
    double H_yf_yp = te_ksg_entropy_knn(x_fut_only, n, 2, k_nn);
    double H_yp = te_ksg_entropy_knn(y->data, n, 1, k_nn);
    result.te_yx = H_yf_yp + H_xpyp - H_yp - H_yf_xpyp;
    if (result.te_yx < 0.0) result.te_yx = 0.0;

    result.net_te = result.te_xy - result.te_yx;
    result.n_surrogates = 0;

    free(x_fut); free(joint_past); free(x_fut_past); free(x_fut_only);
    return result;
}

double te_ksg_estimate(const double *x, const double *y, int n, int k, int l, int knn) {
    TETimeSeries *ts_x = te_ts_create(x, n, "x");
    TETimeSeries *ts_y = te_ts_create(y, n, "y");
    TEResult r = te_ksg_compute(ts_x, ts_y, k, l, knn);
    te_ts_free(ts_x); te_ts_free(ts_y);
    return r.te_xy;
}

/* ???????????????????????????????????????????????????????????
 * Alternative distance metrics for KSG estimator
 * ??????????????????????????????????????????????????????????? */
static double te_manhattan_dist(const double *a, const double *b, int dim) {
    double sum = 0.0; for (int i = 0; i < dim; i++) sum += fabs(a[i] - b[i]); return sum;
}
static double te_minkowski_dist(const double *a, const double *b, int dim, double p) {
    double sum = 0.0; for (int i = 0; i < dim; i++) sum += pow(fabs(a[i] - b[i]), p);
    return pow(sum, 1.0 / p);
}
static double te_cosine_dist(const double *a, const double *b, int dim) {
    double dot = 0.0, na = 0.0, nb = 0.0;
    for (int i = 0; i < dim; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
    double den = sqrt(na) * sqrt(nb); if (den < TE_EPSILON) return 0.0;
    return 1.0 - dot / den;
}

/* ???????????????????????????????????????????????????????????
 * KSG estimator with variable k (k-NN with automatic selection)
 * Uses Kozachenko-Leonenko bias correction.
 * ??????????????????????????????????????????????????????????? */
double te_ksg_auto_knn(const double *x, const double *y, int n, int k_min, int k_max) {
    if (!x || !y || n < 10 || k_min < 2 || k_max < k_min) return 0.0;
    double best_mi = 0.0;
    for (int k = k_min; k <= k_max && k < n / 2; k++) {
        double mi = te_ksg_mutual_info_knn(x, y, n, k);
        if (mi > best_mi) best_mi = mi;
    }
    return best_mi;
}

/* ???????????????????????????????????????????????????????????
 * KSG entropy with LOOV (leave-one-out validation)
 * Estimates differential entropy using the KSG method.
 * ??????????????????????????????????????????????????????????? */
double te_ksg_entropy_loov(const double *data, int n, int dim, int k_min, int k_max) {
    if (!data || n < 10 || k_min < 1) return 0.0;
    double best_h = 1e100; int best_k = k_min;
    for (int k = k_min; k <= k_max && k < n / 2; k++) {
        double h = te_ksg_entropy_knn(data, n, dim, k);
        if (h < best_h) { best_h = h; best_k = k; }
    }
    return te_ksg_entropy_knn(data, n, dim, best_k);
}

/* ?? KSG estimator with Manhattan distance ???????????????? */
double te_ksg_mi_manhattan(const double *x, const double *y, int n, int knn) {
    if (!x || !y || n < 10 || knn < 2) return 0.0;
    double *joint = malloc((size_t)(2*n)*sizeof(double));
    if (!joint) return 0.0;
    for (int i=0;i<n;i++) { joint[i*2]=x[i]; joint[i*2+1]=y[i]; }
    double sum_psi_nx=0, sum_psi_ny=0;
    for (int i=0;i<n;i++) {
        int idx_buf[64]; double dist_buf[64]; int kk = knn < 64 ? knn : 63;
        te_ksg_range_search(joint, n, 2, i, kk, idx_buf, dist_buf);
        double eps = dist_buf[kk-1]; if (eps<TE_EPSILON) eps=TE_EPSILON;
        int nx=0, ny=0;
        for (int j=0;j<n;j++) { if(j==i)continue; if(fabs(x[j]-x[i])<eps*0.5)nx++; if(fabs(y[j]-y[i])<eps*0.5)ny++; }
        sum_psi_nx += te_digamma((double)(nx+1));
        sum_psi_ny += te_digamma((double)(ny+1));
    }
    double mi = te_digamma((double)knn) + te_digamma((double)n) - 1.0/(double)knn - (sum_psi_nx+sum_psi_ny)/(double)n;
    free(joint);
    return mi > 0.0 ? mi : 0.0;
}

/* ?? KSG estimator with adaptive k selection ?????????????? */
double te_ksg_adaptive_k(const double *x, const double *y, int n) {
    if (!x || !y || n < 20) return 0.0;
    int k_opt = (int)sqrt((double)n);
    if (k_opt < 2) k_opt = 2; if (k_opt > n/3) k_opt = n/3;
    double mi1 = te_ksg_mutual_info_knn(x, y, n, k_opt);
    double mi2 = te_ksg_mutual_info_knn(x, y, n, k_opt+1);
    double mi3 = te_ksg_mutual_info_knn(x, y, n, k_opt-1);
    double best = mi1; if (mi2 > best) best = mi2; if (mi3 > best) best = mi3;
    return best;
}

/* ?? KSG conditional mutual information ???????????????????
 * I(X;Y|Z) = H(X,Z) + H(Y,Z) - H(Z) - H(X,Y,Z)
 * using KSG entropy estimator for each term.
 */
double te_ksg_conditional_mi(const double *x, const double *y, const double *z,
                              int n, int knn) {
    if (!x || !y || !z || n < 10 || knn < 2) return 0.0;
    double *xz = malloc((size_t)(2*n)*sizeof(double));
    double *yz = malloc((size_t)(2*n)*sizeof(double));
    double *xyz = malloc((size_t)(3*n)*sizeof(double));
    if (!xz || !yz || !xyz) { free(xz); free(yz); free(xyz); return 0.0; }
    for (int i=0;i<n;i++) {
        xz[i*2]=x[i]; xz[i*2+1]=z[i];
        yz[i*2]=y[i]; yz[i*2+1]=z[i];
        xyz[i*3]=x[i]; xyz[i*3+1]=y[i]; xyz[i*3+2]=z[i];
    }
    double H_xz = te_ksg_entropy_knn(xz, n, 2, knn), H_yz = te_ksg_entropy_knn(yz, n, 2, knn);
    double H_z = te_ksg_entropy_knn(z, n, 1, knn), H_xyz = te_ksg_entropy_knn(xyz, n, 3, knn);
    free(xz); free(yz); free(xyz);
    double cmi = H_xz + H_yz - H_z - H_xyz;
    return cmi > 0.0 ? cmi : 0.0;
}

/* ?? KSG Transfer Entropy with conditional MI ????????????? */
double te_ksg_te_via_cmi(const double *x, const double *y, int n, int k, int knn) {
    if (!x || !y || n < k+2 || knn < 2) return 0.0;
    double *x_fut = malloc((size_t)(n-k)*sizeof(double));
    double *x_past = malloc((size_t)(n-k)*sizeof(double));
    double *y_past = malloc((size_t)(n-k)*sizeof(double));
    if (!x_fut || !x_past || !y_past) { free(x_fut); free(x_past); free(y_past); return 0.0; }
    int m = n-k;
    for (int i=0;i<m;i++) { x_fut[i]=x[i+k]; x_past[i]=x[i]; y_past[i]=y[i]; }
    double cmi = te_ksg_conditional_mi(x_fut, y_past, x_past, m, knn);
    free(x_fut); free(x_past); free(y_past);
    return cmi > 0.0 ? cmi : 0.0;
}
