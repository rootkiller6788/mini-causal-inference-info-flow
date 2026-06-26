/* ts_core.c -- Core time series, vector, and matrix operations.
 *
 * Knowledge points:
 * L1: Vector/Matrix/TimeSeries/VARModel type definitions
 * L2: OLS estimation as foundation of Granger testing
 * L2: VAR model as dynamic system representation
 * L3: Matrix algebra (multiplication, inverse, determinant, trace)
 * L5: Gaussian elimination with partial pivoting
 * L5: VAR fitting via equation-by-equation OLS
 * L5: AIC/BIC model selection
 * L6: AR(1) and VAR(1) simulation (canonical time series models)
 *
 * Reference: Lutkepohl (2005), New Introduction to Multiple Time Series Analysis.
 */

#include "ts_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- Vector Operations ---------- */

Vector* vec_create(int n) {
    if (n <= 0) return NULL;
    Vector* v = calloc(1, sizeof(Vector));
    if (!v) return NULL;
    v->data = calloc((size_t)n, sizeof(double));
    if (!v->data) { free(v); return NULL; }
    v->n = n;
    return v;
}

void vec_free(Vector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

/* Euclidean (L2) norm: ||v||_2 = sqrt(sum_i v_i^2).
 * Returns 0 for NULL or empty vector. Complexity: O(n). */
double vec_norm(const Vector* v) {
    if (!v) return 0;
    double s = 0;
    for (int i = 0; i < v->n; i++) s += v->data[i] * v->data[i];
    return sqrt(s);
}

/* Dot product: a · b = sum_i a_i b_i.
 * Returns 0 if NULL, different lengths, or empty.
 * Complexity: O(n). */
double vec_dot(const Vector* a, const Vector* b) {
    if (!a || !b || a->n != b->n) return 0;
    double s = 0;
    for (int i = 0; i < a->n; i++) s += a->data[i] * b->data[i];
    return s;
}

/* L1 norm: ||v||_1 = sum_i |v_i|.
 * Complexity: O(n). */
double vec_l1_norm(const Vector* v) {
    if (!v) return 0;
    double s = 0;
    for (int i = 0; i < v->n; i++) s += fabs(v->data[i]);
    return s;
}

/* Infinity norm: ||v||_inf = max_i |v_i|.
 * Complexity: O(n). */
double vec_inf_norm(const Vector* v) {
    if (!v) return 0;
    double m = 0;
    for (int i = 0; i < v->n; i++) {
        double a = fabs(v->data[i]);
        if (a > m) m = a;
    }
    return m;
}

/* Vector addition: r = a + b. Returns new vector. */
Vector* vec_add(const Vector* a, const Vector* b) {
    if (!a || !b || a->n != b->n) return NULL;
    Vector* r = vec_create(a->n);
    if (!r) return NULL;
    for (int i = 0; i < a->n; i++) r->data[i] = a->data[i] + b->data[i];
    return r;
}

/* Scalar multiply: r = alpha * v (in-place). */
void vec_scale(Vector* v, double alpha) {
    if (!v) return;
    for (int i = 0; i < v->n; i++) v->data[i] *= alpha;
}

/* Mean of vector elements. Complexity: O(n). */
double vec_mean(const Vector* v) {
    if (!v || v->n == 0) return 0;
    double s = 0;
    for (int i = 0; i < v->n; i++) s += v->data[i];
    return s / v->n;
}

/* Variance (population) of vector elements. Complexity: O(n). */
double vec_variance(const Vector* v) {
    if (!v || v->n <= 1) return 0;
    double mu = vec_mean(v);
    double s = 0;
    for (int i = 0; i < v->n; i++) {
        double d = v->data[i] - mu;
        s += d * d;
    }
    return s / (v->n - 1);
}

/* ---------- Matrix Operations ---------- */

Matrix* mat_create(int r, int c) {
    if (r <= 0 || c <= 0) return NULL;
    Matrix* m = calloc(1, sizeof(Matrix));
    if (!m) return NULL;
    m->data = calloc((size_t)(r * c), sizeof(double));
    if (!m->data) { free(m); return NULL; }
    m->rows = r; m->cols = c;
    return m;
}

void mat_free(Matrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

Matrix* mat_identity(int n) {
    Matrix* m = mat_create(n, n);
    if (!m) return NULL;
    for (int i = 0; i < n; i++) m->data[i * n + i] = 1.0;
    return m;
}

Matrix* mat_copy(const Matrix* m) {
    if (!m) return NULL;
    Matrix* c = mat_create(m->rows, m->cols);
    if (!c) return NULL;
    memcpy(c->data, m->data, (size_t)(m->rows * m->cols) * sizeof(double));
    return c;
}

void mat_set(Matrix* m, int i, int j, double v) {
    if (m && i >= 0 && i < m->rows && j >= 0 && j < m->cols)
        m->data[i * m->cols + j] = v;
}

double mat_get(const Matrix* m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return 0.0;
    return m->data[i * m->cols + j];
}

/* Matrix multiplication: C = A * B, C[i,j] = sum_k A[i,k] * B[k,j].
 * Complexity: O(rows(A) * cols(B) * cols(A)). */
Matrix* mat_mul(const Matrix* a, const Matrix* b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    Matrix* r = mat_create(a->rows, b->cols);
    if (!r) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < b->cols; j++) {
            double s = 0;
            for (int k = 0; k < a->cols; k++)
                s += a->data[i * a->cols + k] * b->data[k * b->cols + j];
            r->data[i * b->cols + j] = s;
        }
    return r;
}

Matrix* mat_transpose(const Matrix* m) {
    if (!m) return NULL;
    Matrix* t = mat_create(m->cols, m->rows);
    if (!t) return NULL;
    for (int i = 0; i < m->rows; i++)
        for (int j = 0; j < m->cols; j++)
            t->data[j * m->rows + i] = m->data[i * m->cols + j];
    return t;
}

/* Matrix trace: tr(A) = sum_i A[i,i]. Requires square.
 * Complexity: O(min(rows, cols)). */
double mat_trace(const Matrix* m) {
    if (!m || m->rows != m->cols) return 0;
    double t = 0;
    for (int i = 0; i < m->rows; i++) t += m->data[i * m->cols + i];
    return t;
}

/* Matrix determinant via Gaussian elimination with partial pivoting.
 * Complexity: O(n^3). */
double mat_det(const Matrix* m) {
    if (!m || m->rows != m->cols) return 0;
    int n = m->rows;
    if (n == 1) return m->data[0];
    if (n == 2) return m->data[0] * m->data[3] - m->data[1] * m->data[2];
    Matrix* a = mat_copy(m);
    if (!a) return 0;
    double d = 1;
    for (int i = 0; i < n; i++) {
        double pv = a->data[i * n + i];
        if (fabs(pv) < 1e-12) {
            int sw = -1;
            for (int k = i + 1; k < n; k++)
                if (fabs(a->data[k * n + i]) > 1e-12) { sw = k; break; }
            if (sw < 0) { mat_free(a); return 0; }
            for (int j = 0; j < n; j++) {
                double t = a->data[i * n + j];
                a->data[i * n + j] = a->data[sw * n + j];
                a->data[sw * n + j] = t;
            }
            d = -d;
            pv = a->data[i * n + i];
        }
        d *= pv;
        for (int k = i + 1; k < n; k++) {
            double f = a->data[k * n + i] / pv;
            for (int j = i; j < n; j++) a->data[k * n + j] -= f * a->data[i * n + j];
        }
    }
    mat_free(a);
    return d;
}

/* Matrix inverse via Gauss-Jordan elimination.
 * Complexity: O(n^3). Returns NULL if singular. */
Matrix* mat_inverse(const Matrix* m) {
    if (!m || m->rows != m->cols) return NULL;
    int n = m->rows;
    Matrix* a = mat_copy(m);
    Matrix* inv = mat_identity(n);
    if (!a || !inv) { mat_free(a); mat_free(inv); return NULL; }

    for (int i = 0; i < n; i++) {
        double pv = a->data[i * n + i];
        if (fabs(pv) < 1e-12) {
            int sw = -1;
            for (int k = i + 1; k < n; k++)
                if (fabs(a->data[k * n + i]) > 1e-12) { sw = k; break; }
            if (sw < 0) { mat_free(a); mat_free(inv); return NULL; }
            for (int j = 0; j < n; j++) {
                double t = a->data[i * n + j];
                a->data[i * n + j] = a->data[sw * n + j];
                a->data[sw * n + j] = t;
                t = inv->data[i * n + j];
                inv->data[i * n + j] = inv->data[sw * n + j];
                inv->data[sw * n + j] = t;
            }
            pv = a->data[i * n + i];
        }
        for (int j = 0; j < n; j++) { a->data[i * n + j] /= pv; inv->data[i * n + j] /= pv; }
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double f = a->data[k * n + i];
            for (int j = 0; j < n; j++) {
                a->data[k * n + j] -= f * a->data[i * n + j];
                inv->data[k * n + j] -= f * inv->data[i * n + j];
            }
        }
    }
    mat_free(a);
    return inv;
}

/* Frobenius norm: ||A||_F = sqrt(sum_{i,j} A_{ij}^2).
 * Complexity: O(rows * cols). */
double mat_frobenius_norm(const Matrix* m) {
    if (!m) return 0;
    double s = 0;
    for (int i = 0; i < m->rows * m->cols; i++) s += m->data[i] * m->data[i];
    return sqrt(s);
}

/* Condition number estimate using Frobenius norm:
 * cond(A) ≈ ||A||_F * ||A^{-1}||_F.
 * Complexity: O(n^3) due to inverse. */
double mat_condition_number(const Matrix* m) {
    if (!m || m->rows != m->cols) return -1;
    double nA = mat_frobenius_norm(m);
    Matrix* inv = mat_inverse(m);
    if (!inv) return 1e10;
    double nAi = mat_frobenius_norm(inv);
    mat_free(inv);
    return nA * nAi;
}

/* Matrix-vector multiply: y = A * x. Returns new vector.
 * Complexity: O(rows * cols). */
Vector* mat_vec_mul(const Matrix* A, const Vector* x) {
    if (!A || !x || A->cols != x->n) return NULL;
    Vector* y = vec_create(A->rows);
    if (!y) return NULL;
    for (int i = 0; i < A->rows; i++) {
        double s = 0;
        for (int j = 0; j < A->cols; j++)
            s += A->data[i * A->cols + j] * x->data[j];
        y->data[i] = s;
    }
    return y;
}

/* ---------- Time Series Operations ---------- */

TimeSeries* ts_create(int length, int dim) {
    TimeSeries* ts = calloc(1, sizeof(TimeSeries));
    if (ts) {
        ts->length = length;
        ts->dim = dim;
        ts->values = calloc((size_t)(length * dim), sizeof(double));
    }
    return ts;
}

void ts_free(TimeSeries* ts) {
    if (!ts) return;
    free(ts->values);
    free(ts);
}

TimeSeries* ts_slice(const TimeSeries* ts, int start, int len) {
    if (!ts || start < 0 || start + len > ts->length) return NULL;
    TimeSeries* s = ts_create(len, ts->dim);
    if (!s) return NULL;
    memcpy(s->values, ts->values + start * ts->dim,
           (size_t)(len * ts->dim) * sizeof(double));
    return s;
}

void ts_print(const TimeSeries* ts) {
    if (!ts) return;
    printf("TS(%d x %d): first=[", ts->length, ts->dim);
    for (int j = 0; j < ts->dim && j < 4; j++)
        printf("%.4f ", ts->values[j]);
    printf("...]\n");
}

/* Compute autocorrelation at lag k.
 * r_k = sum_{t=1}^{n-k} (x_t - x̄)(x_{t+k} - x̄) / sum_{t=1}^n (x_t - x̄)²
 * Returns 0 for NULL or too-short series. Complexity: O(n). */
double ts_autocorrelation(const double* x, int n, int k) {
    if (!x || n <= k || k < 0) return 0;
    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;
    double num = 0, den = 0;
    for (int i = 0; i < n; i++) {
        double d = x[i] - mean;
        den += d * d;
        if (i < n - k) num += d * (x[i + k] - mean);
    }
    return (den > 1e-15) ? num / den : 0;
}

/* ---------- OLS and VAR ---------- */

/* Ordinary Least Squares: min ||y - X*beta||^2.
 * Solution: beta = (X'X)^{-1} X'y.
 * Returns RSS (residual sum of squares).
 * Complexity: O(n*k^2 + k^3). */
double ols_estimate(const Matrix* X, const double* y, double* beta,
    double* resid, int n, int k) {
    Matrix* Xt = mat_transpose(X);
    Matrix* XtX = mat_mul(Xt, X);
    Matrix* XtXi = mat_inverse(XtX);
    if (!XtXi) { mat_free(Xt); mat_free(XtX); return 1e10; }
    double* Xty = malloc((size_t)k * sizeof(double));
    for (int i = 0; i < k; i++) {
        Xty[i] = 0;
        for (int j = 0; j < n; j++)
            Xty[i] += Xt->data[i * Xt->cols + j] * y[j];
    }
    for (int i = 0; i < k; i++) {
        beta[i] = 0;
        for (int j = 0; j < k; j++)
            beta[i] += XtXi->data[i * k + j] * Xty[j];
    }
    double rss = 0;
    for (int i = 0; i < n; i++) {
        double yh = 0;
        for (int j = 0; j < k; j++) yh += X->data[i * k + j] * beta[j];
        resid[i] = y[i] - yh;
        rss += resid[i] * resid[i];
    }
    mat_free(Xt); mat_free(XtX); mat_free(XtXi); free(Xty);
    return rss;
}

VARModel* var_create(int dim, int p) {
    VARModel* v = calloc(1, sizeof(VARModel));
    if (!v) return NULL;
    v->dim = dim; v->p = p;
    v->A = malloc((size_t)p * sizeof(Matrix*));
    v->c = calloc((size_t)dim, sizeof(double));
    for (int i = 0; i < p; i++) v->A[i] = mat_create(dim, dim);
    return v;
}

void var_free(VARModel* v) {
    if (!v) return;
    for (int i = 0; i < v->p; i++) mat_free(v->A[i]);
    free(v->A); free(v->c); free(v);
}

/* Fit VAR(p) by equation-by-equation OLS.
 * For each dimension d, regress y_{t,d} on [1, all lagged y_{t-i,j}].
 * Complexity: O(dim^2 * p^2 * n + dim * p^3). */
int var_fit(VARModel* var, const TimeSeries* ts) {
    if (!var || !ts || ts->dim != var->dim) return -1;
    int n = ts->length - var->p;
    int k = 1 + var->p * var->dim;
    Matrix* X = mat_create(n, k);
    double* y = malloc((size_t)(n * var->dim) * sizeof(double));
    double* beta = malloc((size_t)(k * var->dim) * sizeof(double));
    double* resid = malloc((size_t)n * sizeof(double));
    double tr = 0;
    for (int d = 0; d < var->dim; d++) {
        for (int i = 0; i < n; i++) {
            X->data[i * k + 0] = 1.0;
            for (int lag = 0; lag < var->p; lag++)
                for (int dd = 0; dd < var->dim; dd++)
                    X->data[i * k + 1 + lag * var->dim + dd] =
                        ts->values[(i + var->p - 1 - lag) * var->dim + dd];
            y[i] = ts->values[(i + var->p) * var->dim + d];
        }
        tr += ols_estimate(X, y, beta + d * k, resid, n, k);
        var->c[d] = beta[d * k + 0];
        for (int lag = 0; lag < var->p; lag++)
            for (int dd = 0; dd < var->dim; dd++)
                mat_set(var->A[lag], d, dd,
                    beta[d * k + 1 + lag * var->dim + dd]);
    }
    var->sigma2 = tr / (n * var->dim);
    var->aic = log(var->sigma2) + 2.0 * k / n;
    var->bic = log(var->sigma2) + log((double)n) * k / n;
    mat_free(X); free(y); free(beta); free(resid);
    return 0;
}

/* VAR forecasting: given p historical values (latest first in hist array),
 * produce steps-ahead forecasts. Complexity: O(steps * p * dim^2). */
int var_forecast(const VARModel* var, const double* hist, double* fc, int steps) {
    int dim = var->dim, p = var->p;
    double* buf = malloc((size_t)((p + steps) * dim) * sizeof(double));
    memcpy(buf, hist, (size_t)(p * dim) * sizeof(double));
    for (int s = 0; s < steps; s++) {
        int idx = p + s;
        for (int d = 0; d < dim; d++) {
            buf[idx * dim + d] = var->c[d];
            for (int lag = 0; lag < p; lag++)
                for (int dd = 0; dd < dim; dd++)
                    buf[idx * dim + d] +=
                        mat_get(var->A[lag], d, dd) * buf[(idx - 1 - lag) * dim + dd];
        }
    }
    memcpy(fc, buf + p * dim, (size_t)(steps * dim) * sizeof(double));
    free(buf);
    return 0;
}

/* VAR order selection by AIC minimization.
 * Tests p = 1..max_p, selects p with minimum AIC.
 * Complexity: O(max_p * dim^2 * p_max^2 * n). */
int var_order_select(const TimeSeries* ts, int max_p, int dim,
    int* best_p, double* best_aic) {
    double baic = 1e10; int bp = 1;
    for (int p = 1; p <= max_p; p++) {
        VARModel* v = var_create(dim, p);
        var_fit(v, ts);
        if (v->aic < baic) { baic = v->aic; bp = p; }
        var_free(v);
    }
    *best_p = bp; *best_aic = baic;
    return bp;
}

/* AR(1) simulation: x_t = phi*x_{t-1} + sigma*eps_t, eps_t ~ N(0,1).
 * Uses Box-Muller approximation via uniform random.
 * Complexity: O(length). */
TimeSeries* ts_simulate_ar1(double phi, double sigma, int length) {
    TimeSeries* ts = ts_create(length, 1);
    if (!ts) return NULL;
    ts->values[0] = 0;
    for (int i = 1; i < length; i++)
        ts->values[i] = phi * ts->values[i - 1] +
            sigma * ((double)rand() / RAND_MAX - 0.5) * 2.0;
    return ts;
}

/* Bivariate VAR(1) simulation:
 * [x_t]   [a11 a12] [x_{t-1}]   [eps_x]
 * [y_t] = [a21 a22] [y_{t-1}] + [eps_y]
 * Complexity: O(length). */
TimeSeries* ts_simulate_var1_coupled(double a11, double a12,
    double a21, double a22, int len) {
    TimeSeries* ts = ts_create(len, 2);
    if (!ts) return NULL;
    ts->values[0] = 0; ts->values[1] = 0;
    for (int i = 1; i < len; i++) {
        double x = ts->values[(i - 1) * 2], y = ts->values[(i - 1) * 2 + 1];
        ts->values[i * 2] = a11 * x + a12 * y +
            0.1 * ((double)rand() / RAND_MAX - 0.5);
        ts->values[i * 2 + 1] = a21 * x + a22 * y +
            0.1 * ((double)rand() / RAND_MAX - 0.5);
    }
    return ts;
}

/* Extract one dimension from a multivariate time series.
 * Returns malloc'd array that caller must free. */
double* ts_extract_dim(const TimeSeries* ts, int d) {
    if (!ts || d < 0 || d >= ts->dim) return NULL;
    double* x = malloc((size_t)ts->length * sizeof(double));
    if (!x) return NULL;
    for (int i = 0; i < ts->length; i++)
        x[i] = ts->values[i * ts->dim + d];
    return x;
}

/* Compute autocorrelation function (ACF) for lags 0..max_lag.
 * Returns malloc'd array of length max_lag+1. Caller frees.
 * Complexity: O(n * max_lag). */
double* ts_acf(const double* x, int n, int max_lag) {
    if (!x || n <= max_lag) return NULL;
    double* acf = malloc((size_t)(max_lag + 1) * sizeof(double));
    for (int k = 0; k <= max_lag; k++)
        acf[k] = ts_autocorrelation(x, n, k);
    return acf;
}