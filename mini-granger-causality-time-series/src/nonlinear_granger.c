/* nonlinear_granger.c -- Nonlinear Granger Causality via Kernel Methods.
 *
 * Tests whether X nonlinearly Granger-causes Y using kernel ridge
 * regression in reproducing kernel Hilbert space (RKHS).
 *
 * The fundamental idea (Ancona et al. 2004):
 * 1. Map lagged data into RKHS via Gaussian RBF kernel
 * 2. Fit kernel ridge regression: y_future ~ K(history_Y, history_Y)
 *    (restricted) vs y_future ~ K(history_XY, history_XY) (unrestricted)
 * 3. Compare prediction errors to assess X's contribution
 *
 * Knowledge points:
 * L2: Nonlinear Granger causality concept
 * L3: Reproducing kernel Hilbert space (RKHS) structure
 * L5: Kernel ridge regression algorithm
 * L5: Gaussian RBF kernel computation and median heuristic
 * L8: Nonlinear extension of linear Granger causality
 *
 * Reference: Ancona, Marinazzo, Stramaglia (2004), PRE 70(5):056221.
 *            Marinazzo, Pellicoro, Stramaglia (2008), PRE 77(5):056215.
 */

#include "nonlinear_granger.h"
#include "granger_test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Compute squared Euclidean distance between rows i and j of X.
 * Each row is a d-dimensional embedding vector. */
static double sq_dist(const Matrix* X, int i, int j) {
    double d2 = 0.0;
    for (int k = 0; k < X->cols; k++) {
        double diff = X->data[i * X->cols + k] - X->data[j * X->cols + k];
        d2 += diff * diff;
    }
    return d2;
}

/* Compute Gaussian RBF kernel matrix:
 * K_ij = exp(-||x_i - x_j||^2 / (2 * sigma^2))
 * X is n x d, K is n x n.
 * Complexity: O(n^2 * d). */
Matrix* kernel_rbf(const Matrix* X, double sigma) {
    if (!X || sigma <= 0.0) return NULL;
    int n = X->rows;
    Matrix* K = mat_create(n, n);
    if (!K) return NULL;
    double gamma = 1.0 / (2.0 * sigma * sigma);
    for (int i = 0; i < n; i++) {
        K->data[i * n + i] = 1.0;  /* exp(0) = 1 */
        for (int j = i + 1; j < n; j++) {
            double val = exp(-sq_dist(X, i, j) * gamma);
            K->data[i * n + j] = val;
            K->data[j * n + i] = val;
        }
    }
    return K;
}

/* Median heuristic for RBF kernel bandwidth:
 * sigma = median{||x_i - x_j|| : i < j} / sqrt(2)
 * This is the default bandwidth selector for kernel methods. */
double kernel_median_heuristic(const double* x, int n) {
    if (!x || n < 2) return 1.0;
    /* Sample pairs for efficiency (use O(n) pairs instead of O(n^2)) */
    int n_pairs = n < 100 ? n * (n - 1) / 2 : 500;
    double* dists = malloc((size_t)n_pairs * sizeof(double));
    if (!dists) return 1.0;
    int idx = 0;
    for (int i = 0; i < n && idx < n_pairs; i++) {
        int step = (n - i - 1) > 0 ? ((n - i - 1) / ((n_pairs - idx) > 0 ? (n_pairs - idx) : 1) + 1) : 1;
        for (int j = i + 1; j < n && idx < n_pairs; j += step) {
            dists[idx++] = fabs(x[i] - x[j]);
        }
    }
    /* Sort and take median */
    for (int i = 0; i < idx - 1; i++)
        for (int j = i + 1; j < idx; j++)
            if (dists[i] > dists[j]) {
                double tmp = dists[i]; dists[i] = dists[j]; dists[j] = tmp;
            }
    double median = idx > 0 ? dists[idx / 2] : 1.0;
    free(dists);
    return median / sqrt(2.0);
}

/* Solve (K + lambda*I) * alpha = y for alpha.
 * Uses Gaussian elimination since K+lambda*I is symmetric positive definite.
 * Returns alpha vector of length n, caller must free. */
double* kernel_ridge_regression(const Matrix* K, const double* y, int n, double lambda) {
    if (!K || !y || n <= 0) return NULL;
    Matrix* A = mat_create(n, n);
    if (!A) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A->data[i * n + j] = K->data[i * n + j];
            if (i == j) A->data[i * n + j] += lambda;
        }
    }
    Matrix* A_inv = mat_inverse(A);
    if (!A_inv) { mat_free(A); return NULL; }
    double* alpha = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        alpha[i] = 0.0;
        for (int j = 0; j < n; j++)
            alpha[i] += A_inv->data[i * n + j] * y[j];
    }
    mat_free(A); mat_free(A_inv);
    return alpha;
}

/* Build embedding matrix from scalar time series.
 * For lag p, produce (n-p) x p matrix where row t contains
 * [x_{t}, x_{t+1}, ..., x_{t+p-1}]. */
static Matrix* build_embedding(const double* x, int n, int p) {
    int rows = n - p;
    if (rows <= 0 || p <= 0) return NULL;
    Matrix* M = mat_create(rows, p);
    if (!M) return NULL;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < p; j++)
            M->data[i * p + j] = x[i + j];
    return M;
}

/* Compute kernel regression prediction error (MSE).
 * Given training kernel K_train, training targets y_train,
 * and test kernel K_test (cross-kernel between test and train),
 * predict y_test = K_test * alpha and compute MSE. */
static double kernel_prediction_error(const Matrix* K_train, const double* y_train,
    const Matrix* K_test, const double* y_test, int n_train, int n_test, double lambda) {
    double* alpha = kernel_ridge_regression(K_train, y_train, n_train, lambda);
    if (!alpha) return 1e10;
    double mse = 0.0;
    for (int i = 0; i < n_test; i++) {
        double pred = 0.0;
        for (int j = 0; j < n_train; j++)
            pred += K_test->data[i * n_train + j] * alpha[j];
        double err = y_test[i] - pred;
        mse += err * err;
    }
    free(alpha);
    return mse / n_test;
}

NonlinearGrangerResult* nlg_create(void) {
    return calloc(1, sizeof(NonlinearGrangerResult));
}

void nlg_free(NonlinearGrangerResult* r) { free(r); }

/* Nonlinear Granger causality test.
 *
 * Algorithm:
 * 1. Build embeddings for Y history only (restricted) and XY history (unrestricted)
 * 2. Split into training and testing sets (70/30)
 * 3. Compute RBF kernel matrices for both models
 * 4. Fit kernel ridge regression and compute MSE on test set
 * 5. Compute test statistic: T = (MSE_r - MSE_ur) / MSE_ur
 * 6. Use F-like approximation for p-value
 *
 * Complexity: O(n^3) for kernel matrix inversion. */
int nlg_test(NonlinearGrangerResult* r, const double* x, const double* y,
    int length, int max_lag, double sigma, double alpha) {
    if (!r || !x || !y || length < max_lag + 10) return -1;

    /* Select best lag via AIC-like criterion on simple linear model */
    int lag = 1;
    double best_aic = 1e10;
    for (int p = 1; p <= max_lag && p < length / 4; p++) {
        int n_eff = length - p;
        Matrix* M_ur = build_embedding(y, n_eff, 2 * p);
        if (!M_ur) continue;
        double aic_val = 0.0;
        for (int i = 0; i < n_eff - 1 && i < 50; i++) {
            double diff = y[i + p] - y[i + p - 1];
            aic_val += diff * diff;
        }
        aic_val = log(aic_val / (n_eff > 0 ? n_eff : 1)) + 2.0 * p / (n_eff > 0 ? n_eff : 1);
        mat_free(M_ur);
        if (aic_val < best_aic) { best_aic = aic_val; lag = p; }
    }

    int n = length - lag;
    int n_train = n * 7 / 10;
    int n_test = n - n_train;
    if (n_train < 10 || n_test < 5) { n_train = n / 2; n_test = n - n_train; }

    /* Build embedding matrices */
    Matrix* Emb_y = build_embedding(y, length, lag);
    Matrix* Emb_xy = mat_create(n, 2 * lag);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < lag; j++) {
            Emb_xy->data[i * 2 * lag + j] = y[i + j];
            Emb_xy->data[i * 2 * lag + lag + j] = x[i + j];
        }
    }

    /* Prepare targets (next value of Y) */
    double* targets = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) targets[i] = y[i + lag];

    /* Extract train/test splits */
    Matrix* train_y = mat_create(n_train, lag);
    Matrix* train_xy = mat_create(n_train, 2 * lag);
    Matrix* test_xy = mat_create(n_test, 2 * lag);
    double* train_t = malloc((size_t)n_train * sizeof(double));
    double* test_t = malloc((size_t)n_test * sizeof(double));
    for (int i = 0; i < n_train; i++) {
        for (int j = 0; j < lag; j++) train_y->data[i * lag + j] = Emb_y->data[i * lag + j];
        for (int j = 0; j < 2 * lag; j++) train_xy->data[i * 2 * lag + j] = Emb_xy->data[i * 2 * lag + j];
        train_t[i] = targets[i];
    }
    for (int i = 0; i < n_test; i++) {
        for (int j = 0; j < 2 * lag; j++)
            test_xy->data[i * 2 * lag + j] = Emb_xy->data[(n_train + i) * 2 * lag + j];
        test_t[i] = targets[n_train + i];
    }

    /* Compute kernel matrices */
    Matrix* K_yy = kernel_rbf(train_y, sigma);
    Matrix* K_xy = kernel_rbf(train_xy, sigma);

    /* Cross-kernel matrices for prediction */
    Matrix* K_yy_test = mat_create(n_test, n_train);
    Matrix* K_xy_test = mat_create(n_test, n_train);
    double gamma = 1.0 / (2.0 * sigma * sigma);
    for (int i = 0; i < n_test; i++) {
        for (int j = 0; j < n_train; j++) {
            double d2_y = 0.0, d2_xy = 0.0;
            for (int k = 0; k < lag; k++) {
                double diff_y = test_xy->data[i * 2 * lag + k] - train_y->data[j * lag + k];
                d2_y += diff_y * diff_y;
            }
            for (int k = 0; k < 2 * lag; k++) {
                double diff_xy = test_xy->data[i * 2 * lag + k] - train_xy->data[j * 2 * lag + k];
                d2_xy += diff_xy * diff_xy;
            }
            K_yy_test->data[i * n_train + j] = exp(-d2_y * gamma);
            K_xy_test->data[i * n_train + j] = exp(-d2_xy * gamma);
        }
    }

    double lambda = 1e-4;
    double mse_restricted = kernel_prediction_error(K_yy, train_t, K_yy_test, test_t, n_train, n_test, lambda);
    double mse_unrestricted = kernel_prediction_error(K_xy, train_t, K_xy_test, test_t, n_train, n_test, lambda);

    if (mse_unrestricted < 1e-15) mse_unrestricted = 1e-15;
    r->test_statistic = (mse_restricted - mse_unrestricted) / mse_unrestricted;
    if (r->test_statistic < 0.0) r->test_statistic = 0.0;
    r->sigma = sigma;
    r->lag_x = lag;
    r->lag_y = lag;
    r->n_obs = n_test;
    r->significance_level = alpha;
    /* Approximate p-value via normal approximation for prediction improvement */
    double z = sqrt((double)n_test) * r->test_statistic / (1.0 + r->test_statistic);
    r->p_value = 2.0 * (1.0 - 0.5 * (1.0 + erf(-fabs(z) / sqrt(2.0))));
    if (r->p_value > 1.0) r->p_value = 1.0;
    if (r->p_value < 0.0) r->p_value = 0.0;
    r->is_significant = (r->p_value < alpha);

    mat_free(Emb_y); mat_free(Emb_xy); mat_free(train_y); mat_free(train_xy);
    mat_free(test_xy); mat_free(K_yy); mat_free(K_xy);
    mat_free(K_yy_test); mat_free(K_xy_test);
    free(targets); free(train_t); free(test_t);
    return 0;
}

void nlg_print(const NonlinearGrangerResult* r) {
    if (!r) return;
    printf("Nonlinear Granger: T=%.4f p=%.6f sigma=%.3f lag=%d sig=%s\n",
        r->test_statistic, r->p_value, r->sigma, r->lag_x,
        r->is_significant ? "YES" : "NO");
}