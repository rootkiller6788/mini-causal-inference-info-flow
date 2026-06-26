/* fevd.c -- Forecast Error Variance Decomposition for VAR models.
 *
 * The FEVD quantifies how much of the h-step-ahead forecast error
 * variance of each variable can be attributed to shocks in each
 * variable. Based on Sims (1980) and Lutkepohl (2005).
 *
 * Knowledge points:
 * L3: Cholesky decomposition of covariance matrix
 * L3: Impulse response function via MA(inf) representation
 * L4: Wold representation theorem
 * L5: FEVD algorithm from VAR coefficients
 * L6: Interpreting FEVD as causal contribution
 */

#include "fevd.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

Matrix* cholesky_decompose(const Matrix* A) {
    if (!A || A->rows != A->cols || A->rows == 0) return NULL;
    int n = A->rows;
    Matrix* L = mat_create(n, n);
    if (!L) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += L->data[i * n + k] * L->data[j * n + k];
            if (i == j) {
                double val = A->data[i * n + i] - sum;
                if (val <= 1e-15) { mat_free(L); return NULL; }
                L->data[i * n + i] = sqrt(val);
            } else {
                if (fabs(L->data[j * n + j]) < 1e-15) { mat_free(L); return NULL; }
                L->data[i * n + j] = (A->data[i * n + j] - sum) / L->data[j * n + j];
            }
        }
    }
    return L;
}

Matrix** compute_irf(const VARModel* var, int horizon) {
    if (!var || horizon <= 0) return NULL;
    int dim = var->dim, p = var->p;
    Matrix** irf = malloc((size_t)horizon * sizeof(Matrix*));
    if (!irf) return NULL;
    for (int h = 0; h < horizon; h++) {
        irf[h] = mat_create(dim, dim);
        if (!irf[h]) {
            for (int k = 0; k < h; k++) mat_free(irf[k]);
            free(irf); return NULL;
        }
    }
    for (int i = 0; i < dim; i++)
        irf[0]->data[i * dim + i] = 1.0;
    for (int h = 1; h < horizon; h++) {
        for (int i = 1; i <= p && i <= h; i++) {
            const Matrix* Ai = var->A[i - 1];
            const Matrix* Phi = irf[h - i];
            for (int r = 0; r < dim; r++)
                for (int c = 0; c < dim; c++)
                    for (int k = 0; k < dim; k++)
                        irf[h]->data[r * dim + c] +=
                            Ai->data[r * dim + k] * Phi->data[k * dim + c];
        }
    }
    return irf;
}

FEVDResult* fevd_compute(const VARModel* var, int horizon) {
    if (!var || horizon <= 0) return NULL;
    int dim = var->dim;
    FEVDResult* f = calloc(1, sizeof(FEVDResult));
    if (!f) return NULL;
    f->n_vars = dim;
    f->horizon = horizon;
    f->matrix = malloc((size_t)horizon * sizeof(double*));
    f->own_shocks = malloc((size_t)horizon * sizeof(double));
    f->causal_ordering = malloc((size_t)dim * sizeof(int));
    if (!f->matrix || !f->own_shocks || !f->causal_ordering) {
        fevd_free(f); return NULL;
    }
    for (int h = 0; h < horizon; h++) {
        f->matrix[h] = calloc((size_t)(dim * dim), sizeof(double));
        if (!f->matrix[h]) { fevd_free(f); return NULL; }
    }
    Matrix* Sigma = mat_create(dim, dim);
    if (!Sigma) { fevd_free(f); return NULL; }
    for (int i = 0; i < dim; i++)
        Sigma->data[i * dim + i] = var->sigma2;
    Matrix* P = cholesky_decompose(Sigma);
    if (!P) { mat_free(Sigma); fevd_free(f); return NULL; }
    Matrix** irf = compute_irf(var, horizon);
    if (!irf) { mat_free(Sigma); mat_free(P); fevd_free(f); return NULL; }
    for (int h_step = 1; h_step <= horizon; h_step++) {
        double* mse = calloc((size_t)(dim * dim), sizeof(double));
        for (int l = 0; l < h_step; l++) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    double theta = 0.0;
                    for (int k = 0; k < dim; k++)
                        theta += irf[l]->data[i * dim + k] * P->data[k * dim + j];
                    mse[i * dim + j] += theta * theta;
                }
            }
        }
        for (int i = 0; i < dim; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < dim; j++) row_sum += mse[i * dim + j];
            if (row_sum < 1e-15) row_sum = 1.0;
            for (int j = 0; j < dim; j++)
                f->matrix[h_step - 1][i * dim + j] = mse[i * dim + j] / row_sum;
        }
        f->own_shocks[h_step - 1] = 0.0;
        for (int i = 0; i < dim; i++)
            f->own_shocks[h_step - 1] += f->matrix[h_step - 1][i * dim + i];
        f->own_shocks[h_step - 1] /= (double)dim;
        free(mse);
    }
    for (int i = 0; i < dim; i++)
        f->causal_ordering[i] = i;
    for (int l = 0; l < horizon; l++) mat_free(irf[l]);
    free(irf);
    mat_free(Sigma);
    mat_free(P);
    return f;
}

FEVDResult* fevd_generalized(const VARModel* var, int horizon) {
    if (!var || horizon <= 0) return NULL;
    int dim = var->dim;
    FEVDResult* f = calloc(1, sizeof(FEVDResult));
    if (!f) return NULL;
    f->n_vars = dim;
    f->horizon = horizon;
    f->matrix = malloc((size_t)horizon * sizeof(double*));
    f->own_shocks = malloc((size_t)horizon * sizeof(double));
    f->causal_ordering = malloc((size_t)dim * sizeof(int));
    for (int h = 0; h < horizon; h++)
        f->matrix[h] = calloc((size_t)(dim * dim), sizeof(double));
    Matrix* Sigma = mat_create(dim, dim);
    for (int i = 0; i < dim; i++) Sigma->data[i * dim + i] = var->sigma2;
    Matrix** irf = compute_irf(var, horizon);
    if (!irf) { mat_free(Sigma); fevd_free(f); return NULL; }
    for (int h_step = 1; h_step <= horizon; h_step++) {
        double* mse = calloc((size_t)(dim * dim), sizeof(double));
        for (int l = 0; l < h_step; l++) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    double theta = 0.0;
                    for (int k = 0; k < dim; k++)
                        theta += irf[l]->data[i * dim + k] * Sigma->data[k * dim + j];
                    double sigma_jj = Sigma->data[j * dim + j];
                    if (sigma_jj < 1e-15) sigma_jj = 1.0;
                    theta /= sqrt(sigma_jj);
                    mse[i * dim + j] += theta * theta;
                }
            }
        }
        for (int i = 0; i < dim; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < dim; j++) row_sum += mse[i * dim + j];
            if (row_sum < 1e-15) row_sum = 1.0;
            for (int j = 0; j < dim; j++)
                f->matrix[h_step - 1][i * dim + j] = mse[i * dim + j] / row_sum;
        }
        free(mse);
    }
    for (int i = 0; i < dim; i++) f->causal_ordering[i] = i;
    for (int l = 0; l < horizon; l++) mat_free(irf[l]);
    free(irf);
    mat_free(Sigma);
    return f;
}

void fevd_free(FEVDResult* f) {
    if (!f) return;
    if (f->matrix) {
        for (int h = 0; h < f->horizon; h++) free(f->matrix[h]);
        free(f->matrix);
    }
    free(f->own_shocks);
    free(f->causal_ordering);
    free(f);
}

void fevd_print(const FEVDResult* f, int step) {
    if (!f || step < 0 || step >= f->horizon) return;
    printf("FEVD (horizon %d):\n", step + 1);
    for (int i = 0; i < f->n_vars; i++) {
        printf("  Var %d: ", i);
        for (int j = 0; j < f->n_vars; j++)
            printf("%.3f ", f->matrix[step][i * f->n_vars + j]);
        printf("\n");
    }
}