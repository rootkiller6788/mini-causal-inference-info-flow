/* partial_granger.c -- Partial Granger Causality and Partial Directed Coherence.
 *
 * Partial Granger causality extends conditional Granger by considering
 * the influence of exogenous inputs and latent confounding variables.
 *
 * Knowledge points:
 * L2: Partial Granger causality concept (controlling for exogenous variables)
 * L3: Partial directed coherence (PDC) in frequency domain
 * L5: PDC computation from VAR coefficients
 * L5: Statistical significance of PDC via surrogate data
 * L8: Advanced causal inference with confounders
 *
 * Reference: Guo, Seth, Kendrick, Zhou, Feng (2008), J. Neurosci. Methods, 172(1):79-93.
 *            Baccala & Sameshima (2001), Biol. Cybern., 84(6):463-474.
 */

#include "partial_granger.h"
#include "ts_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

PDCResult* pdc_create(int dim, int n_freqs) {
    PDCResult* p = calloc(1, sizeof(PDCResult));
    if (!p) return NULL;
    p->n_vars = dim;
    p->n_freqs = n_freqs;
    p->frequencies = malloc((size_t)n_freqs * sizeof(double));
    p->pdc = malloc((size_t)n_freqs * sizeof(double*));
    for (int f = 0; f < n_freqs; f++)
        p->pdc[f] = calloc((size_t)(dim * dim), sizeof(double));
    return p;
}

void pdc_free(PDCResult* p) {
    if (!p) return;
    free(p->frequencies);
    for (int f = 0; f < p->n_freqs; f++) free(p->pdc[f]);
    free(p->pdc);
    free(p);
}

/* Compute PDC from a fitted VAR model.
 * For each frequency omega, compute Ā(omega), then normalize columns.
 * Complexity: O(n_freqs * dim^3). */
int pdc_compute(PDCResult* pdc, const VARModel* var, int n_freqs) {
    if (!pdc || !var) return -1;
    int d = var->dim, p = var->p;
    for (int fidx = 0; fidx < n_freqs; fidx++) {
        double omega = M_PI * fidx / (n_freqs - 1 + 1e-10);
        pdc->frequencies[fidx] = omega;
        /* Build Ā(omega) = I - sum A_k e^{-i k omega} */
        double* A_real = calloc((size_t)(d * d), sizeof(double));
        double* A_imag = calloc((size_t)(d * d), sizeof(double));
        for (int i = 0; i < d; i++) A_real[i * d + i] = 1.0;
        for (int k = 0; k < p; k++) {
            double cr = cos(omega * (k + 1));
            double ci = -sin(omega * (k + 1));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++) {
                    double a = mat_get(var->A[k], i, j);
                    A_real[i * d + j] -= a * cr;
                    A_imag[i * d + j] -= a * ci;
                }
        }
        /* For each source j, compute PDC_{i<-j} */
        for (int j = 0; j < d; j++) {
            double denom = 0.0;
            for (int i = 0; i < d; i++) {
                double ar = A_real[i * d + j], ai = A_imag[i * d + j];
                denom += ar * ar + ai * ai;
            }
            if (denom < 1e-15) denom = 1.0;
            double inv_sqrt = 1.0 / sqrt(denom);
            for (int i = 0; i < d; i++) {
                double ar = A_real[i * d + j], ai = A_imag[i * d + j];
                pdc->pdc[fidx][i * d + j] = sqrt(ar * ar + ai * ai) * inv_sqrt;
            }
        }
        free(A_real); free(A_imag);
    }
    pdc->max_lag = p;
    return 0;
}

void pdc_print(const PDCResult* pdc, int freq_idx) {
    if (!pdc || freq_idx < 0 || freq_idx >= pdc->n_freqs) return;
    printf("PDC at omega=%.3f:\n", pdc->frequencies[freq_idx]);
    for (int i = 0; i < pdc->n_vars; i++) {
        printf("  Row %d: ", i);
        for (int j = 0; j < pdc->n_vars; j++)
            printf("%.3f ", pdc->pdc[freq_idx][i * pdc->n_vars + j]);
        printf("\n");
    }
}

PDCBands* pdc_band_average(const PDCResult* pdc) {
    if (!pdc) return NULL;
    PDCBands* b = calloc(1, sizeof(PDCBands));
    b->dim = pdc->n_vars;
    int d = pdc->n_vars;
    double band_limits[] = {0.05, 0.15, 0.30, 0.50, M_PI};
    int counts[5] = {0};
    for (int f = 0; f < pdc->n_freqs; f++) {
        double omega = pdc->frequencies[f];
        int band = 0;
        for (int bi = 1; bi < 5; bi++)
            if (omega <= band_limits[bi]) { band = bi - 1; break; }
            else if (bi == 4) band = 4;
        double* target = NULL;
        if (band == 0) target = b->delta;
        else if (band == 1) target = b->theta;
        else if (band == 2) target = b->alpha;
        else if (band == 3) target = b->beta;
        else target = b->gamma;
        for (int i = 0; i < d * d; i++)
            target[i] += pdc->pdc[f][i];
        counts[band]++;
    }
    for (int i = 0; i < d * d; i++) {
        if (counts[0] > 0) b->delta[i] /= counts[0];
        if (counts[1] > 0) b->theta[i] /= counts[1];
        if (counts[2] > 0) b->alpha[i] /= counts[2];
        if (counts[3] > 0) b->beta[i] /= counts[3];
        if (counts[4] > 0) b->gamma[i] /= counts[4];
    }
    return b;
}

void pdc_bands_free(PDCBands* b) { free(b); }

void pdc_bands_print(const PDCBands* b) {
    if (!b) return;
    printf("PDC Band Averages:\n");
    printf("  Delta: "); for(int i=0;i<b->dim;i++){for(int j=0;j<b->dim;j++)printf("%.2f ",b->delta[i*b->dim+j]);printf("| ");}
    printf("\n  Theta: "); for(int i=0;i<b->dim;i++){for(int j=0;j<b->dim;j++)printf("%.2f ",b->theta[i*b->dim+j]);printf("| ");}
    printf("\n  Alpha: "); for(int i=0;i<b->dim;i++){for(int j=0;j<b->dim;j++)printf("%.2f ",b->alpha[i*b->dim+j]);printf("| ");}
    printf("\n  Beta:  "); for(int i=0;i<b->dim;i++){for(int j=0;j<b->dim;j++)printf("%.2f ",b->beta[i*b->dim+j]);printf("| ");}
    printf("\n  Gamma: "); for(int i=0;i<b->dim;i++){for(int j=0;j<b->dim;j++)printf("%.2f ",b->gamma[i*b->dim+j]);printf("| ");}
    printf("\n");
}

PartialGrangerResult* pg_create(void) { return calloc(1, sizeof(PartialGrangerResult)); }
void pg_free(PartialGrangerResult* pg) { free(pg); }

/* Remove Z's influence from X via linear regression:
 * X* = X - Z*(Z'Z)^{-1}Z'X (residuals after projecting X onto Z). */
static void orthogonalize(const double* x, const double* z, int n,
    double* x_star) {
    /* Simple: regress x on [1, z] and take residuals */
    double sx = 0, sz = 0, sxz = 0, szz = 0;
    for (int i = 0; i < n; i++) {
        sx += x[i]; sz += z[i]; sxz += x[i] * z[i]; szz += z[i] * z[i];
    }
    double denom = n * szz - sz * sz;
    if (fabs(denom) < 1e-15) {
        memcpy(x_star, x, (size_t)n * sizeof(double));
        return;
    }
    double b1 = (n * sxz - sx * sz) / denom;
    double b0 = (sx - b1 * sz) / n;
    for (int i = 0; i < n; i++)
        x_star[i] = x[i] - (b0 + b1 * z[i]);
}

int pg_test(PartialGrangerResult* pg, const double* x, const double* y,
    const double* z, int length, int max_lag, double alpha) {
    if (!pg || !x || !y || !z || length < max_lag + 4) return -1;
    /* Orthogonalize X and Y w.r.t. Z */
    double* x_star = malloc((size_t)length * sizeof(double));
    double* y_star = malloc((size_t)length * sizeof(double));
    orthogonalize(x, z, length, x_star);
    orthogonalize(y, z, length, y_star);
    /* Now use standard Granger test on orthogonalized series */
    int p = max_lag, n = length - p, k_r = 1 + p, k_ur = 1 + 2 * p;
    Matrix* X_r = mat_create(n, k_r);
    Matrix* X_ur = mat_create(n, k_ur);
    double* Y_vec = malloc((size_t)n * sizeof(double));
    double* beta = malloc((size_t)k_ur * sizeof(double));
    double* resid = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        X_r->data[i * k_r + 0] = 1.0; X_ur->data[i * k_ur + 0] = 1.0;
        for (int lag = 0; lag < p; lag++) {
            X_r->data[i * k_r + 1 + lag] = y_star[i + p - 1 - lag];
            X_ur->data[i * k_ur + 1 + lag] = y_star[i + p - 1 - lag];
            X_ur->data[i * k_ur + 1 + p + lag] = x_star[i + p - 1 - lag];
        }
        Y_vec[i] = y_star[i + p];
    }
    pg->rss_restricted = ols_estimate(X_r, Y_vec, beta, resid, n, k_r);
    pg->rss_unrestricted = ols_estimate(X_ur, Y_vec, beta, resid, n, k_ur);
    if (pg->rss_unrestricted < 1e-15) pg->rss_unrestricted = 1e-15;
    double f = ((pg->rss_restricted - pg->rss_unrestricted) / p) /
               (pg->rss_unrestricted / (n - k_ur));
    if (f < 0) f = 0;
    pg->f_statistic = f;
    pg->lag_order = p;
    pg->n_obs = n;
    pg->significance_level = alpha;
    /* F-distribution p-value using approximation */
    double xv = (double)(n - k_ur) / (n - k_ur + p * f);
    pg->p_value = 1.0 - pow(xv, (n - k_ur) / 2.0);
    if (pg->p_value > 1.0) pg->p_value = 1.0;
    pg->is_significant = (pg->p_value < alpha);
    mat_free(X_r); mat_free(X_ur);
    free(Y_vec); free(beta); free(resid);
    free(x_star); free(y_star);
    return 0;
}

void pg_print(const PartialGrangerResult* pg) {
    if (!pg) return;
    printf("Partial Granger: F=%.4f p=%.4f sig=%s (lag=%d)\n",
        pg->f_statistic, pg->p_value, pg->is_significant ? "YES" : "NO",
        pg->lag_order);
}