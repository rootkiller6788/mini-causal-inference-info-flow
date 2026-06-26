/* spectral_granger.c -- Frequency-domain Granger Causality.
 *
 * Decomposes Granger causality into frequency components using
 * the Fourier transform of VAR coefficients. Based on Geweke (1982).
 *
 * The spectral Granger causality from X to Y at frequency omega is:
 *   f_{X->Y}(omega) = log( |S_yy(omega)| / |S_yy - H_yx Sigma_xx H_yx*| )
 * where H(omega) = (I - sum_{k=1}^p A_k e^{-i k omega})^{-1}
 *
 * Knowledge points:
 * L2: Frequency-domain causality (spectral Granger)
 * L3: Transfer function matrix H(omega) from VAR coefficients
 * L3: Complex matrix operations for spectral analysis
 * L4: Geweke decomposition theorem
 * L5: FFT-based spectral computation
 * L6: Frequency-specific causal patterns
 *
 * Reference: Geweke, J. (1982). JASA, 77(378):304-313.
 *            Geweke, J. (1984). JASA, 79(388):907-915.
 */

#include "spectral_granger.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double re, im; } Complex;

static Complex __attribute__((unused)) c_add(Complex a, Complex b) { Complex r = { a.re + b.re, a.im + b.im }; return r; }
static Complex c_sub(Complex a, Complex b) { Complex r = { a.re - b.re, a.im - b.im }; return r; }
static Complex c_mul(Complex a, Complex b) {
    Complex r = { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
    return r;
}
static Complex c_exp_i(double theta) {
    Complex r = { cos(theta), sin(theta) };
    return r;
}
static double c_abs(Complex a) { return sqrt(a.re * a.re + a.im * a.im); }
static Complex c_inv(Complex a) {
    double den = a.re * a.re + a.im * a.im;
    if (den < 1e-15) { Complex z = { 1e10, 0 }; return z; }
    Complex r = { a.re / den, -a.im / den };
    return r;
}

/* Compute H(omega) = (I - sum_{k=1}^p A_k * e^{-i k omega})^{-1}
 * This is the transfer function matrix of the VAR at frequency omega.
 * dim: dimension of VAR, p: lag order, A: coefficient matrices.
 * Returns dim x dim complex matrix (stored as 2 * dim * dim doubles: [re, im] interleaved). */
static Complex* transfer_function(const VARModel* var, double omega) {
    int d = var->dim, p = var->p;
    Complex* H = calloc((size_t)(d * d), sizeof(Complex));
    /* I_d - sum_{k=1}^p A_k * e^{-i k omega} */
    for (int i = 0; i < d; i++) {
        H[i * d + i].re = 1.0;
        H[i * d + i].im = 0.0;
    }
    for (int k = 0; k < p; k++) {
        Complex e_ik = c_exp_i(-omega * (k + 1));
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++) {
                Complex a_ij = { mat_get(var->A[k], i, j), 0.0 };
                Complex term = c_mul(a_ij, e_ik);
                H[i * d + j] = c_sub(H[i * d + j], term);
            }
    }
    /* Invert the matrix (complex Gaussian elimination) */
    Complex* inv = calloc((size_t)(d * d), sizeof(Complex));
    /* Build augmented matrix [H | I] */
    int n = d;
    Complex* aug = calloc((size_t)(n * 2 * n), sizeof(Complex));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i * (2 * n) + j] = H[i * n + j];
        aug[i * (2 * n) + n + i].re = 1.0;
    }
    /* Gauss-Jordan elimination */
    for (int i = 0; i < n; i++) {
        Complex piv = aug[i * (2 * n) + i];
        if (c_abs(piv) < 1e-12) {
            int sw = -1;
            for (int k = i + 1; k < n; k++)
                if (c_abs(aug[k * (2 * n) + i]) > 1e-12) { sw = k; break; }
            if (sw < 0) { free(aug); free(H); return NULL; }
            for (int j = 0; j < 2 * n; j++) {
                Complex tmp = aug[i * (2 * n) + j];
                aug[i * (2 * n) + j] = aug[sw * (2 * n) + j];
                aug[sw * (2 * n) + j] = tmp;
            }
            piv = aug[i * (2 * n) + i];
        }
        Complex inv_piv = c_inv(piv);
        for (int j = 0; j < 2 * n; j++)
            aug[i * (2 * n) + j] = c_mul(aug[i * (2 * n) + j], inv_piv);
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            Complex factor = aug[k * (2 * n) + i];
            for (int j = 0; j < 2 * n; j++)
                aug[k * (2 * n) + j] = c_sub(aug[k * (2 * n) + j],
                    c_mul(factor, aug[i * (2 * n) + j]));
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i * n + j] = aug[i * (2 * n) + n + j];
    free(aug); free(H);
    return inv;
}

SpectralGranger* sg_create(int n_freqs) {
    SpectralGranger* sg = calloc(1, sizeof(SpectralGranger));
    if (sg) {
        sg->n_freqs = n_freqs;
        sg->frequencies = malloc((size_t)n_freqs * sizeof(double));
        sg->gc_x_to_y = malloc((size_t)n_freqs * sizeof(double));
        sg->gc_y_to_x = malloc((size_t)n_freqs * sizeof(double));
    }
    return sg;
}

void sg_free(SpectralGranger* sg) {
    if (!sg) return;
    free(sg->frequencies); free(sg->gc_x_to_y); free(sg->gc_y_to_x);
    free(sg);
}

/* Fit a bivariate VAR to x and y for spectral Granger computation.
 * Returns a VARModel* or NULL on failure. */
static VARModel* fit_bivariate_var(const double* x, const double* y, int length, int p) {
    TimeSeries* ts = ts_create(length, 2);
    if (!ts) return NULL;
    for (int i = 0; i < length; i++) {
        ts->values[i * 2] = x[i];
        ts->values[i * 2 + 1] = y[i];
    }
    VARModel* var = var_create(2, p);
    if (!var) { ts_free(ts); return NULL; }
    var_fit(var, ts);
    ts_free(ts);
    return var;
}

/* Geweke decomposition: separate total interdependence into
 * X->Y, Y->X, and instantaneous components in frequency domain.
 *
 * Total: f_{X,Y}(omega) = log(S_xx * S_yy / |S|)
 * Directional: f_{X->Y} = log(S_yy / (S_yy - H_yx Sigma_xx H_yx*))
 *
 * For a simple implementation, we use the VAR-based approach:
 * compute transfer function H(omega), then derive spectral GC. */
int sg_geweke_decomposition(SpectralGranger* sg,
    const double* x, const double* y, int length, int n_freqs, int p) {
    if (!sg || !x || !y || length < p + 2 || n_freqs < 1) return -1;

    VARModel* var = fit_bivariate_var(x, y, length, p);
    if (!var) return -1;

    int d = var->dim;
    double sigma_xx = var->sigma2;
    double sigma_yy = var->sigma2;
    (void)sigma_yy; /* reserved for full Geweke decomposition */

    for (int fidx = 0; fidx < n_freqs; fidx++) {
        double omega = M_PI * fidx / (n_freqs - 1 + 1e-10);
        sg->frequencies[fidx] = omega;
        Complex* H = transfer_function(var, omega);
        if (!H) {
            sg->gc_x_to_y[fidx] = 0.0;
            sg->gc_y_to_x[fidx] = 0.0;
            continue;
        }
        /* Spectral density components via H(omega) * Sigma * H*(omega) */
        Complex H_xy = H[0 * d + 1];  /* transfer from y to x */
        Complex H_yx = H[1 * d + 0];  /* transfer from x to y */
        /* S_yy(omega) approx = sigma_yy * |H_yy|^2 at this frequency */
        double S_yy = sigma_yy * c_abs(H[1 * d + 1]) * c_abs(H[1 * d + 1]) + 1e-10;
        double H_yx_sq = c_abs(H_yx) * c_abs(H_yx);
        double f_xy = log(S_yy / (S_yy - sigma_xx * H_yx_sq + 1e-10));
        double H_xy_sq = c_abs(H_xy) * c_abs(H_xy);
        double f_yx = log(sigma_xx * c_abs(H[0 * d + 0]) * c_abs(H[0 * d + 0]) /
            (sigma_xx * c_abs(H[0 * d + 0]) * c_abs(H[0 * d + 0]) -
             sigma_yy * H_xy_sq + 1e-10));
        sg->gc_x_to_y[fidx] = f_xy > 0 ? f_xy : 0.0;
        sg->gc_y_to_x[fidx] = f_yx > 0 ? f_yx : 0.0;
        free(H);
    }
    var_free(var);
    return 0;
}

int sg_compute(SpectralGranger* sg, const double* x, const double* y,
    int length, int n_freqs, int max_lag) {
    return sg_geweke_decomposition(sg, x, y, length, n_freqs, max_lag);
}

void sg_print(const SpectralGranger* sg) {
    if (!sg) return;
    printf("Spectral Granger: %d frequencies\n", sg->n_freqs);
    for (int i = 0; i < sg->n_freqs && i < 10; i++)
        printf("  omega=%.3f: gc(x->y)=%.4f gc(y->x)=%.4f\n",
            sg->frequencies[i], sg->gc_x_to_y[i], sg->gc_y_to_x[i]);
    if (sg->n_freqs > 10) printf("  ... (%d more)\n", sg->n_freqs - 10);
}