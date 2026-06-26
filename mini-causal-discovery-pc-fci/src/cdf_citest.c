/* cdf_citest.c — Conditional Independence Tests
 *
 * Implements Fisher's z-test for partial correlation and G² test
 * for discrete data. These are the statistical engines of PC/FCI.
 *
 * Fisher's z-test:
 *   rho = partial correlation ρ_{XY|Z}
 *   z = 0.5 * ln((1+ρ)/(1-ρ)) ~ N(0, 1/(N-|Z|-3))
 *   test statistic = z * sqrt(N - |Z| - 3) ~ N(0,1)
 *
 * Partial correlation via precision matrix:
 *   Let Σ be the (2+|Z|)×(2+|Z|) covariance matrix of (X,Y,Z).
 *   Let Ω = Σ^{-1} be the precision matrix.
 *   Then ρ_{XY|Z} = -Ω_{01} / sqrt(Ω_{00} * Ω_{11})
 */

#include "cdf_citest.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ── Helper: standard normal CDF (Abramowitz & Stegun approximation) ─ */
static double normal_cdf(double x)
{
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    int sign = (x < 0) ? -1 : 1;
    x = fabs(x) / sqrt(2.0);

    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * exp(-x*x);

    return 0.5 * (1.0 + sign * y);
}

/* ──────────────── Correlation Matrix ─────────────────────────────── */

void cdf_citest_corr_matrix(const CdfDataset *ds, double *corr)
{
    if (!ds || !corr || ds->p < 1) return;

    int N = ds->N, p = ds->p;

    /* Compute means */
    double *means = (double*)calloc((size_t)p, sizeof(double));
    if (!means) return;

    for (int i = 0; i < p; i++)
        for (int j = 0; j < N; j++)
            means[i] += ds->data[i * N + j];
    for (int i = 0; i < p; i++)
        means[i] /= (double)N;

    /* Compute standard deviations */
    double *stds = (double*)calloc((size_t)p, sizeof(double));
    if (!stds) { free(means); return; }

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < N; j++) {
            double d = ds->data[i * N + j] - means[i];
            stds[i] += d * d;
        }
        stds[i] = sqrt(stds[i] / (double)(N - 1));
    }

    /* Compute correlations */
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            if (i == j) {
                corr[i * p + j] = 1.0;
                continue;
            }
            if (stds[i] < CDF_EPS || stds[j] < CDF_EPS) {
                corr[i * p + j] = 0.0;
                continue;
            }
            double cov = 0.0;
            for (int k = 0; k < N; k++)
                cov += (ds->data[i * N + k] - means[i]) *
                       (ds->data[j * N + k] - means[j]);
            cov /= (double)(N - 1);
            corr[i * p + j] = cov / (stds[i] * stds[j]);
        }
    }

    free(means); free(stds);
}

double cdf_citest_correlation(const CdfDataset *ds, int x_idx, int y_idx)
{
    if (!ds || x_idx < 0 || y_idx < 0 || x_idx >= ds->p || y_idx >= ds->p)
        return 0.0;

    int N = ds->N;
    double mx = 0.0, my = 0.0;
    for (int i = 0; i < N; i++) {
        mx += ds->data[x_idx * N + i];
        my += ds->data[y_idx * N + i];
    }
    mx /= N; my /= N;

    double sx = 0.0, sy = 0.0, sxy = 0.0;
    for (int i = 0; i < N; i++) {
        double dx = ds->data[x_idx * N + i] - mx;
        double dy = ds->data[y_idx * N + i] - my;
        sx += dx * dx;
        sy += dy * dy;
        sxy += dx * dy;
    }

    if (sx < CDF_EPS || sy < CDF_EPS) return 0.0;
    return sxy / sqrt(sx * sy);
}

/* ──────────────── Partial Correlation ────────────────────────────── */

/*
 * Partial correlation via precision matrix method.
 *
 * Let V = {X, Y} ∪ Z have size m = 2 + nZ.
 * Compute the m×m correlation matrix corr_m, then invert it.
 * The partial correlation is: -Ω_{0,1} / sqrt(Ω_{0,0} * Ω_{1,1})
 *
 * For numerical stability with large nZ, we use the correlation matrix
 * of the subset and invert via Gaussian elimination.
 */
double cdf_citest_partial_corr(const CdfDataset *ds,
                                int x_idx, int y_idx,
                                const int *Z, int nZ)
{
    if (!ds || nZ < 0 || nZ > CDF_MAX_NODES) return 0.0;

    int m = 2 + nZ;
    if (m > CDF_MAX_NODES) return 0.0;

    /* Build index list: [x_idx, y_idx, Z[0], ..., Z[nZ-1]] */
    int *idx = (int*)malloc((size_t)m * sizeof(int));
    if (!idx) return 0.0;
    idx[0] = x_idx;
    idx[1] = y_idx;
    for (int i = 0; i < nZ; i++)
        idx[2 + i] = Z[i];

    /* Build correlation submatrix */
    double *sub_corr = (double*)malloc((size_t)m * m * sizeof(double));
    if (!sub_corr) { free(idx); return 0.0; }

    int N = ds->N;
    /* Compute means and stds for the m variables */
    double *means = (double*)calloc((size_t)m, sizeof(double));
    double *stds  = (double*)calloc((size_t)m, sizeof(double));
    if (!means || !stds) { free(idx); free(sub_corr); free(means); free(stds); return 0.0; }

    for (int a = 0; a < m; a++) {
        int vi = idx[a];
        double sum = 0.0;
        for (int k = 0; k < N; k++)
            sum += ds->data[vi * N + k];
        means[a] = sum / N;
    }

    for (int a = 0; a < m; a++) {
        int vi = idx[a];
        double ss = 0.0;
        for (int k = 0; k < N; k++) {
            double d = ds->data[vi * N + k] - means[a];
            ss += d * d;
        }
        stds[a] = sqrt(ss / (N - 1));
    }

    for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) {
            if (a == b) { sub_corr[a * m + b] = 1.0; continue; }
            if (stds[a] < CDF_EPS || stds[b] < CDF_EPS) {
                sub_corr[a * m + b] = 0.0; continue;
            }
            int va = idx[a], vb = idx[b];
            double cov = 0.0;
            for (int k = 0; k < N; k++)
                cov += (ds->data[va * N + k] - means[a]) *
                       (ds->data[vb * N + k] - means[b]);
            cov /= (N - 1);
            sub_corr[a * m + b] = cov / (stds[a] * stds[b]);
        }
    }

    /* Invert sub_corr via Gaussian elimination with partial pivoting */
    /* Augment with identity: [corr | I] → [I | corr^{-1}] */
    double *aug = (double*)calloc((size_t)m * (2 * m), sizeof(double));
    if (!aug) { free(idx); free(sub_corr); free(means); free(stds); return 0.0; }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            aug[i * (2 * m) + j] = sub_corr[i * m + j];
        aug[i * (2 * m) + m + i] = 1.0;
    }

    /* Forward elimination */
    for (int col = 0; col < m; col++) {
        /* Pivot */
        int pivot = col;
        double maxv = fabs(aug[col * (2 * m) + col]);
        for (int row = col + 1; row < m; row++) {
            double v = fabs(aug[row * (2 * m) + col]);
            if (v > maxv) { maxv = v; pivot = row; }
        }
        if (maxv < CDF_EPS) { free(aug); free(idx); free(sub_corr); free(means); free(stds); return 0.0; }

        if (pivot != col) {
            for (int j = 0; j < 2 * m; j++) {
                double t = aug[col * (2 * m) + j];
                aug[col * (2 * m) + j] = aug[pivot * (2 * m) + j];
                aug[pivot * (2 * m) + j] = t;
            }
        }

        double piv = aug[col * (2 * m) + col];
        /* Normalize */
        for (int j = 0; j < 2 * m; j++)
            aug[col * (2 * m) + j] /= piv;

        /* Eliminate other rows */
        for (int row = 0; row < m; row++) {
            if (row == col) continue;
            double factor = aug[row * (2 * m) + col];
            for (int j = 0; j < 2 * m; j++)
                aug[row * (2 * m) + j] -= factor * aug[col * (2 * m) + j];
        }
    }

    /* Extract Ω = corr^{-1} from the right half */
    double omega_01 = aug[0 * (2 * m) + m + 1];
    double omega_00 = aug[0 * (2 * m) + m + 0];
    double omega_11 = aug[1 * (2 * m) + m + 1];

    double rho = -omega_01 / sqrt(fmax(omega_00 * omega_11, CDF_EPS));

    free(aug); free(idx); free(sub_corr); free(means); free(stds);
    return rho;
}

/* ──────────────── Fisher Z-Transform ─────────────────────────────── */

double cdf_citest_fisher_z_transform(double r)
{
    /* Clamp r to avoid log domain errors */
    if (r >= 1.0 - CDF_EPS) r = 1.0 - CDF_EPS;
    if (r <= -1.0 + CDF_EPS) r = -1.0 + CDF_EPS;
    return 0.5 * log((1.0 + r) / (1.0 - r));
}

double cdf_citest_fisher_z_inverse(double z)
{
    double e2z = exp(2.0 * z);
    return (e2z - 1.0) / (e2z + 1.0);
}

double cdf_citest_normal_pvalue(double z_score)
{
    /* Two-sided p-value */
    double p = 2.0 * (1.0 - normal_cdf(fabs(z_score)));
    return fmax(p, CDF_EPS);
}

/* ──────────────── Fisher's Z-Test ────────────────────────────────── */

CdfCITestResult cdf_citest_fisher_z(const CdfDataset *ds,
                                     int x_idx, int y_idx,
                                     const int *Z, int nZ, double alpha)
{
    CdfCITestResult res = {0.0, 1.0, 0.0, true};

    if (!ds || nZ < 0) return res;

    /* Compute partial correlation */
    double rho = cdf_citest_partial_corr(ds, x_idx, y_idx, Z, nZ);

    /* Fisher z-transform */
    double z_fisher = cdf_citest_fisher_z_transform(rho);

    /* Standard error: 1 / sqrt(N - nZ - 3) */
    int N = ds->N;
    double df_denom = (double)(N - nZ - 3);
    if (df_denom < 1.0) df_denom = 1.0;

    double se = 1.0 / sqrt(df_denom);
    double z_stat = z_fisher / se;  /* ~ N(0,1) under H0 */

    res.stat = z_stat;
    res.p_value = cdf_citest_normal_pvalue(z_stat);
    res.df = df_denom;
    res.is_independent = (res.p_value > alpha);

    return res;
}

/* ──────────────── G² Test (Simplified) ───────────────────────────── */

CdfCITestResult cdf_citest_g2(const CdfDataset *ds,
                               const int *levels,
                               int x_idx, int y_idx,
                               const int *Z, int nZ, double alpha)
{
    CdfCITestResult res = {0.0, 1.0, 0.0, true};

    if (!ds || !levels || nZ < 0) return res;

    int N = ds->N;
    int lx = levels[x_idx], ly = levels[y_idx];
    if (lx < 2 || ly < 2) return res;

    /* G² test for discrete data.
     * For nZ == 0: standard mutual information G² = 2N * I(X;Y).
     * For nZ > 0: conditional MI via stratification.
     * Uses chi-square approximation: 2N * I ~ χ²(df) where
     * df = (|X|-1)(|Y|-1) for marginal, scaled for conditional.
     */

    /* Compute total Z levels product for conditioning stratification */
    int nz_configs = 1;
    for (int i = 0; i < nZ; i++) {
        int iz = Z[i];
        if (iz < 0 || iz >= ds->p) continue;
        if (levels[iz] >= 2) {
            nz_configs *= levels[iz];
            if (nz_configs > 500) { nz_configs = 500; break; }
        }
    }

    if (nZ == 0) {
        /* Zero conditioning: standard mutual information G² */
        int *joint = (int*)calloc((size_t)lx * ly, sizeof(int));
        int *marg_x = (int*)calloc((size_t)lx, sizeof(int));
        int *marg_y = (int*)calloc((size_t)ly, sizeof(int));
        if (!joint || !marg_x || !marg_y) {
            free(joint); free(marg_x); free(marg_y); return res;
        }

        for (int k = 0; k < N; k++) {
            int xi = (int)ds->data[x_idx * N + k];
            int yi = (int)ds->data[y_idx * N + k];
            if (xi < 0 || xi >= lx || yi < 0 || yi >= ly) continue;
            joint[xi * ly + yi]++;
            marg_x[xi]++;
            marg_y[yi]++;
        }

        double mi = 0.0;
        for (int xi = 0; xi < lx; xi++) {
            for (int yi = 0; yi < ly; yi++) {
                int n_xy = joint[xi * ly + yi];
                if (n_xy == 0) continue;
                double expected = (double)marg_x[xi] * (double)marg_y[yi] / (double)N;
                if (expected < 1e-9) continue;
                mi += (double)n_xy * log((double)n_xy / expected);
            }
        }
        mi /= (double)N;

        double g2 = 2.0 * (double)N * mi;
        double df = (double)((lx - 1) * (ly - 1));

        res.stat = g2;
        res.df = df;
        res.p_value = cdf_citest_chi2_survival(g2, df);
        res.is_independent = (res.p_value > alpha);

        free(joint); free(marg_x); free(marg_y);
        return res;
    }

    /* For nZ > 0: stratify by Z and accumulate weighted G².
     * Encode Z-values into an integer index for bucketing.
     * Accumulate: G² = Σ_{z} 2 * n_z * I(X;Y|Z=z) */
    if (nz_configs > 100) {
        /* Too many Z configs: fall back to partial-correlation-based
         * approximation using Fisher z-transform */
        double rho = cdf_citest_partial_corr(ds, x_idx, y_idx, Z, nZ);
        double z_f = cdf_citest_fisher_z_transform(rho);
        double df_denom = (double)(N - nZ - 3);
        if (df_denom < 1.0) df_denom = 1.0;
        double se = 1.0 / sqrt(df_denom);
        double z_stat = z_f / se;
        res.stat = z_stat;
        res.p_value = cdf_citest_normal_pvalue(z_stat);
        res.df = df_denom;
        res.is_independent = (res.p_value > alpha);
        return res;
    }

    /* Build per-Z-configuration tables */
    int *z_idx_arr = (int*)malloc((size_t)N * sizeof(int));
    int *z_cnt = (int*)calloc((size_t)nz_configs, sizeof(int));
    if (!z_idx_arr || !z_cnt) {
        free(z_idx_arr); free(z_cnt); return res;
    }

    /* Assign each sample to a Z-configuration bucket */
    for (int k = 0; k < N; k++) {
        int zid = 0, stride = 1;
        int valid = 1;
        for (int zi = 0; zi < nZ && valid; zi++) {
            int iz = Z[zi];
            if (iz < 0 || iz >= ds->p || levels[iz] < 2) continue;
            int val = (int)ds->data[iz * N + k];
            if (val < 0 || val >= levels[iz]) { valid = 0; break; }
            zid += val * stride;
            stride *= levels[iz];
        }
        z_idx_arr[k] = valid ? zid : -1;
        if (valid && zid < nz_configs) z_cnt[zid]++;
    }

    /* Accumulate G² across Z-strata */
    double g2_total = 0.0;
    double df_total = 0.0;

    for (int zc = 0; zc < nz_configs; zc++) {
        int nz_samples = z_cnt[zc];
        if (nz_samples < 5) continue; /* too sparse */

        int *joint_z = (int*)calloc((size_t)lx * ly, sizeof(int));
        int *mx = (int*)calloc((size_t)lx, sizeof(int));
        int *my = (int*)calloc((size_t)ly, sizeof(int));
        if (!joint_z || !mx || !my) {
            free(joint_z); free(mx); free(my); continue;
        }

        for (int k = 0; k < N; k++) {
            if (z_idx_arr[k] != zc) continue;
            int xi = (int)ds->data[x_idx * N + k];
            int yi = (int)ds->data[y_idx * N + k];
            if (xi < 0 || xi >= lx || yi < 0 || yi >= ly) continue;
            joint_z[xi * ly + yi]++;
            mx[xi]++; my[yi]++;
        }

        double mi_z = 0.0;
        for (int xi = 0; xi < lx; xi++) {
            for (int yi = 0; yi < ly; yi++) {
                int n_xy = joint_z[xi * ly + yi];
                if (n_xy == 0) continue;
                double expected = (double)mx[xi] * (double)my[yi] / (double)nz_samples;
                if (expected < 1e-9) continue;
                mi_z += (double)n_xy * log((double)n_xy / expected);
            }
        }
        mi_z /= (double)nz_samples;

        g2_total += 2.0 * (double)nz_samples * mi_z;
        df_total += (double)((lx - 1) * (ly - 1));

        free(joint_z); free(mx); free(my);
    }

    if (df_total < 0.5)
        df_total = (double)((lx - 1) * (ly - 1));

    res.stat = g2_total;
    res.df = df_total;
    res.p_value = cdf_citest_chi2_survival(g2_total, df_total);
    res.is_independent = (res.p_value > alpha);

    free(z_idx_arr); free(z_cnt);
    return res;
}

/* ──────────────── Chi-square Survival ────────────────────────────── */

double cdf_citest_chi2_survival(double x, double df)
{
    /* Wilson-Hilferty approximation for chi-square survival */
    if (x < 0.0 || df < 0.1) return 1.0;

    double z = (pow(x / df, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df))) /
               sqrt(2.0 / (9.0 * df));

    return 1.0 - normal_cdf(z);
}

void cdf_citest_result_free(CdfCITestResult *res)
{
    (void)res;  /* no dynamic allocation */
}

void cdf_citest_print(const CdfCITestResult *res)
{
    if (!res) return;
    printf("CI Test: stat=%.4f, p=%.4f, df=%.1f, %s\n",
           res->stat, res->p_value, res->df,
           res->is_independent ? "INDEPENDENT" : "DEPENDENT");
}