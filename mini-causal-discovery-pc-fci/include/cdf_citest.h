#ifndef CDF_CITEST_H
#define CDF_CITEST_H

/*
 * cdf_citest.h — Conditional Independence Tests
 *
 * The PC/FCI algorithms rely on statistical tests of conditional
 * independence (CI) from observational data.
 *
 * Tests implemented:
 *   1. Fisher's z-test for partial correlation (Gaussian data)
 *      H0: ρ_{XY|Z} = 0  vs  H1: ρ_{XY|Z} ≠ 0
 *
 *   2. Mutual information based CI test (discrete/discretized data)
 *      Using G² likelihood ratio test for contingency tables.
 *
 * For Gaussian data, the partial correlation ρ_{XY|Z} measures the
 * correlation between X and Y after removing the linear effects of Z.
 * Fisher's z-transform: z = 0.5·ln((1+ρ)/(1-ρ)) → N(0, 1/(N-|Z|-3))
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "cdf_core.h"

/* ── Test Result ─────────────────────────────────────────────────── */

typedef struct {
    double stat;         /* test statistic */
    double p_value;      /* p-value */
    double df;           /* degrees of freedom */
    bool   is_independent; /* true if H0 not rejected at level alpha */
} CdfCITestResult;

/* ── Fisher's Z-Test (Gaussian Partial Correlation) ──────────────── */

/**
 * Test X ⊥ Y | Z using partial correlation and Fisher's z-transform.
 *
 * @param ds     dataset (p variables × N samples, column-major)
 * @param x_idx  index of variable X
 * @param y_idx  index of variable Y
 * @param Z      [|Z|] indices of conditioning variables
 * @param nZ     size of conditioning set
 * @param alpha  significance level
 * @return       test result
 */
CdfCITestResult cdf_citest_fisher_z(const CdfDataset *ds,
                                     int x_idx, int y_idx,
                                     const int *Z, int nZ, double alpha);

/* ── Mutual Information / G² Test ────────────────────────────────── */

/**
 * Test X ⊥ Y | Z using mutual information (for discrete data).
 * Uses the G² likelihood ratio test.
 *
 * @param ds      dataset (integer-valued for discrete variables)
 * @param levels  [p] number of levels for each variable
 * @param x_idx   index of X
 * @param y_idx   index of Y
 * @param Z       conditioning set indices
 * @param nZ      conditioning set size
 * @param alpha   significance level
 */
CdfCITestResult cdf_citest_g2(const CdfDataset *ds,
                               const int *levels,
                               int x_idx, int y_idx,
                               const int *Z, int nZ, double alpha);

/* ── Partial Correlation ─────────────────────────────────────────── */

/** Compute the partial correlation ρ_{XY|Z} between X and Y given Z.
 *
 *  Uses the correlation matrix inversion method:
 *    ρ_{XY|Z} = -Ω_{XY} / sqrt(Ω_{XX} · Ω_{YY})
 *  where Ω = Σ^{-1} is the precision matrix of (X, Y, Z).
 *
 *  For efficient computation with large nZ, uses the sweep operator
 *  or recursive formula from the correlation matrix.
 */
double cdf_citest_partial_corr(const CdfDataset *ds,
                                int x_idx, int y_idx,
                                const int *Z, int nZ);

/* ── Correlation Matrix ──────────────────────────────────────────── */

/** Compute the p×p correlation matrix from the dataset.
 *  Uses column-major storage: corr[i*p + j] = corr(X_i, X_j). */
void cdf_citest_corr_matrix(const CdfDataset *ds, double *corr);

/** Compute correlation between two variables (marginal). */
double cdf_citest_correlation(const CdfDataset *ds, int x_idx, int y_idx);

/* ── Utilities ───────────────────────────────────────────────────── */

/** Fisher z-transform: z = 0.5 * ln((1+r)/(1-r)) */
double cdf_citest_fisher_z_transform(double r);

/** Inverse Fisher z-transform: r = (e^{2z} - 1) / (e^{2z} + 1) */
double cdf_citest_fisher_z_inverse(double z);

/** Compute p-value from standard normal z-score (two-sided). */
double cdf_citest_normal_pvalue(double z_score);

/** Chi-square survival function (approximation). */
double cdf_citest_chi2_survival(double x, double df);

/** Free test result (no dynamic allocation currently). */
void cdf_citest_result_free(CdfCITestResult *res);

/** Print test result. */
void cdf_citest_print(const CdfCITestResult *res);

#endif /* CDF_CITEST_H */