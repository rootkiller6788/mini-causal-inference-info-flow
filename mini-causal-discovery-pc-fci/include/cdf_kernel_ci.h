#ifndef CDF_KERNEL_CI_H
#define CDF_KERNEL_CI_H

/*
 * cdf_kernel_ci.h — Kernel-Based Conditional Independence Tests
 *
 * Traditional parametric CI tests (Fisher's z, G^2) assume Gaussian or
 * discrete distributions. Kernel-based tests use reproducing kernel
 * Hilbert space (RKHS) embeddings to detect arbitrary non-linear
 * dependencies without distributional assumptions.
 *
 * Implemented methods:
 *   1. HSIC (Hilbert-Schmidt Independence Criterion):
 *      Tests marginal independence X ⟂ Y by comparing kernel
 *      covariance to zero. Gretton et al. (2005, 2008).
 *
 *   2. KCI (Kernel Conditional Independence):
 *      Tests X ⟂ Y | Z using conditional kernel mean embeddings.
 *      Zhang et al. (2011), "Kernel-based Conditional Independence
 *      Test and Application in Causal Discovery". UAI.
 *
 * Key insight (KCI):
 *   Under H0: X ⟂ Y | Z, the conditional cross-covariance operator
 *   Σ_{XY|Z} = 0. The test statistic is based on the normalized
 *   Hilbert-Schmidt norm of the empirical estimate.
 *
 *   Test statistic: T_n = (1/n) Tr(K_{X|Z} K_{Y|Z})
 *   where K_{X|Z} = (I - K_Z(K_Z + n*ε*I)^{-1}) K_X (I - ...)
 *   with K_X, K_Y, K_Z being kernel matrices on X, Y, Z respectively.
 *
 * Reference: Gretton, Bousquet, Smola & Scholkopf (2005),
 *   "Measuring Statistical Dependence with Hilbert-Schmidt Norms". ALT.
 * Reference: Zhang, Peters, Janzing & Scholkopf (2011),
 *   "Kernel-based Conditional Independence Test". UAI.
 */

#include <stddef.h>
#include <stdbool.h>
#include "cdf_core.h"
#include "cdf_citest.h"

/* ── Kernel Types ──────────────────────────────────────────────────── */

typedef enum {
    CDF_KERNEL_GAUSSIAN = 0,  /* exp(-||x-y||^2 / (2*sigma^2)) */
    CDF_KERNEL_LAPLACE  = 1,  /* exp(-||x-y|| / sigma)           */
    CDF_KERNEL_LINEAR   = 2,  /* <x, y>                         */
    CDF_KERNEL_POLY     = 3   /* (<x,y> + c)^d                  */
} CdfKernelType;

/* ── Kernel Configuration ──────────────────────────────────────────── */

typedef struct {
    CdfKernelType type;     /* kernel function type */
    double        sigma;    /* bandwidth (Gaussian/Laplace), 0=median heuristic */
    double        poly_c;   /* constant c for polynomial kernel */
    double        poly_d;   /* degree d for polynomial kernel */
    double        reg_eps;  /* regularization epsilon for KCI */
    int           n_perm;   /* number of permutations for null distribution */
    bool          use_perm; /* use permutation test (true) or gamma approx (false) */
} CdfKernelCIConfig;

/* ── Kernel Matrix ─────────────────────────────────────────────────── */

typedef struct {
    double *data;   /* [n * n] kernel gram matrix, row-major */
    int     n;      /* number of samples */
} CdfKernelMatrix;

/* ── KCI API ───────────────────────────────────────────────────────── */

/**
 * Create default kernel CI configuration.
 * Defaults: Gaussian kernel, median heuristic bandwidth, epsilon=1e-3,
 * 200 permutations, permutation test enabled.
 */
CdfKernelCIConfig cdf_kernel_ci_config_default(void);

/**
 * Build a kernel Gram matrix for the given variable.
 *
 * @param ds     dataset
 * @param var    variable index
 * @param config kernel configuration
 * @return       kernel matrix (caller must cdf_kernel_matrix_free)
 */
CdfKernelMatrix* cdf_kernel_build_matrix(const CdfDataset *ds,
                                          int var,
                                          const CdfKernelCIConfig *config);

/** Free kernel matrix. */
void cdf_kernel_matrix_free(CdfKernelMatrix *K);

/**
 * HSIC test: test marginal independence X ⟂ Y.
 * Uses the Hilbert-Schmidt Independence Criterion (Gretton et al., 2005).
 *
 * HSIC(X,Y) = (1/n^2) Tr(K_X H K_Y H) where H = I - (1/n)11^T.
 * Under H0, n*HSIC converges to a weighted sum of chi-squares.
 * p-value computed via permutation or gamma approximation.
 *
 * Complexity: O(n^3) for matrix operations or O(B*n^2) for permutation.
 *
 * @param ds      dataset
 * @param x_idx   index of X
 * @param y_idx   index of Y
 * @param config  kernel config
 * @param alpha   significance level
 * @return        CI test result
 */
CdfCITestResult cdf_kernel_hsic_test(const CdfDataset *ds,
                                      int x_idx, int y_idx,
                                      const CdfKernelCIConfig *config,
                                      double alpha);

/**
 * KCI test: test conditional independence X ⟂ Y | Z.
 * Uses the Kernel Conditional Independence test (Zhang et al., 2011).
 *
 * The test statistic is computed from the centered kernel matrices
 * of X and Y conditioned on Z:
 *   Σ_{X|Z} = K_X - K_X K_Z(K_Z + n*ε*I)^{-1}K_Z + ...
 *
 * Complexity: O(n^3) primarily due to matrix inversion of K_Z.
 *
 * @param ds      dataset
 * @param x_idx   X
 * @param y_idx   Y
 * @param Z       conditioning set indices
 * @param nZ      size of Z
 * @param config  kernel config
 * @param alpha   significance level
 * @return        CI test result
 */
CdfCITestResult cdf_kernel_kci_test(const CdfDataset *ds,
                                     int x_idx, int y_idx,
                                     const int *Z, int nZ,
                                     const CdfKernelCIConfig *config,
                                     double alpha);

/**
 * Median distance heuristic for Gaussian kernel bandwidth.
 *
 * σ = median{ ||x_i - x_j|| : i < j } / sqrt(2)
 *
 * This heuristic (Gretton et al.) adapts the bandwidth to the
 * data scale and is the default when sigma=0 is specified.
 *
 * @param ds    dataset
 * @param var   variable index
 * @return      bandwidth sigma
 */
double cdf_kernel_median_heuristic(const CdfDataset *ds, int var);

/**
 * Centering matrix: compute H K H where H = I - (1/n)11^T.
 * Used in HSIC test statistic computation.
 *
 * @param K  [n*n] kernel matrix, modified in-place
 * @param n  number of samples
 */
void cdf_kernel_center_matrix(double *K, int n);

/**
 * Permutation-based p-value for kernel CI tests.
 * Recomputes test statistic under random permutations of X or Y.
 *
 * @param ds          dataset
 * @param x_idx       X variable
 * @param y_idx       Y variable
 * @param Z           conditioning set
 * @param nZ          size of Z
 * @param config      kernel config
 * @param orig_stat   test statistic on original (unpermuted) data
 * @param is_cond     true=KCI, false=HSIC
 * @return            p-value from permutation test
 */
double cdf_kernel_permutation_pvalue(const CdfDataset *ds,
                                      int x_idx, int y_idx,
                                      const int *Z, int nZ,
                                      const CdfKernelCIConfig *config,
                                      double orig_stat, bool is_cond);

#endif /* CDF_KERNEL_CI_H */
