#ifndef NONLINEAR_GRANGER_H
#define NONLINEAR_GRANGER_H
#include "ts_core.h"

/* Nonlinear Granger Causality via Kernel Methods
 *
 * Definition (Ancona et al. 2004): X nonlinearly Granger-causes Y if
 * past values of X contain information about future Y beyond what is
 * contained in past Y alone, where the relationship may be nonlinear.
 *
 * Approach: Map data into reproducing kernel Hilbert space (RKHS),
 * compute restricted and unrestricted regressions, compare residuals.
 *
 * Kernel: Gaussian RBF K(x,y) = exp(-||x-y||² / (2σ²))
 *
 * Reference: Ancona, Marinazzo, Stramaglia (2004), PRE 70(5):056221.
 */

typedef struct {
    double test_statistic;
    double p_value;
    double sigma;
    int lag_x, lag_y;
    int n_obs;
    bool is_significant;
    double significance_level;
} NonlinearGrangerResult;

/* Compute Gaussian RBF kernel matrix K_ij = exp(-||x_i - x_j||² / (2σ²)).
 * Complexity: O(n² * dim). Returns n×n matrix. */
Matrix* kernel_rbf(const Matrix* X, double sigma);

/* Compute kernel ridge regression: α = (K + λI)⁻¹ y.
 * Returns coefficient vector α of length n.
 * Complexity: O(n³) due to matrix inverse. */
double* kernel_ridge_regression(const Matrix* K, const double* y, int n, double lambda);

/* Test nonlinear Granger causality X→Y using kernel methods.
 * Fits kernel ridge regressions with and without X's history,
 * then computes an F-like test statistic on the residuals.
 * Complexity: O(max_lag * n³) for kernel matrix inversion. */
NonlinearGrangerResult* nlg_create(void);
void nlg_free(NonlinearGrangerResult* r);
int nlg_test(NonlinearGrangerResult* r, const double* x, const double* y,
    int length, int max_lag, double sigma, double alpha);
void nlg_print(const NonlinearGrangerResult* r);

/* Select optimal kernel bandwidth σ via median heuristic:
 * σ = median{||x_i - x_j|| : i < j} / √2
 * Complexity: O(n²) pairwise distances. */
double kernel_median_heuristic(const double* x, int n);

#endif
