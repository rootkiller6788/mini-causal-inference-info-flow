#ifndef FEVD_H
#define FEVD_H
#include "ts_core.h"

/* Forecast Error Variance Decomposition (FEVD)
 * Theorem (Sims 1980): For a stationary VAR(p), the h-step-ahead
 * forecast error variance of variable i attributable to shock j is:
 *   θ_{ij}(h) = Σ_{l=0}^{h-1} (e_i' Φ_l P e_j)² / Σ_{l=0}^{h-1} e_i' Φ_l Σ Φ_l' e_i
 * where Φ_l are MA(∞) coefficient matrices, P is the Cholesky factor of Σ.
 *
 * Reference: Sims, C.A. (1980). Econometrica, 48(1):1-48.
 */

typedef struct {
    int n_vars;
    int horizon;
    double** matrix;
    double* own_shocks;
    int* causal_ordering;
} FEVDResult;

/* Cholesky decomposition: A = L * L' where L is lower triangular.
 * Complexity: O(n³).  Theorem: Golub & Van Loan (2013), Alg 4.2.1. */
Matrix* cholesky_decompose(const Matrix* A);

/* Orthogonalized impulse response: compute MA(∞) coefficients Φ_h.
 * Φ_0 = I, Φ_h = Σ_{i=1}^{p} A_i Φ_{h-i} for h≥1.
 * Complexity: O(horizon * p * dim³). */
Matrix** compute_irf(const VARModel* var, int horizon);

/* Compute FEVD for a fitted VAR model via Cholesky decomposition.
 * Complexity: O(horizon * p * dim³ + dim³). */
FEVDResult* fevd_compute(const VARModel* var, int horizon);

void fevd_free(FEVDResult* f);
void fevd_print(const FEVDResult* f, int step);

/* Generalized FEVD (Pesaran & Shin 1998): order-invariant decomposition. */
FEVDResult* fevd_generalized(const VARModel* var, int horizon);

#endif
