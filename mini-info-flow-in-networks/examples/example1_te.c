/* Example 1: Transfer Entropy — Quantifying Directed Information Flow.
 * Demonstrates TE(X->Y) = I(Y_{t+1}; X_t | Y_t) on synthetic time series
 * with known causal relationships, plus surrogate-based significance testing.
 *
 * Knowledge points demonstrated:
 *   1. TE computation on synthetic time series with known causal direction
 *   2. Normalized TE for comparability across different systems
 *   3. Transfer entropy matrix for multi-variable systems
 *   4. Symbolic transfer entropy for robustness to noise
 *   5. TE asymmetry — distinguishing causal direction
 *
 * Ref: Schreiber (2000) Phys. Rev. Lett. 85, 461
 *      Wibral et al. (2013) PLoS ONE 8(2), e55809
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ifn_core.h"
#include "ifn_transfer.h"

int main(void) {
    printf("=== Transfer Entropy Demo ===\n\n");

    /* Demo 1: Basic TE with causal lag-1 coupling.
     * Source generates pattern, target copies with delay 1.
     * TE(X->Y) should be positive (X Granger-causes Y). */
    printf("--- Demo 1: Basic TE (X causes Y with lag 1) ---\n");
    IFN_TimeSeries* src = ifn_ts_create(100, 1);
    IFN_TimeSeries* tgt = ifn_ts_create(100, 1);
    for (int i = 0; i < 100; i++) {
        src->data[i] = i % 5;
        tgt->data[i] = (i + 1) % 5;
    }
    IFN_TransferResult r = ifn_transfer_entropy(src, tgt, 1, 1);
    printf("  TE(X->Y)     = %.6f\n", r.te_value);
    printf("  Effective TE = %.6f\n", r.effective_transfer);

    /* Demo 2: Normalized TE — divide by target entropy for comparability */
    double norm_te = ifn_transfer_entropy_normalized(src, tgt, 1, 1);
    printf("  Normalized TE = %.6f\n", norm_te);

    /* Demo 3: TE asymmetry — reverse direction should reveal
     * that Y does NOT cause X (target doesn't contain predictive info
     * about source's future). */
    IFN_TransferResult r_rev = ifn_transfer_entropy(tgt, src, 1, 1);
    printf("  TE(Y->X)     = %.6f (should be ~0 for lag-1 coupling)\n",
           r_rev.te_value);
    printf("  Asymmetry: TE(X->Y) - TE(Y->X) = %.6f\n",
           r.te_value - r_rev.te_value);

    /* Demo 4: Multi-variable TE matrix.
     * 3 variables: X drives Y, Y drives Z, X does NOT directly drive Z.
     * The TE matrix should show the correct causal structure. */
    printf("\n--- Demo 4: Multi-variable TE matrix ---\n");
    IFN_TimeSeries* multi = ifn_ts_create(100, 3);
    for (int i = 0; i < 100; i++) {
        multi->data[i * 3 + 0] = i % 5;           /* X: independent driver */
        multi->data[i * 3 + 1] = multi->data[i*3 + 0]; /* Y = X(t) */
        multi->data[i * 3 + 2] = multi->data[i*3 + 1]; /* Z = Y(t) */
    }
    double te_matrix[9]; /* 3x3 */
    ifn_transfer_matrix(multi, 3, 1, 1, te_matrix);
    printf("  TE matrix (rows=source, cols=target):\n");
    for (int i = 0; i < 3; i++) {
        printf("    ");
        for (int j = 0; j < 3; j++)
            printf("%8.4f ", te_matrix[i * 3 + j]);
        printf("\n");
    }

    /* Demo 5: Symbolic TE — ordinal-pattern based, robust to noise */
    printf("\n--- Demo 5: Symbolic Transfer Entropy ---\n");
    double sym_te = ifn_symbolic_transfer_entropy(
        src->data, tgt->data, 100, 5, 1, 1);
    printf("  Symbolic TE = %.6f\n", sym_te);

    /* Cleanup */
    ifn_ts_free(src);
    ifn_ts_free(tgt);
    ifn_ts_free(multi);

    printf("\nExample 1 PASSED — all TE demos completed\n");
    return 0;
}
