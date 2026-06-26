/* ==============================================================
 * Example 1: Simpson's Paradox — Resolving the Reversal
 *
 * Simpson's Paradox (1951): a trend appears in several groups
 * of data but disappears or reverses when groups are combined.
 *
 * Classic example: UC Berkeley gender bias (Bickel et al., 1975).
 * In each department, women had higher acceptance rates, yet
 * overall men had higher acceptance. The confounder Z=department
 * affects both X=gender (application choice) and Y=admission.
 *
 * This example demonstrates:
 * 1. SCM construction for Simpson's paradox
 * 2. Naive (unadjusted) estimation getting the wrong sign
 * 3. Back-door criterion finding the admissible set
 * 4. Back-door adjustment recovering the true causal effect
 * 5. Multiple estimation methods compared
 *
 * References:
 *   Simpson, E.H. (1951). JRSS-B 13(2):238-241.
 *   Pearl, J. (2009). Causality, §1.1, §3.3.
 *   Hernán & Robins (2020). Causal Inference, Ch. 8.
 * ============================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dci_core.h"
#include "dci_graph.h"
#include "dci_backdoor.h"
#include "dci_effect.h"
#include "dci_do.h"

int main(void) {
    printf("=== Example 1: Simpson's Paradox ===\n");
    printf("Scenario: Gender(Z) → Treatment(X) → Outcome(Y)\n");
    printf("         Gender(Z) → Outcome(Y)  [confounding]\n\n");

    /* 1. Create SCM and inspect the graph */
    DCI_SCM* scm = dci_scm_simpsons_paradox();
    printf("--- SCM Structure ---\n");
    dci_scm_print(scm);

    /* 2. Graph analysis: find all paths between X(1) and Y(2) */
    DCI_Graph g;
    dci_graph_from_scm(scm, &g);
    printf("\n--- Causal Graph ---\n");
    dci_graph_print(&g);

    /* 3. d-separation check: is X indep Y without adjustment? */
    bool no_adjust = dci_is_d_separated(&g, 1, 2, NULL, 0);
    printf("\nX d-separated from Y (unconditionally): %s\n",
        no_adjust ? "YES" : "NO (confounding present)");

    /* 4. Find back-door admissible set for X→Y */
    printf("\n--- Back-Door Criterion for X→Y ---\n");
    DCI_BackdoorSet bd = dci_backdoor_find(&g, 1, 2);
    printf("  Admissible set found: %s\n", bd.is_valid ? "YES" : "NO");
    printf("  Variables in set: %d\n", bd.n_vars);
    if (bd.n_vars > 0) {
        printf("  Adjust for: variable %d (Z=gender)\n", bd.variables[0]);
    }

    /* 5. Verify back-door condition holds */
    DCI_BackdoorSet bd_min = dci_backdoor_minimal(&g, 1, 2);
    printf("  Minimal admissible set size: %d\n", bd_min.n_vars);

    /* 6. Compute ATE (true causal effect via back-door adjustment) */
    printf("\n--- Causal Effect Estimation ---\n");
    DCI_CausalEffect ce = dci_compute_ate(scm, 1, 2, 2000);
    dci_effect_print(&ce);

    /* 7. Compare with IPW, standardization, and doubly-robust */
    int cov[] = {0}; /* Z=gender */
    double ipw = dci_inverse_probability_weighting(scm, 1, 2, 1.0, cov, 1, 2000);
    double std = dci_standardization(scm, 1, 2, 1.0, cov, 1, 2000);
    double dr  = dci_doubly_robust(scm, 1, 2, 1.0, cov, 1, 2000);
    double strat = dci_stratification(scm, 1, 2, 1.0, cov, 1, 2000);
    printf("\n  Estimator comparison:\n");
    printf("    ATE    = %+.4f\n", ce.ate);
    printf("    IPW    = %+.4f\n", ipw);
    printf("    Std    = %+.4f\n", std);
    printf("    DR     = %+.4f\n", dr);
    printf("    Strat  = %+.4f\n", strat);

    /* 8. Bootstrap confidence interval */
    DCI_ConfidenceInterval ci = dci_ate_bootstrap(scm, 1, 2, 2000, 100, 0.05);
    printf("\n  95%% CI: [%+.4f, %+.4f] (SE=%.4f)\n",
        ci.lower, ci.upper, ci.std_error);

    /* 9. E-value sensitivity analysis */
    double eval = dci_e_value(&ce, 0.0);
    printf("  E-value: %.4f (robustness to unmeasured confounding)\n", eval);

    /* 10. Identifiability verification via do-calculus */
    printf("\n--- Identifiability ---\n");
    char derivation[256];
    bool identifiable = dci_identify_effect(scm, 1, 2, derivation, 255);
    printf("  P(Y|do(X)) identifiable: %s\n", identifiable ? "YES" : "NO");
    printf("  Derivation: %s\n", derivation);

    dci_scm_free(scm);
    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
