/* ==============================================================
 * Example 3: Counterfactual Analysis — PN, PS, PNS & Mediation
 *
 * Counterfactuals answer "what if" questions at the individual
 * level, going beyond population-level causal effects.
 *
 * Three-step counterfactual computation (Pearl, 2000):
 *   1. Abduction: update P(U) given evidence E
 *   2. Action:    apply intervention do(X=x)
 *   3. Prediction: compute probability of consequent Y=y
 *
 * Key quantities:
 *   PN  = P(Y_{X=0}=0 | X=1, Y=1) — probability of necessity
 *   PS  = P(Y_{X=1}=1 | X=0, Y=0) — probability of sufficiency
 *   PNS = P(Y_{X=1}=1, Y_{X=0}=0) — probability of nec & suff
 *
 * References:
 *   Pearl, J. (2000). Causality, Ch. 7.
 *   Balke & Pearl (1994). UAI 10:46-54.
 *   Tian & Pearl (2000). UAI 16:569-577.
 * ============================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dci_core.h"
#include "dci_graph.h"
#include "dci_counterfactual.h"
#include "dci_do.h"

int main(void) {
    printf("=== Example 3: Counterfactual Analysis ===\n\n");

    /* 1. Build SCM for counterfactual analysis */
    DCI_SCM* scm = dci_scm_create("counterfactual_example");
    dci_scm_add_variable(scm, "X", DCI_BINARY, NULL, 0, NULL, true);
    dci_scm_add_variable(scm, "Y", DCI_CONTINUOUS, NULL, 0, NULL, false);
    dci_scm_add_variable(scm, "M", DCI_CONTINUOUS, NULL, 0, NULL, false);
    /* X → M → Y, X → Y (mediation structure) */
    dci_scm_add_edge(scm, 0, 1);  /* X → M */
    dci_scm_add_edge(scm, 0, 2);  /* X → Y */
    dci_scm_add_edge(scm, 1, 2);
    printf("SCM: %s (%d variables)\n", scm->name, scm->n_vars);
    printf("  V0: X (treatment, binary)\n");
    printf("  V1: M (mediator)\n");
    printf("  V2: Y (outcome)\n\n");

    /* 2. Probability of Necessity */
    printf("--- Probability of Necessity (PN) ---\n");
    printf("PN = P(Y_{X=0}=0 | X=1, Y=1)\n");
    printf("  \"Would the patient have recovered without treatment?\"\n");
    double pn = dci_probability_of_necessity(scm, 0, 2, 1000);
    printf("  PN = %.4f\n\n", pn);

    /* 3. Probability of Sufficiency */
    printf("--- Probability of Sufficiency (PS) ---\n");
    printf("PS = P(Y_{X=1}=1 | X=0, Y=0)\n");
    printf("  \"Would the treatment alone cause recovery?\"\n");
    double ps = dci_probability_of_sufficiency(scm, 0, 2, 1000);
    printf("  PS = %.4f\n\n", ps);

    /* 4. Probability of Necessity & Sufficiency */
    printf("--- Probability of Necessity & Sufficiency (PNS) ---\n");
    printf("PNS = P(Y_{X=1}=1, Y_{X=0}=0)\n");
    printf("  \"Is X both necessary and sufficient for Y?\"\n");
    double pns = dci_probability_of_necessity_sufficiency(scm, 0, 2, 1000);
    printf("  PNS = %.4f\n", pns);
    printf("  Bounds: [max(0, P(y1)-P(y0)), min(P(y1), 1-P(y0))]\n\n");

    /* 5. Counterfactual probability via twin network */
    printf("--- Counterfactual Probability ---\n");
    DCI_Counterfactual cf = {
        .antecedent_var = 0,
        .antecedent_value = 0.0,
        .consequent_var = 2,
        .consequent_value = 1.0,
    };
    snprintf(cf.description, 255, "P(Y_{X=0} = 1)");
    DCI_CounterfactualResult r = dci_counterfactual_probability(
        scm, &cf, NULL, NULL, 0, 500);
    printf("Query: %s\n", cf.description);
    printf("  Result: %.4f (identifiable: %s, samples: %d)\n\n",
        r.probability, r.is_identifiable ? "YES" : "NO", r.n_samples);

    /* 6. Counterfactual with evidence (3-step process) */
    printf("--- Counterfactual with Evidence ---\n");
    printf("P(Y_{X=1}=5 | evidence: X=0, Y=3)\n");
    double evidence_vals[] = {0.0, 3.0};
    int evidence_vars[] = {0, 2};
    DCI_Counterfactual cf2 = {.antecedent_var=0, .antecedent_value=1.0,
        .consequent_var=2, .consequent_value=5.0};
    DCI_CounterfactualResult r2 = dci_counterfactual_probability(
        scm, &cf2, evidence_vals, evidence_vars, 2, 500);
    printf("  Result: %.4f\n\n", r2.probability);

    /* 7. Individual-level counterfactual */
    printf("--- Individual-Level Counterfactual ---\n");
    printf("Given a specific patient (X=1, M=0.7, Y=4.2),\n");
    printf("what would Y be had X been 0?\n");
    double obs_full[] = {1.0, 0.7, 4.2};
    int obs_vars[] = {0, 1, 2};
    double ind_cf = dci_individual_counterfactual(scm,
        obs_full, obs_vars, 3, 0, 0.0, 2);
    printf("  Y_{X=0}(u) = %.4f\n\n", ind_cf);

    /* 8. Causal mediation via counterfactuals */
    printf("--- Causal Mediation ---\n");
    printf("X → M → Y, X → Y\n");
    double mediation = dci_causal_mediation(scm, 0, 1, 2, 1000);
    printf("  Mediation effect: %.4f\n\n", mediation);

    /* 9. Path-specific effect (X→M→Y path) */
    printf("--- Path-Specific Counterfactual Effect ---\n");
    int path_nodes[] = {0, 1, 2};
    double pse = dci_path_specific_effect(scm, 0, 2, path_nodes, 3, 1000);
    printf("  Effect via X→M→Y: %.4f\n\n", pse);

    /* 10. Attributable fraction */
    printf("--- Attributable Fraction ---\n");
    double af = dci_attributable_fraction(scm, 0, 2, 1000);
    printf("  AF (population attributable fraction): %.4f\n\n", af);

    dci_scm_free(scm);
    printf("=== Example 3 Complete ===\n");
    return 0;
}
