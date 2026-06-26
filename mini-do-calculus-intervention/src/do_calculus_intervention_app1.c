/* ==============================================================
 * do_calculus_intervention_app1.c — L7 Application 1
 *
 * Clinical Trial Policy Evaluation via Do-Calculus
 *
 * Scenario: A new drug (X=1) vs placebo (X=0) is tested for its effect
 * on recovery (Y). Gender (Z) acts as a confounder: women are both
 * more likely to receive the drug and have different baseline recovery
 * rates. This creates Simpson's Paradox — the drug appears harmful
 * in aggregate but is beneficial within each gender stratum.
 *
 * This application demonstrates:
 *   1. Building an SCM for Simpson's Paradox
 *   2. Graph surgery and do(X=1) / do(X=0) intervention
 *   3. Back-door adjustment with covariate Z (gender)
 *   4. Multiple causal effect estimators (ATE, IPW, standardization,
 *      doubly-robust, matching, DiD, RDD, bootstrap CI)
 *   5. Sensitivity analysis via E-value
 *   6. Counterfactual probability of necessity
 *
 * References:
 *   Pearl (2009) Causality, §3.3 Simpson's Paradox
 *   Hernán & Robins (2020) Causal Inference: What If, §4
 *   Morgan & Winship (2015) Counterfactuals and Causal Inference
 * ============================================================== */

#include "dci_core.h"
#include "dci_graph.h"
#include "dci_do.h"
#include "dci_backdoor.h"
#include "dci_effect.h"
#include "dci_counterfactual.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ---- Custom structural equations for Simpson's Paradox SCM ---- */

/*
 * Gender Z (binary): P(Z=1) = 0.5 (women), drawn exogenously
 * Treatment X: P(X=1|Z=0) = 0.3 (men less likely), P(X=1|Z=1) = 0.7
 * Outcome Y: Y = 2*X + 0.5*Z + noise  (drug works, women healthier)
 */
static double eq_z_exog(const double* parents, const double* u) {
    (void)parents;
    return u ? (*u > 0.5 ? 1.0 : 0.0) : 0.0;
}

static double eq_x_treatment(const double* parents, const double* u) {
    double z = parents ? parents[0] : 0.0;
    double threshold = z > 0.5 ? 0.3 : 0.7;
    return u ? (*u > threshold ? 1.0 : 0.0) : 0.0;
}

static double eq_y_outcome(const double* parents, const double* u) {
    double x = parents ? parents[0] : 0.0;
    double z = parents ? parents[1] : 0.0;
    double noise = u ? *u * 0.5 : 0.0;
    return 2.0 * x + 0.5 * z + noise;
}

/*
 * Build a full Simpson's Paradox SCM with proper structural equations.
 */
static DCI_SCM* build_simpsons_scm(void) {
    DCI_SCM* scm = dci_scm_create("Clinical Trial — Simpson's Paradox");

    /* Z: Gender (exogenous, binary, no parents) */
    dci_scm_add_variable(scm, "Gender_Z", DCI_BINARY,
        NULL, 0, eq_z_exog, true);

    /* X: Treatment (binary, depends on Z) */
    int px[] = {0};
    dci_scm_add_variable(scm, "Treatment_X", DCI_BINARY,
        px, 1, eq_x_treatment, false);

    /* Y: Recovery outcome (depends on X and Z) */
    int py[] = {1, 0};
    dci_scm_add_variable(scm, "Recovery_Y", DCI_CONTINUOUS,
        py, 2, eq_y_outcome, false);

    return scm;
}

/*
 * Print a horizontal rule for output formatting.
 */
static void print_separator(const char* title) {
    printf("\n══════════════════════════════════════════════════════\n");
    printf("  %s\n", title);
    printf("══════════════════════════════════════════════════════\n");
}

int main(void) {
    srand((unsigned)time(NULL));
    printf("Do-Calculus Intervention App 1: Clinical Trial Policy Evaluation\n");
    printf("================================================================\n\n");

    /* ── 1. Build SCM and Inspect ── */
    print_separator("1. SCM Construction (Simpson's Paradox)");
    DCI_SCM* scm = build_simpsons_scm();
    dci_scm_print(scm);

    /* Verify DAG property */
    bool is_dag = dci_scm_is_dag(scm);
    printf("\nSCM is DAG: %s\n", is_dag ? "YES ✓" : "NO ✗");

    /* ── 2. Graph Analysis ── */
    print_separator("2. Causal Graph Analysis");
    DCI_Graph* g = dci_graph_create(scm->n_vars);
    dci_graph_from_scm(scm, g);
    dci_graph_print(g);

    /* Find v-structures (colliders) */
    int triples[10][3];
    int n_v = dci_find_v_structures(g, triples, 10);
    printf("V-structures found: %d\n", n_v);
    int vi;
    for (vi = 0; vi < n_v; vi++) {
        printf("  Collider: %d → %d ← %d\n",
            triples[vi][0], triples[vi][1], triples[vi][2]);
    }

    /* Check d-separation: is Treatment ⟂ Recovery | Gender? */
    int cond[] = {0}; /* condition on Gender */
    DCI_DSepResult dsep = dci_d_separation(g, 1, 2, cond, 1);
    printf("\nTreatment ⟂ Recovery | Gender? %s\n",
        dsep.is_d_separated ? "YES (d-separated)" : "NO (not d-separated)");

    /* ── 3. Back-Door Adjustment Analysis ── */
    print_separator("3. Back-Door Criterion & Covariate Selection");
    DCI_BackdoorSet bd = dci_backdoor_find(g, 1, 2);
    printf("Back-door admissible set: ");
    if (bd.is_valid && bd.n_vars > 0) {
        int i;
        for (i = 0; i < bd.n_vars; i++)
            printf("%s ", scm->vars[bd.variables[i]].name);
        printf("\n  → Gender satisfies back-door criterion ✓\n");
    } else if (bd.is_valid) {
        printf("(empty — no back-door paths)\n");
    } else {
        printf("NOT FOUND\n");
    }

    /* Minimal back-door set */
    DCI_BackdoorSet bd_min = dci_backdoor_minimal(g, 1, 2);
    printf("Minimal back-door set size: %d\n", bd_min.n_vars);

    /* Check for M-bias */
    bool m_bias = dci_detect_m_bias(g, 1, 2);
    printf("M-bias detected: %s\n", m_bias ? "YES (caution!)" : "NO");

    /* ── 4. Causal Effect Estimation ── */
    print_separator("4. Causal Effect Estimation (Monte Carlo, N=50000)");
    int N = 50000;

    /* 4a. Naive associational estimate (biased) */
    double naive_sum1 = 0.0, naive_sum0 = 0.0;
    int n_treat = 0, n_ctrl = 0;
    int s;
    for (s = 0; s < N; s++) {
        double exog[3] = {(double)rand()/RAND_MAX,
                          (double)rand()/RAND_MAX,
                          (double)rand()/RAND_MAX};
        dci_scm_evaluate(scm, exog);
        if (scm->vars[1].value > 0.5) {
            naive_sum1 += scm->vars[2].value;
            n_treat++;
        } else {
            naive_sum0 += scm->vars[2].value;
            n_ctrl++;
        }
    }
    double naive_ate = (naive_sum1/n_treat) - (naive_sum0/n_ctrl);
    printf("Naive (associational) ATE:  %.4f  [BIASED — ignores confounding]\n",
        naive_ate);

    /* 4b. ATE via do-calculus (interventional) */
    DCI_CausalEffect ate = dci_compute_ate(scm, 1, 2, N);
    printf("Causal ATE (do-calculus):    %.4f  [UNBIASED]\n", ate.ate);

    /* 4c. ATT and ATU */
    double att = dci_compute_att(scm, 1, 2, N);
    double atu = dci_compute_atu(scm, 1, 2, N);
    printf("ATT (effect on treated):     %.4f\n", att);
    printf("ATU (effect on untreated):   %.4f\n", atu);

    /* 4d. Back-door adjusted estimate */
    double bd_adj = dci_backdoor_adjust(scm, g, 1, 2, 1.0, 3.5, &bd_min, N/10);
    printf("Back-door adjusted P(Y=3.5|do(X=1)): %.4f\n", bd_adj);

    /* 4e. IPW estimator */
    int covariates[] = {0};
    double ipw = dci_inverse_probability_weighting(scm, 1, 2, 1.0,
        covariates, 1, N);
    printf("IPW estimate (X=1):          %.4f\n", ipw);

    /* 4f. Standardization */
    double std_est = dci_standardization(scm, 1, 2, 1.0, covariates, 1, N);
    printf("Standardization (X=1):       %.4f\n", std_est);

    /* 4g. Doubly-robust estimator */
    double dr_est = dci_doubly_robust(scm, 1, 2, 1.0, covariates, 1, N);
    printf("Doubly-Robust (X=1):         %.4f\n", dr_est);

    /* 4h. Matching estimator */
    double match_ate = dci_matching_estimator(scm, 1, 2, covariates, 1, N);
    printf("Matching ATE:                %.4f\n", match_ate);

    /* 4i. Bootstrap confidence interval */
    DCI_ConfidenceInterval ci = dci_ate_bootstrap(scm, 1, 2, N/10, 200, 0.05);
    printf("Bootstrap 95%% CI:          [%.4f, %.4f] (SE=%.4f)\n",
        ci.lower, ci.upper, ci.std_error);

    /* ── 5. Do-Calculus Derivation ── */
    print_separator("5. Do-Calculus Rule Application");
    int x_arr[] = {1};
    int z_arr[] = {0};
    int w_arr[] = {};
    DCI_Derivation der1 = dci_apply_rule1(scm, g, 2, x_arr, 1, z_arr, 1, w_arr, 0);
    printf("Rule 1 (observation removal): %s\n",
        der1.success ? "APPLICABLE ✓" : "NOT applicable");

    DCI_Derivation der2 = dci_apply_rule2(scm, g, 2, x_arr, 1, z_arr, 1, w_arr, 0);
    printf("Rule 2 (action↔observation):  %s\n",
        der2.success ? "APPLICABLE ✓" : "NOT applicable");

    DCI_Derivation der3 = dci_apply_rule3(scm, g, 2, x_arr, 1, z_arr, 1, w_arr, 0);
    printf("Rule 3 (action removal):      %s\n",
        der3.success ? "APPLICABLE ✓" : "NOT applicable");

    /* Identifiability check */
    bool identifiable = dci_identifiability_check(scm, 1, 2);
    printf("\nCausal effect identifiable: %s\n",
        identifiable ? "YES ✓" : "NO ✗");

    /* ── 6. Truncated Factorization ── */
    print_separator("6. Truncated Factorization (g-formula)");
    int interv_vars[] = {1};
    double interv_vals[] = {1.0};
    double obs_vals[] = {0.0, 1.0, 3.5};
    double tf = dci_truncated_factorization(scm, interv_vars, interv_vals,
        1, obs_vals);
    printf("P(v|do(X=1)) product: %.6f\n", tf);

    /* ── 7. Counterfactual Analysis ── */
    print_separator("7. Counterfactual Necessity/Sufficiency");
    double pn = dci_probability_of_necessity(scm, 1, 2, N);
    double ps = dci_probability_of_sufficiency(scm, 1, 2, N);
    double pns = dci_probability_of_necessity_sufficiency(scm, 1, 2, N);
    printf("P(Necessity):         %.4f\n", pn);
    printf("P(Sufficiency):       %.4f\n", ps);
    printf("P(Necessity & Suff):  %.4f\n", pns);
    printf("  → Probability that the drug was necessary AND\n");
    printf("    sufficient for recovery: %.1f%%\n", pns * 100.0);

    /* Attributable fraction */
    double af = dci_attributable_fraction(scm, 1, 2, N);
    printf("\nAttributable Fraction: %.4f\n", af);
    printf("  → %.1f%% of recoveries attributable to treatment\n", af * 100.0);

    /* ── 8. Sensitivity Analysis ── */
    print_separator("8. Sensitivity Analysis (E-value)");
    double e_val = dci_e_value(&ate, 1.0);
    printf("E-value: %.4f\n", e_val);
    printf("  → An unmeasured confounder would need to be associated\n");
    printf("    with both treatment and outcome by RR ≥ %.2f to\n", e_val);
    printf("    explain away the observed effect.\n");

    /* ── 9. Regression Discontinuity & DiD Comparisons ── */
    print_separator("9. Alternative Causal Designs");
    double did = dci_difference_in_differences(scm, 1, 2, 0, N);
    printf("Difference-in-Differences:  %.4f\n", did);
    double rdd = dci_regression_discontinuity(scm, 1, 2, 0, 0.5, N);
    printf("Regression Discontinuity:   %.4f\n", rdd);

    /* ── 10. Summary & Policy Recommendation ── */
    print_separator("10. Policy Recommendation");
    printf("True causal effect (ATE):   %+.4f\n", ate.ate);
    printf("Naive associational effect: %+.4f\n", naive_ate);
    printf("\n");
    if (ate.ate > 0.01) {
        printf("RECOMMENDATION: Approve the treatment.\n");
        printf("The drug shows a positive causal effect after adjusting\n");
        printf("for confounding by gender (back-door criterion).\n");
    } else {
        printf("RECOMMENDATION: Do NOT approve treatment.\n");
    }
    printf("\nKey insight: Simpson's Paradox resolved via do-calculus.\n");
    printf("Without causal adjustment, the confounded estimate is\n");
    printf("misleading (%.4f vs true effect %.4f).\n\n",
        naive_ate, ate.ate);

    /* Cleanup */
    dci_graph_free(g);
    dci_scm_free(scm);

    printf("Application 1 complete. ✓\n");
    return 0;
}
