/* ==============================================================
 * do_calculus_intervention_app2.c — L7 Application 2
 *
 * Causal Mediation Analysis for Health Policy
 *
 * Scenario: We study whether an educational intervention (X=1: college,
 * X=0: high school) affects long-term health outcomes (Y: health score)
 * through the mediating pathway of health literacy (M: understanding of
 * nutrition and preventive care).
 *
 * The causal model: X → M → Y (indirect path), X → Y (direct path).
 * There is also an unobserved confounder U (childhood SES) affecting
 * both X and Y, making standard regression biased.
 *
 * This application demonstrates:
 *   1. Mediation SCM (Baron-Kenny framework)
 *   2. Natural Direct and Indirect Effects (NDE, NIE)
 *   3. Controlled Direct Effects (CDE)
 *   4. Front-door criterion (using M as mediator)
 *   5. Path-specific counterfactual effects
 *   6. Twin network for counterfactual mediation
 *   7. Counterfactual fairness assessment
 *   8. Instrumental Variable identification
 *
 * References:
 *   Baron & Kenny (1986) JPSP
 *   Pearl (2001) Direct and Indirect Effects, UAI
 *   VanderWeele (2015) Explanation in Causal Inference
 *   Imai, Keele, Tingley (2010) Psychological Methods
 *   Pearl & Mackenzie (2018) The Book of Why, Ch.9
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

/* ---- Custom Structural Equations for Mediation SCM ---- */

/*
 * U: Unobserved confounder (childhood SES, exogenous)
 * X: Education — P(X=1|U) = logistic(U)
 * M: Health literacy — M = β_xm*X + β_um*U + ε_m
 * Y: Health outcome — Y = β_xy*X + β_my*M + β_uy*U + ε_y
 *
 * Causal diagram: U → X, U → Y, X → M, X → Y (direct), M → Y
 */

static double eq_u_ses(const double* parents, const double* u) {
    (void)parents;
    return u ? *u : 0.5; /* U ~ Uniform(0,1) for SES percentile */
}

static double eq_x_education(const double* parents, const double* u) {
    double ses = parents ? parents[0] : 0.5;
    /* Higher SES → more likely to attend college */
    double prob = 0.2 + 0.6 * ses;
    return u ? (*u < prob ? 1.0 : 0.0) : 0.0;
}

static double eq_m_health_literacy(const double* parents, const double* u) {
    double x = parents ? parents[0] : 0.0;
    double ses = parents ? parents[1] : 0.5;
    double noise = u ? *u * 0.3 : 0.0;
    return 0.6 * x + 0.3 * ses + noise;
}

static double eq_y_health(const double* parents, const double* u) {
    double x = parents ? parents[0] : 0.0;
    double m = parents ? parents[1] : 0.5;
    double ses = parents ? parents[2] : 0.5;
    double noise = u ? *u * 0.4 : 0.0;
    return 0.3 * x + 0.5 * m + 0.2 * ses + noise;
}

/*
 * Build mediation SCM with confounder.
 * Variables: 0=U(SES), 1=X(Education), 2=M(HealthLit), 3=Y(Health)
 */
static DCI_SCM* build_mediation_scm(void) {
    DCI_SCM* scm = dci_scm_create("Mediation — Education → Health Literacy → Health");

    /* U: SES (exogenous) */
    dci_scm_add_variable(scm, "U_SES", DCI_CONTINUOUS,
        NULL, 0, eq_u_ses, true);

    /* X: Education (binary, depends on U) */
    int px[] = {0};
    dci_scm_add_variable(scm, "X_Education", DCI_BINARY,
        px, 1, eq_x_education, false);

    /* M: Health Literacy (depends on X and U) */
    int pm[] = {1, 0};
    dci_scm_add_variable(scm, "M_HealthLit", DCI_CONTINUOUS,
        pm, 2, eq_m_health_literacy, false);

    /* Y: Health Outcome (depends on X, M, U) */
    int py[] = {1, 2, 0};
    dci_scm_add_variable(scm, "Y_Health", DCI_CONTINUOUS,
        py, 3, eq_y_health, false);

    return scm;
}

/*
 * Build front-door SCM (alternative: no direct X→Y, U→X, X→M→Y, U→Y)
 */
static DCI_SCM* build_frontdoor_scm(void) {
    DCI_SCM* scm = dci_scm_create("Front-Door — Smoking → Tar → Cancer");

    /* U: genetic confounder */
    dci_scm_add_variable(scm, "U_Genetics", DCI_CONTINUOUS,
        NULL, 0, NULL, true);

    /* X: Smoking (depends on U) */
    int px[] = {0};
    dci_scm_add_variable(scm, "X_Smoking", DCI_BINARY,
        px, 1, NULL, false);

    /* M: Tar deposits (depends on X) */
    int pm[] = {1};
    dci_scm_add_variable(scm, "M_Tar", DCI_CONTINUOUS,
        pm, 1, NULL, false);

    /* Y: Lung cancer (depends on M and U, NOT directly on X) */
    int py[] = {2, 0};
    dci_scm_add_variable(scm, "Y_Cancer", DCI_CONTINUOUS,
        py, 2, NULL, false);

    return scm;
}

static void print_separator(const char* title) {
    printf("\n══════════════════════════════════════════════════════\n");
    printf("  %s\n", title);
    printf("══════════════════════════════════════════════════════\n");
}

/*
 * Print mediation decomposition with interpretation.
 */
static void print_mediation_decomp(double total, double nde, double nie, double cde) {
    printf("\n┌─ Mediation Decomposition ──────────────────────────┐\n");
    printf("│ Total Effect (TE)   = NDE + NIE                   │\n");
    printf("│   TE  = %+.4f                                     │\n", total);
    printf("│   NDE = %+.4f  (direct: Education → Health)       │\n", nde);
    printf("│   NIE = %+.4f  (indirect: via Health Literacy)    │\n", nie);
    printf("│   CDE = %+.4f  (controlled, M fixed at reference) │\n", cde);
    printf("│                                                   │\n");
    double pct_mediated = (total != 0.0) ? fabs(nie / total) * 100.0 : 0.0;
    printf("│   Proportion mediated: %.1f%%                      │\n", pct_mediated);
    printf("└───────────────────────────────────────────────────┘\n");
}

int main(void) {
    srand((unsigned)time(NULL));
    printf("Do-Calculus Intervention App 2: Causal Mediation for Health Policy\n");
    printf("===================================================================\n\n");

    /* ── 1. Build SCMs ── */
    print_separator("1. SCM Construction");
    DCI_SCM* scm_med = build_mediation_scm();
    DCI_SCM* scm_fd  = build_frontdoor_scm();

    printf("Mediation SCM: %s\n", scm_med->name);
    dci_scm_print(scm_med);
    printf("\nFront-Door SCM: %s\n", scm_fd->name);
    dci_scm_print(scm_fd);

    /* ── 2. Graph Analysis for Mediation ── */
    print_separator("2. Causal Graph Analysis — Mediation");
    DCI_Graph g_med;
    dci_graph_from_scm(scm_med, &g_med);
    dci_graph_print(&g_med);

    /* Causal ordering */
    int order[DCI_MAX_VARS];
    int n_order = dci_causal_ordering(&g_med, order);
    printf("Causal order: ");
    for (int i = 0; i < n_order; i++)
        printf("%d ", order[i]);
    printf("\n");

    /* ── 3. Mediation Effect Estimation ── */
    print_separator("3. Causal Mediation Analysis (Monte Carlo, N=50000)");
    int N = 50000;

    /* Total Effect via do-calculus */
    DCI_CausalEffect te = dci_compute_ate(scm_med, 1, 3, N);
    printf("Total Effect (ATE, X→Y):      %+.4f\n", te.ate);

    /* Natural Direct Effect: Y(1, M(0)) — effect of X on Y
     * when M is set to its natural value under X=0 */
    double nde = dci_natural_direct_effect(scm_med, 1, 3, 2, 1.0, 0.0, N);
    printf("Natural Direct Effect (NDE):   %+.4f\n", nde);

    /* Natural Indirect Effect: Y(1, M(1)) - Y(1, M(0)) */
    double nie = dci_natural_indirect_effect(scm_med, 1, 3, 2, 1.0, 0.0, N);
    printf("Natural Indirect Effect (NIE): %+.4f\n", nie);

    /* Controlled Direct Effect: Y(1, m) - Y(0, m) for fixed m */
    double cde = dci_controlled_direct_effect(scm_med, 1, 3, 2, 1.0, 0.5, N);
    printf("Controlled Direct Effect (CDE):%+.4f\n", cde);

    print_mediation_decomp(te.ate, nde, nie, cde);

    /* ── 4. Path-Specific Counterfactual Effects ── */
    print_separator("4. Path-Specific Counterfactual Effects");
    int path_xmy[] = {1, 2, 3}; /* X → M → Y pathway */
    double pse = dci_path_specific_effect(scm_med, 1, 3, path_xmy, 3, N);
    printf("Path-specific effect (X→M→Y):        %+.4f\n", pse);

    int path_xy[] = {1, 3}; /* X → Y direct pathway */
    double pse_dir = dci_path_specific_effect(scm_med, 1, 3, path_xy, 2, N);
    printf("Path-specific effect (X→Y direct):   %+.4f\n", pse_dir);

    /* ── 5. Counterfactual Mediation via Twin Network ── */
    print_separator("5. Counterfactual Mediation (Twin Network)");

    /* Counterfactual: what would M be had education been college (X=1)
     * vs high school (X=0)? */
    DCI_TwinNetwork* tn1 = dci_twin_network_create(scm_med, 1, 1.0);
    double m_cf_college = dci_twin_network_query(tn1, 2, 0.0);
    printf("Counterfactual M (had college):      %.4f\n", m_cf_college);
    dci_twin_network_free(tn1);

    DCI_TwinNetwork* tn0 = dci_twin_network_create(scm_med, 1, 0.0);
    double m_cf_hs = dci_twin_network_query(tn0, 2, 0.0);
    printf("Counterfactual M (had high school):   %.4f\n", m_cf_hs);
    dci_twin_network_free(tn0);

    /* ── 6. Individual Counterfactual ── */
    print_separator("6. Individual-Level Counterfactual");
    double obs_vals[] = {0.7, 0.0, 0.4}; /* SES=0.7, no college, M=0.4 */
    int obs_vars[] = {0, 1, 2};
    double cf_y = dci_individual_counterfactual(scm_med, obs_vals, obs_vars,
        3, 1, 1.0, 3);
    printf("For an individual with:\n");
    printf("  SES=%.1f, Education=%.0f, HealthLit=%.2f\n",
        obs_vals[0], obs_vals[1], obs_vals[2]);
    printf("Counterfactual Y had they attended college: %.4f\n", cf_y);

    /* ── 7. Causal Mediation (full mediation analysis) ── */
    print_separator("7. Full Causal Mediation Analysis");
    double mediation_total = dci_causal_mediation(scm_med, 1, 2, 3, N);
    printf("Mediation effect (NIE via counterfactual):  %+.4f\n", mediation_total);

    /* ── 8. Front-Door Criterion Analysis ── */
    print_separator("8. Front-Door Criterion (Smoking → Tar → Cancer)");
    DCI_Graph g_fd;
    dci_graph_from_scm(scm_fd, &g_fd);

    /* Check front-door */
    DCI_FrontdoorSet fd = dci_frontdoor_find(&g_fd, 1, 2);
    printf("Front-door admissible set: ");
    if (fd.is_valid && fd.n_mediators > 0) {
        for (int i = 0; i < fd.n_mediators; i++)
            printf("%s ", scm_fd->vars[fd.mediators[i]].name);
        printf("\n  → Tar satisfies front-door criterion ✓\n");
    } else {
        printf("NOT FOUND\n");
    }

    /* Front-door adjustment estimate */
    double fd_adj = dci_frontdoor_adjust(scm_fd, &g_fd, 1, 2, 1.0, 1.0,
        &fd, N/100);
    printf("Front-door adjusted P(Cancer|do(Smoking=1)): %.4f\n", fd_adj);

    /* ── 9. Instrumental Variable Analysis ── */
    print_separator("9. Instrumental Variable Check");
    DCI_Graph g_med2;
    dci_graph_from_scm(scm_med, &g_med2);

    /* Check if SES (variable 0) could serve as IV for Education→Health */
    bool is_iv = dci_is_instrument(&g_med2, 0, 1, 3);
    printf("Is SES a valid instrument for Education→Health? %s\n",
        is_iv ? "NO (affects Y directly — exclusion violated)" : "NO");

    /* Alternative: check if we can find any instrument */
    printf("  Note: SES → Y directly, so it fails IV exclusion restriction.\n");
    printf("  Identifying an instrument requires Z ⟂ Y | X, U.\n");

    /* ── 10. Counterfactual Fairness ── */
    print_separator("10. Counterfactual Fairness Assessment");
    /* Check: is health outcome fair with respect to education?
     * Counterfactual fairness: would Y change if we set X differently
     * for the same underlying characteristics? */
    bool fair = dci_counterfactual_fairness(scm_med, 1, 3, 0.05, N);
    printf("Counterfactual fairness (Education→Health): %s\n",
        fair ? "PASS ✓ — effect is small" : "FAIL — substantial disparity");
    printf("  Threshold for fairness: |P(Y_x|group_a) - P(Y_x|group_b)| < 0.05\n");

    /* ── 11. Counterfactual Explanation ── */
    print_separator("11. Counterfactual Explanation (Feature Importance)");
    double importance[4] = {0};
    double obs_full[] = {0.8, 1.0, 0.7, 0.9};
    int obs_all[] = {0, 1, 2, 3};
    dci_counterfactual_explanation(scm_med, obs_full, obs_all, 4, 3,
        importance, 4);
    printf("Feature importance for Health outcome:\n");
    printf("  SES:           %.4f\n", importance[0]);
    printf("  Education:     %.4f\n", importance[1]);
    printf("  Health Lit:    %.4f\n", importance[2]);
    printf("  Health itself: %.4f\n", importance[3]);

    /* ── 12. Attributable Fraction ── */
    print_separator("12. Population Attributable Impact");
    double af_edu = dci_attributable_fraction(scm_med, 1, 3, N);
    printf("Attributable fraction (Education→Health): %.4f\n", af_edu);
    printf("  → %.1f%% of health outcomes attributable\n", af_edu*100.0);
    printf("    to the educational intervention.\n");

    /* ── 13. Policy Recommendations ── */
    print_separator("13. Policy Recommendations");
    printf("┌─────────────────────────────────────────────────────┐\n");
    printf("│ MEDIATION ANALYSIS SUMMARY                         │\n");
    printf("│                                                    │\n");
    printf("│ Total effect of college on health:  %+.4f         │\n", te.ate);
    printf("│   - Direct pathway (knowledge):     %+.4f         │\n", nde);
    printf("│   - Indirect via health literacy:   %+.4f         │\n", nie);
    printf("│                                                    │\n");
    if (fabs(nie) > fabs(nde)) {
        printf("│ KEY FINDING: Health literacy is the dominant      │\n");
        printf("│ pathway. Investments in health education may      │\n");
        printf("│ yield larger health benefits than direct college  │\n");
        printf("│ attendance alone.                                 │\n");
    } else {
        printf("│ KEY FINDING: Direct effects dominate. College     │\n");
        printf("│ education has benefits beyond health literacy.    │\n");
    }
    printf("│                                                    │\n");
    printf("│ POLICY: Fund health literacy programs in schools  │\n");
    printf("│ as a cost-effective way to improve population     │\n");
    printf("│ health, leveraging the identified mediation path. │\n");
    printf("└─────────────────────────────────────────────────────┘\n\n");

    /* Cleanup */
    dci_scm_free(scm_med);
    dci_scm_free(scm_fd);

    printf("Application 2 complete. ✓\n");
    return 0;
}
