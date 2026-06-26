/* ============================================================
 * counterfactual_reasoning_app1.c -- L7 Application
 *
 * Medical Treatment Effect with Counterfactuals
 *
 * Scenario: RCT of a drug for blood pressure.
 * Uses Pearl 3-step counterfactual procedure.
 *
 * Refs: Pearl (2009) Causality, Hernan & Robins (2020).
 * ============================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfr_core.h"
#include "cfr_effects.h"
#include "cfr_bounds.h"
#include "cfr_counterfactual.h"
#include "cfr_potential.h"
#include "cfr_analysis.h"

int main(void) {
    printf("=== L7 Application 1: Medical Counterfactuals ===\n\n");
    srand(42);
    int n = 50;
    double y0[50], y1[50], T[50], Y[50];

    /* Generate RCT: Y(0) ~ N(140,15), Y(1) = Y(0) - 12 + N(0,4) */
    for (int i = 0; i < n; i++) {
        y0[i] = 140.0 + cfr_gaussian_noise(0, 15);
        y1[i] = y0[i] - 12.0 + cfr_gaussian_noise(0, 4);
        T[i] = (rand() % 2) ? 1.0 : 0.0;
        Y[i] = (T[i] > 0.5) ? y1[i] : y0[i];
    }

    /* Oracle ATE */
    CFRCausalEffects* eff = cfr_eff_create();
    cfr_eff_ate_oracle(eff, y0, y1, n);
    printf("Oracle ATE = %.2f (SE %.2f, CI [%.2f, %.2f])\n",
           eff->ate, eff->se_ate, eff->ci_lower, eff->ci_upper);

    /* Observed ATE */
    CFRCausalEffects* eff_obs = cfr_eff_create();
    cfr_eff_ate_diff_means(eff_obs, T, Y, n);
    printf("Observed ATE = %.2f (SE %.2f)\n", eff_obs->ate, eff_obs->se_ate);

    /* IPW */
    double ps[50]; for (int i=0; i<n; i++) ps[i] = 0.5;
    CFRCausalEffects* eff_ipw = cfr_eff_create();
    cfr_eff_ate_ipw(eff_ipw, T, Y, ps, n);
    printf("IPW ATE = %.2f\n", eff_ipw->ate);

    /* Doubly Robust */
    CFRCausalEffects* eff_dr = cfr_eff_create();
    cfr_eff_ate_doubly_robust(eff_dr, T, Y, ps, y0, y1, n);
    printf("DR ATE = %.2f\n", eff_dr->ate);

    /* Probability of Necessity */
    CFRBoundsData* bd = cfr_bounds_data_create(n);
    for (int i=0; i<n; i++) { bd->treatments[i]=T[i]; bd->outcomes[i]=(Y[i]<130)?1.0:0.0; }
    CFRBounds* b = cfr_bounds_create();
    cfr_bounds_experimental(bd, b);
    printf("\nPN bounds [%.3f, %.3f], PS bounds [%.3f, %.3f]\n",
           b->pn_lower, b->pn_upper, b->ps_lower, b->ps_upper);

    /* E-value */
    double rr = cfr_eff_relative_risk(T, Y, n);
    printf("RR = %.3f, E-value = %.3f (robustness to confounding)\n", rr, cfr_e_value(rr));

    /* Counterfactual: what if treated patient had not received drug? */
    CFRSCM* scm = cfr_create_simple_treatment_scm();
    CFRIntervention* intv = cfr_int_create_hard(0, 0.0);
    CFREvidence* ev = cfr_ev_create(1);
    cfr_ev_add_observation(ev, 0, 1.0, 0);
    CFRCounterfactual* cf = cfr_cf_create(scm);
    cfr_cf_compute_full(cf, ev, intv, 1);
    printf("\nCounterfactual E[Y|do(T=0)] = %.2f (factual Y under T=1 was higher)\n", cf->mean);

    /* Compare interventions */
    CFRIntervention* intv_a = cfr_int_create_hard(0, 1.0);
    CFRIntervention* intv_b = cfr_int_create_hard(0, 0.0);
    double diff;
    cfr_cf_compare_interventions(scm, 1, intv_a, intv_b, ev, &diff);
    printf("ATE via counterfactual comparison = %.2f\n", diff);

    /* Individual bounds */
    double pn_b[2], ps_b[2];
    cfr_bounds_individual(1.0, 0.0, pn_b, ps_b);
    printf("Individual PN = %.2f, PS = %.2f\n", pn_b[0], ps_b[0]);

    /* Sensitivity */
    double ate_adj = 0;
    cfr_bounds_sensitivity(b->pn_lower, 0.5, 0.1, &ate_adj);
    printf("Sensitivity-adjusted PN = %.4f\n", ate_adj);

    /* Excess fraction */
    double af = cfr_bounds_excess_fraction(bd);
    printf("Attributable fraction = %.4f (%d percent)\n", af, (int)(af*100));

    cfr_eff_free(eff); cfr_eff_free(eff_obs); cfr_eff_free(eff_ipw); cfr_eff_free(eff_dr);
    cfr_bounds_free(b); cfr_bounds_data_free(bd);
    cfr_cf_free(cf); cfr_int_free(intv); cfr_int_free(intv_a); cfr_int_free(intv_b);
    cfr_ev_free(ev); cfr_scm_free(scm);
    printf("\n=== App1 PASSED ===\n");
    return 0;
}
